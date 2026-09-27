import argparse
from collections import deque
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

from RWLS_3D import loop_soup_3d, loop_stats
from view_rwls3d_two_clusters import (
    CLUSTER_COLORS,
    UnionFind,
    annotated_clusters,
    bbox_trace,
    build_intersection_clusters,
    closest_cluster_pair,
    closest_distance_trace,
    cluster_trace,
    cluster_vertices,
    cube_trace,
    eligible_clusters,
    infer_N,
    load_payload,
    loop_measurements,
)


DUST_COLOR = "#7c3aed"
CONTEXT_DUST_COLOR = "rgba(79, 70, 229, 0.34)"


def cluster_bbox(macro_records, selected_clusters):
    mins = np.array([np.inf, np.inf, np.inf])
    maxs = np.array([-np.inf, -np.inf, -np.inf])

    for cluster in selected_clusters:
        mins = np.minimum(mins, cluster["stats"]["bbox_min"])
        maxs = np.maximum(maxs, cluster["stats"]["bbox_max"])

    return mins.astype(int), maxs.astype(int)


def local_box_from_clusters(macro_records, selected_clusters, N, padding):
    mins, maxs = cluster_bbox(macro_records, selected_clusters)
    lower = mins - padding
    upper = maxs + padding
    side = int(max(upper - lower + 1))
    side = max(3, min(side, N))

    center = (lower + upper) / 2.0
    lower = np.floor(center - (side - 1) / 2.0).astype(int)
    lower = np.maximum(0, np.minimum(lower, N - side))
    upper = lower + side - 1
    side_lengths = np.array([side, side, side])

    return tuple(lower.tolist()), tuple(upper.tolist()), tuple(side_lengths.tolist())


def translate_loops(loops, offset):
    ox, oy, oz = offset
    return [
        [(x + ox, y + oy, z + oz) for x, y, z in loop]
        for loop in loops
    ]


def loop_intersects_box(loop, lower, upper):
    lx, ly, lz = lower
    ux, uy, uz = upper

    for x, y, z in loop:
        if lx <= x <= ux and ly <= y <= uy and lz <= z <= uz:
            return True

    return False


def generate_local_dust_soup(box_shape, c, seed, progress_every):
    sx, sy, sz = box_shape
    if sx != sy or sy != sz:
        raise ValueError(f"Local dust domain must be cubic, got shape={box_shape}.")

    side = sx

    if side < 3:
        return []

    return loop_soup_3d(side, c, seed=seed, progress_every=progress_every)


def clipped_translated_local_dust(box_shape, offset, c, seed, progress_every):
    loops = generate_local_dust_soup(box_shape, c, seed, progress_every)
    translated = translate_loops(loops, offset)
    upper = tuple(offset[i] + box_shape[i] - 1 for i in range(3))
    return [
        loop
        for loop in translated
        if loop_intersects_box(loop, offset, upper)
    ]


def scale_value(length, diameter, cutoff_metric):
    if cutoff_metric == "diameter":
        return float(diameter)
    return float(length)


def scale_label(cutoff_metric):
    if cutoff_metric == "diameter":
        return "lattice diameter"
    return "time length"


def make_records(loops, min_size=None, max_size=None, cutoff_metric="length"):
    records = []

    for source_index, loop in enumerate(loops):
        length, diameter = loop_measurements(loop)
        size = scale_value(length, diameter, cutoff_metric)

        if min_size is not None and size < min_size:
            continue
        if max_size is not None and size >= max_size:
            continue

        records.append(
            {
                "source_index": source_index,
                "loop": loop,
                "length": length,
                "diameter": diameter,
                "size": size,
            }
        )

    return records


def load_or_generate_soup(path, N, c, seed, progress_every, label):
    if path is not None:
        payload = load_payload(path)
        return payload["loopsoup"], payload.get("N"), payload.get("c"), payload.get("seed")

    print(f"Generating independent {label} RWLS: N={N}, c={c}, seed={seed}")
    loopsoup = loop_soup_3d(N, c, seed=seed, progress_every=progress_every)
    return loopsoup, N, c, seed


def choose_closest_macro_clusters(macro_records, min_macro_cluster_loops):
    clusters = build_intersection_clusters(macro_records)
    annotated = annotated_clusters(macro_records, clusters)
    eligible = [
        cluster
        for cluster in annotated
        if cluster["stats"]["loop_count"] >= min_macro_cluster_loops
    ]
    if len(eligible) < 2:
        raise ValueError(
            f"Only found {len(eligible)} macro cluster(s) with at least "
            f"{min_macro_cluster_loops} loop(s). Lower Min Macro Cluster Loops or Macro Cutoff L."
        )

    selected, distance = closest_cluster_pair(macro_records, eligible)

    if len(selected) < 2:
        raise ValueError("Could not find two macro clusters. Lower the macro cutoff.")

    return selected, distance, len(clusters), len(eligible)


def unique_loop_vertices(loop):
    return {tuple(point) for point in loop}


def add_edge(adjacency, left, right):
    adjacency.setdefault(left, set()).add(right)
    adjacency.setdefault(right, set()).add(left)


def shortest_graph_path(adjacency, start, target):
    queue = deque([start])
    previous = {start: None}

    while queue:
        node = queue.popleft()
        if node == target:
            break

        for neighbor in adjacency.get(node, ()):
            if neighbor in previous:
                continue
            previous[neighbor] = node
            queue.append(neighbor)

    if target not in previous:
        return []

    path = []
    node = target
    while node is not None:
        path.append(node)
        node = previous[node]

    return path[::-1]


def sprinkle_until_glued(dust_records, macro_records, cluster_a, cluster_b, epsilon_min):
    dust_candidates = [
        record for record in dust_records if record["size"] >= epsilon_min
    ]
    dust_candidates.sort(key=lambda record: (record["size"], record["length"]), reverse=True)

    a_node = 0
    b_node = 1
    union_find = UnionFind(2 + len(dust_candidates))
    adjacency = {a_node: set(), b_node: set()}
    vertex_owner = {}

    for vertex in map(tuple, cluster_vertices(macro_records, cluster_a["indices"])):
        vertex_owner[vertex] = a_node

    for vertex in map(tuple, cluster_vertices(macro_records, cluster_b["indices"])):
        existing = vertex_owner.get(vertex)
        if existing is None:
            vertex_owner[vertex] = b_node
        else:
            union_find.union(existing, b_node)
            add_edge(adjacency, existing, b_node)

    activated_nodes = []
    node_to_dust_index = {}
    gluing_node = None
    gluing_cutoff = None

    for dust_index, record in enumerate(dust_candidates):
        node = 2 + dust_index
        node_to_dust_index[node] = dust_index
        activated_nodes.append(node)
        adjacency.setdefault(node, set())
        touched_owners = set()

        for vertex in unique_loop_vertices(record["loop"]):
            owner = vertex_owner.get(vertex)
            if owner is None:
                vertex_owner[vertex] = node
            else:
                touched_owners.add(owner)

        for owner in touched_owners:
            union_find.union(node, owner)
            add_edge(adjacency, node, owner)

        if union_find.find(a_node) == union_find.find(b_node):
            gluing_node = node
            gluing_cutoff = record["size"]
            break

    path = shortest_graph_path(adjacency, a_node, b_node) if gluing_node is not None else []
    bridge_nodes = [node for node in path if node >= 2]
    bridge_dust_indices = [node_to_dust_index[node] for node in bridge_nodes]
    activated_dust_indices = [node_to_dust_index[node] for node in activated_nodes]

    return {
        "glued": gluing_node is not None,
        "gluing_cutoff": gluing_cutoff,
        "dust_candidates": dust_candidates,
        "bridge_dust_indices": bridge_dust_indices,
        "activated_dust_indices": activated_dust_indices,
        "bridge_path": path,
    }


def dust_trace(records, indices, N, name, color, width, opacity=0.82):
    if not indices:
        return None

    scale = max(1, N - 1)
    x, y, z = [], [], []

    for index in indices:
        for px, py, pz in records[index]["loop"]:
            x.append(px / scale)
            y.append(py / scale)
            z.append(pz / scale)
        x.append(None)
        y.append(None)
        z.append(None)

    return go.Scatter3d(
        x=x,
        y=y,
        z=z,
        mode="lines",
        line=dict(color=color, width=width),
        opacity=opacity,
        name=name,
        hoverinfo="skip",
    )


def capped_context_indices(records, indices, max_points):
    selected = []
    points = 0

    for index in indices:
        loop_points = len(records[index]["loop"])
        if points + loop_points > max_points:
            continue
        selected.append(index)
        points += loop_points

    return selected


def build_figure(
    payload_meta,
    macro_records,
    selected_clusters,
    sprinkle_result,
    N,
    macro_cutoff,
    epsilon_min,
    cutoff_metric,
    show_boxes,
    show_context_dust,
    max_context_points,
    height,
):
    cluster_a, cluster_b = selected_clusters
    dust_records = sprinkle_result["dust_candidates"]
    bridge_indices = sprinkle_result["bridge_dust_indices"]
    activated_indices = sprinkle_result["activated_dust_indices"]

    traces = [cube_trace()]

    if show_context_dust and activated_indices:
        context_indices = capped_context_indices(dust_records, activated_indices, max_context_points)
        context = dust_trace(
            dust_records,
            context_indices,
            N,
            f"sprinkled dust shown: {len(context_indices):,} loops",
            CONTEXT_DUST_COLOR,
            width=3,
            opacity=0.48,
        )
        if context is not None:
            traces.append(context)

    traces.append(
        cluster_trace(
            macro_records,
            cluster_a["indices"],
            N,
            CLUSTER_COLORS[0],
            f"macro cluster A: {cluster_a['stats']['loop_count']:,} loops",
            width=6,
        )
    )
    traces.append(
        cluster_trace(
            macro_records,
            cluster_b["indices"],
            N,
            CLUSTER_COLORS[1],
            f"macro cluster B: {cluster_b['stats']['loop_count']:,} loops",
            width=6,
        )
    )

    bridge = dust_trace(
        dust_records,
        bridge_indices,
        N,
        f"gluing dust bridge: {len(bridge_indices):,} loops",
        DUST_COLOR,
        width=7,
        opacity=0.90,
    )
    if bridge is not None:
        traces.append(bridge)

    distance = closest_distance_trace(selected_clusters, N)
    if distance is not None:
        traces.append(distance)

    if show_boxes:
        traces.append(bbox_trace(cluster_a["stats"], N, CLUSTER_COLORS[0], "macro A"))
        traces.append(bbox_trace(cluster_b["stats"], N, CLUSTER_COLORS[1], "macro B"))

    gluing_cutoff = sprinkle_result["gluing_cutoff"]
    if sprinkle_result["glued"]:
        if cutoff_metric == "diameter":
            gluing_text = (
                f"CONNECTED: glues when epsilon reaches "
                f"{gluing_cutoff / max(1, N - 1):.4f} rescaled "
                f"({gluing_cutoff:.2f} lattice diameter)"
            )
        else:
            gluing_text = (
                f"CONNECTED: glues when epsilon reaches "
                f"{gluing_cutoff:.0f} time steps"
            )
    else:
        gluing_text = f"NOT CONNECTED for epsilon >= {epsilon_min:.2f}"

    summary = [
        f"Cutoff metric: {scale_label(cutoff_metric)}",
        f"Macro cutoff L: {macro_cutoff:.2f}",
        f"Dust window: epsilon <= {scale_label(cutoff_metric)} < L, epsilon_min={epsilon_min:.2f}",
        gluing_text,
        f"Bridge dust loops: {len(bridge_indices):,}",
        f"Dust candidates inspected: {len(activated_indices):,} / {len(dust_records):,}",
        "Macro clusters intersect by shared lattice vertices; dust loops are independent.",
    ]

    title = (
        f"Loop-Soup Sprinkling Gluing "
        f"(N={N}, c={payload_meta['macro_c']}, dust c={payload_meta['dust_c']})"
    )

    fig = go.Figure(data=traces)
    fig.update_layout(
        title=dict(text=title, x=0.02, y=0.97, xanchor="left", font=dict(size=20)),
        paper_bgcolor="#f7f5ef",
        plot_bgcolor="#f7f5ef",
        margin=dict(l=0, r=0, t=58, b=0),
        height=height,
        legend=dict(
            x=0.02,
            y=0.88,
            bgcolor="rgba(255,255,255,0.84)",
            bordercolor="rgba(24,24,27,0.16)",
            borderwidth=1,
            font=dict(size=12),
        ),
        scene=dict(
            xaxis=dict(title="x", range=[0, 1], gridcolor="rgba(24,24,27,0.14)"),
            yaxis=dict(title="y", range=[0, 1], gridcolor="rgba(24,24,27,0.14)"),
            zaxis=dict(title="z", range=[0, 1], gridcolor="rgba(24,24,27,0.14)"),
            aspectmode="cube",
            camera=dict(eye=dict(x=1.55, y=1.42, z=1.16), up=dict(x=0, y=0, z=1)),
        ),
        annotations=[
            dict(
                text="<br>".join(summary),
                x=0.02,
                y=0.02,
                xref="paper",
                yref="paper",
                showarrow=False,
                align="left",
                bgcolor="rgba(255,255,255,0.84)",
                bordercolor="rgba(24,24,27,0.16)",
                borderwidth=1,
                borderpad=8,
                font=dict(size=13, color="#374151"),
            )
        ],
    )

    return fig


def run_experiment(args):
    dust_c = args.dust_c if args.dust_c is not None else args.c
    macro_loops, macro_N, macro_c, macro_seed = load_or_generate_soup(
        args.macro_pickle,
        args.N,
        args.c,
        args.macro_seed,
        args.progress_every,
        "macro",
    )

    N = macro_N or args.N
    macro_c = macro_c if macro_c is not None else args.c

    macro_records = make_records(
        macro_loops,
        min_size=args.macro_cutoff,
        cutoff_metric=args.cutoff_metric,
    )

    if len(macro_records) < 2:
        raise ValueError("Not enough macro loops. Lower --macro-cutoff.")

    selected_clusters, closest_distance, cluster_count, eligible_count = choose_closest_macro_clusters(
        macro_records,
        args.min_macro_cluster_loops,
    )

    local_box = None

    if args.same_soup_dust:
        dust_loops = macro_loops
        dust_N = N
        dust_seed = macro_seed
        loaded_dust_c = macro_c
    elif args.local_dust:
        lower, upper, box_shape = local_box_from_clusters(
            macro_records,
            selected_clusters,
            int(N),
            args.local_padding,
        )
        print(
            f"Generating local dust RWLS in box lower={lower}, upper={upper}, "
            f"shape={box_shape}, c={dust_c}, seed={args.dust_seed}"
        )
        dust_loops = clipped_translated_local_dust(
            box_shape,
            lower,
            dust_c,
            args.dust_seed,
            args.progress_every,
        )
        dust_N = N
        dust_seed = args.dust_seed
        loaded_dust_c = dust_c
        local_box = {"lower": lower, "upper": upper, "shape": box_shape}
    else:
        dust_loops, dust_N, loaded_dust_c, dust_seed = load_or_generate_soup(
            args.dust_pickle,
            N,
            dust_c,
            args.dust_seed,
            args.progress_every,
            "dust",
        )

    if dust_N is not None and int(dust_N) != int(N):
        raise ValueError(f"Macro N={N} and dust N={dust_N} do not match.")

    dust_records = make_records(
        dust_loops,
        min_size=args.epsilon_min,
        max_size=args.macro_cutoff,
        cutoff_metric=args.cutoff_metric,
    )

    if not dust_records:
        raise ValueError("No dust loops in the requested scale window.")

    sprinkle_result = sprinkle_until_glued(
        dust_records,
        macro_records,
        selected_clusters[0],
        selected_clusters[1],
        args.epsilon_min,
    )

    payload_meta = {
        "macro_c": macro_c,
        "dust_c": loaded_dust_c if loaded_dust_c is not None else dust_c,
        "macro_seed": macro_seed,
        "dust_seed": dust_seed,
    }

    fig = build_figure(
        payload_meta,
        macro_records,
        selected_clusters,
        sprinkle_result,
        int(N),
        args.macro_cutoff,
        args.epsilon_min,
        args.cutoff_metric,
        args.show_boxes,
        args.show_context_dust,
        args.max_context_points,
        args.height,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    pio.write_html(
        fig,
        file=args.output,
        include_plotlyjs="cdn" if args.cdn else True,
        full_html=True,
        config={
            "responsive": True,
            "displaylogo": False,
            "scrollZoom": True,
            "modeBarButtonsToRemove": ["lasso2d", "select2d"],
        },
    )

    stats = {
        "N": int(N),
        "macro_loop_count": len(macro_records),
        "macro_cluster_count": cluster_count,
        "eligible_macro_cluster_count": eligible_count,
        "closest_macro_distance": closest_distance,
        "dust_loop_count": len(dust_records),
        "glued": sprinkle_result["glued"],
        "gluing_cutoff": sprinkle_result["gluing_cutoff"],
        "bridge_dust_loop_count": len(sprinkle_result["bridge_dust_indices"]),
        "activated_dust_loop_count": len(sprinkle_result["activated_dust_indices"]),
        "cutoff_metric": args.cutoff_metric,
        "local_dust": args.local_dust,
        "local_box": local_box,
        "macro_soup_stats": loop_stats(macro_loops),
        "dust_soup_stats": loop_stats(dust_loops),
    }

    return stats


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Sprinkle an independent small-loop soup onto two closest macro loop clusters "
            "and visualize the dust scale at which they glue."
        )
    )
    parser.add_argument("--macro-pickle", type=Path, default=None)
    parser.add_argument("--dust-pickle", type=Path, default=None)
    parser.add_argument("--same-soup-dust", action="store_true")
    parser.add_argument("--N", type=int, default=80)
    parser.add_argument("--c", type=float, default=0.5)
    parser.add_argument("--dust-c", type=float, default=None)
    parser.add_argument("--macro-seed", type=int, default=11)
    parser.add_argument("--dust-seed", type=int, default=29)
    parser.add_argument("--local-dust", action="store_true")
    parser.add_argument("--local-padding", type=int, default=8)
    parser.add_argument("--macro-cutoff", type=float, default=8.0)
    parser.add_argument("--epsilon-min", type=float, default=2.0)
    parser.add_argument(
        "--cutoff-metric",
        choices=["length", "diameter"],
        default="length",
        help="Use loop time length or spatial diameter for macro/dust cutoffs.",
    )
    parser.add_argument("--min-macro-cluster-loops", type=int, default=1)
    parser.add_argument("--progress-every", type=int, default=20)
    parser.add_argument("--show-boxes", action="store_true")
    parser.add_argument("--show-context-dust", action="store_true")
    parser.add_argument("--max-context-points", type=int, default=100000)
    parser.add_argument("--height", type=int, default=900)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("rwls3d/output/sprinkling-gluing.html"),
    )
    parser.add_argument(
        "--cdn",
        action="store_true",
        help="Use Plotly from CDN instead of embedding it. Smaller file, needs internet in browser.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    stats = run_experiment(args)

    print(f"N: {stats['N']}")
    print(f"Cutoff metric: {stats['cutoff_metric']}")
    print(f"Macro loops above cutoff: {stats['macro_loop_count']}")
    print(f"Macro clusters: {stats['macro_cluster_count']}")
    print(f"Dust loops in epsilon window: {stats['dust_loop_count']}")
    if stats["local_dust"]:
        print(f"Local dust box: {stats['local_box']}")
    print(f"Closest macro distance: {stats['closest_macro_distance']:.2f} lattice units")

    if stats["glued"]:
        print(f"Glued at epsilon*: {stats['gluing_cutoff']:.2f}")
        print(f"Bridge dust loops: {stats['bridge_dust_loop_count']}")
        print(f"Activated dust loops before gluing: {stats['activated_dust_loop_count']}")
    else:
        print("No gluing occurred inside the requested dust window.")

    print(f"Saved sprinkling viewer to {args.output}")


if __name__ == "__main__":
    main()
