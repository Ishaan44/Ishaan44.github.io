import argparse
import pickle
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
from scipy.spatial import cKDTree


CLUSTER_COLORS = ["#0f766e", "#c2410c", "#4338ca", "#a16207"]


class UnionFind:
    def __init__(self, size):
        self.parent = list(range(size))
        self.rank = [0] * size

    def find(self, item):
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, left, right):
        root_left = self.find(left)
        root_right = self.find(right)

        if root_left == root_right:
            return

        if self.rank[root_left] < self.rank[root_right]:
            self.parent[root_left] = root_right
        elif self.rank[root_left] > self.rank[root_right]:
            self.parent[root_right] = root_left
        else:
            self.parent[root_right] = root_left
            self.rank[root_left] += 1


def load_payload(path):
    with Path(path).open("rb") as f:
        payload = pickle.load(f)

    if isinstance(payload, dict) and "loopsoup" in payload:
        return payload

    return {
        "N": None,
        "c": None,
        "seed": None,
        "dimension": 3,
        "loopsoup": payload,
        "stats": {},
    }


def infer_N(loops):
    maximum = 0
    for loop in loops:
        for x, y, z in loop:
            maximum = max(maximum, x, y, z)
    return maximum + 1


def loop_measurements(loop):
    points = np.asarray(loop, dtype=float)
    length = len(loop) - 1
    span = points.max(axis=0) - points.min(axis=0)
    diameter = float(np.linalg.norm(span))
    return length, diameter


def filtered_records(loops, min_length, min_diameter):
    records = []

    for source_index, loop in enumerate(loops):
        length, diameter = loop_measurements(loop)
        if length < min_length or diameter < min_diameter:
            continue

        records.append(
            {
                "source_index": source_index,
                "loop": loop,
                "length": length,
                "diameter": diameter,
            }
        )

    return records


def build_intersection_clusters(records):
    union_find = UnionFind(len(records))
    first_loop_at_vertex = {}

    for loop_index, record in enumerate(records):
        for point in record["loop"]:
            vertex = tuple(point)
            previous_loop = first_loop_at_vertex.get(vertex)

            if previous_loop is None:
                first_loop_at_vertex[vertex] = loop_index
            else:
                union_find.union(loop_index, previous_loop)

    clusters = {}
    for loop_index in range(len(records)):
        root = union_find.find(loop_index)
        clusters.setdefault(root, []).append(loop_index)

    return list(clusters.values())


def cluster_stats(records, cluster_indices):
    mins = np.array([np.inf, np.inf, np.inf])
    maxs = np.array([-np.inf, -np.inf, -np.inf])
    total_steps = 0
    max_loop_length = 0

    for loop_index in cluster_indices:
        record = records[loop_index]
        points = np.asarray(record["loop"], dtype=float)
        mins = np.minimum(mins, points.min(axis=0))
        maxs = np.maximum(maxs, points.max(axis=0))
        total_steps += record["length"]
        max_loop_length = max(max_loop_length, record["length"])

    diameter = float(np.linalg.norm(maxs - mins))

    return {
        "loop_count": len(cluster_indices),
        "total_steps": total_steps,
        "diameter": diameter,
        "max_loop_length": max_loop_length,
        "bbox_min": mins,
        "bbox_max": maxs,
    }


def annotated_clusters(records, clusters):
    annotated = []

    for cluster_indices in clusters:
        stats = cluster_stats(records, cluster_indices)
        annotated.append({"indices": cluster_indices, "stats": stats})

    return annotated


def eligible_clusters(annotated, min_cluster_loops):
    eligible = [
        cluster for cluster in annotated if cluster["stats"]["loop_count"] >= min_cluster_loops
    ]
    if len(eligible) < 2:
        return annotated

    return eligible


def ranked_clusters(annotated, rank_by):
    if rank_by == "closest":
        raise ValueError("Use closest_cluster_pair for closest-pair selection.")

    eligible = annotated

    def score(cluster):
        stats = cluster["stats"]
        if rank_by == "loops":
            return stats["loop_count"]
        if rank_by == "diameter":
            return stats["diameter"]
        if rank_by == "max-loop":
            return stats["max_loop_length"]
        return stats["total_steps"]

    return sorted(eligible, key=score, reverse=True)


def cluster_vertices(records, cluster_indices):
    vertices = set()

    for loop_index in cluster_indices:
        for point in records[loop_index]["loop"]:
            vertices.add(tuple(point))

    return np.array(list(vertices), dtype=float)


def closest_cluster_pair(records, annotated):
    if len(annotated) < 2:
        return annotated, None

    point_blocks = []
    owners = []

    for cluster_index, cluster in enumerate(annotated):
        vertices = cluster_vertices(records, cluster["indices"])
        cluster["vertices"] = vertices
        point_blocks.append(vertices)
        owners.extend([cluster_index] * len(vertices))

    points = np.vstack(point_blocks)
    owners = np.asarray(owners, dtype=int)
    tree = cKDTree(points)
    point_count = len(points)
    k = min(8, point_count)
    best = None
    best_detail = None

    while True:
        distances, neighbors = tree.query(points, k=k)
        if k == 1:
            distances = distances[:, None]
            neighbors = neighbors[:, None]

        unresolved_bounds = []

        for point_index in range(point_count):
            owner = owners[point_index]
            found_other_cluster = False

            for neighbor_distance, neighbor_index in zip(
                distances[point_index, 1:],
                neighbors[point_index, 1:],
            ):
                neighbor_owner = owners[neighbor_index]
                if neighbor_owner == owner:
                    continue

                left, right = sorted((int(owner), int(neighbor_owner)))
                combined_steps = (
                    annotated[left]["stats"]["total_steps"]
                    + annotated[right]["stats"]["total_steps"]
                )
                combined_diameter = (
                    annotated[left]["stats"]["diameter"]
                    + annotated[right]["stats"]["diameter"]
                )
                candidate = (
                    float(neighbor_distance),
                    -combined_steps,
                    -combined_diameter,
                    left,
                    right,
                )

                if best is None or candidate < best:
                    best = candidate
                    best_detail = (
                        int(owner),
                        int(neighbor_owner),
                        points[point_index].copy(),
                        points[neighbor_index].copy(),
                    )

                found_other_cluster = True
                break

            if not found_other_cluster:
                unresolved_bounds.append(float(distances[point_index, -1]))

        if best is None:
            if k == point_count:
                break
            k = min(point_count, k * 2)
            continue

        if not unresolved_bounds or best[0] <= min(unresolved_bounds) or k == point_count:
            break

        k = min(point_count, k * 2)

    if best is None:
        return annotated[:2], None

    distance, _, _, left, right = best
    selected = [annotated[left], annotated[right]]
    owner_a, owner_b, point_a, point_b = best_detail

    if owner_a == left:
        annotated[left]["closest_point"] = point_a
        annotated[right]["closest_point"] = point_b
    else:
        annotated[left]["closest_point"] = point_b
        annotated[right]["closest_point"] = point_a

    for cluster in selected:
        cluster["selection_distance"] = distance

    return selected, distance


def loops_to_render(records, cluster_indices, max_points):
    ordered = sorted(cluster_indices, key=lambda index: records[index]["length"], reverse=True)
    selected = []
    point_count = 0

    for loop_index in ordered:
        loop = records[loop_index]["loop"]
        if max_points is not None and point_count + len(loop) > max_points:
            if selected:
                continue
            selected.append(loop_index)
            point_count += len(loop)
            break
        selected.append(loop_index)
        point_count += len(loop)

    return selected, point_count


def cube_trace():
    vertices = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ],
        dtype=float,
    )
    edges = [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 0),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 4),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    ]

    x, y, z = [], [], []
    for start, end in edges:
        x.extend([vertices[start, 0], vertices[end, 0], None])
        y.extend([vertices[start, 1], vertices[end, 1], None])
        z.extend([vertices[start, 2], vertices[end, 2], None])

    return go.Scatter3d(
        x=x,
        y=y,
        z=z,
        mode="lines",
        line=dict(color="rgba(24, 24, 27, 0.42)", width=3),
        hoverinfo="skip",
        showlegend=False,
    )


def cluster_trace(records, loop_indices, N, color, name, width):
    scale = max(1, N - 1)
    x, y, z = [], [], []

    for loop_index in loop_indices:
        for px, py, pz in records[loop_index]["loop"]:
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
        opacity=0.78,
        name=name,
        hoverinfo="skip",
    )


def bbox_trace(stats, N, color, name):
    scale = max(1, N - 1)
    mins = stats["bbox_min"] / scale
    maxs = stats["bbox_max"] / scale

    vertices = np.array(
        [
            [mins[0], mins[1], mins[2]],
            [maxs[0], mins[1], mins[2]],
            [maxs[0], maxs[1], mins[2]],
            [mins[0], maxs[1], mins[2]],
            [mins[0], mins[1], maxs[2]],
            [maxs[0], mins[1], maxs[2]],
            [maxs[0], maxs[1], maxs[2]],
            [mins[0], maxs[1], maxs[2]],
        ]
    )
    edges = [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 0),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 4),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    ]

    x, y, z = [], [], []
    for start, end in edges:
        x.extend([vertices[start, 0], vertices[end, 0], None])
        y.extend([vertices[start, 1], vertices[end, 1], None])
        z.extend([vertices[start, 2], vertices[end, 2], None])

    return go.Scatter3d(
        x=x,
        y=y,
        z=z,
        mode="lines",
        line=dict(color=color, width=2, dash="dash"),
        opacity=0.32,
        name=f"{name} bbox",
        hoverinfo="skip",
        showlegend=False,
    )


def closest_distance_trace(selected_clusters, N):
    if len(selected_clusters) < 2:
        return None

    point_a = selected_clusters[0].get("closest_point")
    point_b = selected_clusters[1].get("closest_point")
    if point_a is None or point_b is None:
        return None

    scale = max(1, N - 1)
    point_a = np.asarray(point_a, dtype=float) / scale
    point_b = np.asarray(point_b, dtype=float) / scale

    return go.Scatter3d(
        x=[point_a[0], point_b[0]],
        y=[point_a[1], point_b[1]],
        z=[point_a[2], point_b[2]],
        mode="lines+markers",
        line=dict(color="#111827", width=9),
        marker=dict(size=5, color="#111827"),
        name="shortest distance",
        hovertemplate="shortest distance<extra></extra>",
    )


def make_figure(payload, args):
    loops = payload["loopsoup"]
    N = payload.get("N") or infer_N(loops)
    records = filtered_records(loops, args.min_length, args.min_diameter)

    if not records:
        raise ValueError("No loops survived the length/diameter filters.")

    clusters = build_intersection_clusters(records)
    annotated = eligible_clusters(annotated_clusters(records, clusters), args.min_cluster_loops)
    selection_distance = None

    if args.rank_by == "closest":
        selected_clusters, selection_distance = closest_cluster_pair(records, annotated)
    else:
        ranked = ranked_clusters(annotated, args.rank_by)
        selected_clusters = ranked[: args.clusters]

    if len(selected_clusters) < args.clusters:
        raise ValueError(f"Only found {len(selected_clusters)} cluster(s) after filtering.")

    traces = [cube_trace()]
    summary_lines = []

    distance_trace = closest_distance_trace(selected_clusters, N)
    if distance_trace is not None:
        traces.append(distance_trace)

    for display_index, cluster in enumerate(selected_clusters, start=1):
        color = CLUSTER_COLORS[(display_index - 1) % len(CLUSTER_COLORS)]
        render_indices, render_points = loops_to_render(
            records,
            cluster["indices"],
            args.max_points_per_cluster,
        )
        stats = cluster["stats"]
        name = (
            f"cluster {display_index}: "
            f"{stats['loop_count']:,} loops, {stats['total_steps']:,} steps"
        )

        traces.append(
            cluster_trace(
                records,
                render_indices,
                N,
                color,
                name,
                width=max(3, args.line_width - display_index + 1),
            )
        )

        if args.show_boxes:
            traces.append(bbox_trace(stats, N, color, f"cluster {display_index}"))

        summary_lines.append(
            f"Cluster {display_index}: {stats['loop_count']:,} loops, "
            f"{stats['total_steps']:,} total steps, "
            f"diameter {stats['diameter'] / max(1, N - 1):.3f}, "
            f"rendered {len(render_indices):,} loops / {render_points:,} points"
        )

    if selection_distance is not None:
        summary_lines.append(
            f"Closest distance: {selection_distance / max(1, N - 1):.4f} rescaled "
            f"({selection_distance:.2f} lattice units)"
        )

    title_prefix = (
        "Two Closest Intersection Clusters"
        if args.rank_by == "closest"
        else "Two Intersection Clusters"
    )
    title = (
        f"{title_prefix} in 3D RWLS "
        f"(N={N}, c={payload.get('c')}, filtered loops={len(records):,})"
    )

    fig = go.Figure(data=traces)
    fig.update_layout(
        title=dict(text=title, x=0.02, y=0.97, xanchor="left", font=dict(size=20)),
        paper_bgcolor="#f7f5ef",
        plot_bgcolor="#f7f5ef",
        margin=dict(l=0, r=0, t=58, b=0),
        height=args.height,
        legend=dict(
            x=0.02,
            y=0.88,
            bgcolor="rgba(255,255,255,0.82)",
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
                text=(
                    "<br>".join(summary_lines)
                    + f"<br>Intersection rule: loops share a lattice vertex. "
                    + f"Filters: length >= {args.min_length}, diameter >= {args.min_diameter}."
                ),
                x=0.02,
                y=0.02,
                xref="paper",
                yref="paper",
                showarrow=False,
                align="left",
                bgcolor="rgba(255,255,255,0.82)",
                bordercolor="rgba(24,24,27,0.16)",
                borderwidth=1,
                borderpad=8,
                font=dict(size=13, color="#374151"),
            )
        ],
    )

    return fig, len(records), len(clusters), selected_clusters


def parse_args():
    parser = argparse.ArgumentParser(
        description="Render only a small number of loop-intersection clusters from a 3D RWLS pickle."
    )
    parser.add_argument("pickle_path", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--clusters", type=int, default=2)
    parser.add_argument("--min-length", type=int, default=8)
    parser.add_argument("--min-diameter", type=float, default=0.0)
    parser.add_argument("--min-cluster-loops", type=int, default=1)
    parser.add_argument(
        "--rank-by",
        choices=["closest", "steps", "loops", "diameter", "max-loop"],
        default="closest",
    )
    parser.add_argument("--max-points-per-cluster", type=int, default=180000)
    parser.add_argument("--line-width", type=float, default=6.0)
    parser.add_argument("--height", type=int, default=900)
    parser.add_argument("--show-boxes", action="store_true")
    parser.add_argument(
        "--cdn",
        action="store_true",
        help="Use Plotly from CDN instead of embedding it. Smaller file, needs internet in browser.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    payload = load_payload(args.pickle_path)
    output = args.output or args.pickle_path.with_name(args.pickle_path.stem + "-two-clusters.html")

    fig, filtered_count, cluster_count, selected_clusters = make_figure(payload, args)
    output.parent.mkdir(parents=True, exist_ok=True)
    pio.write_html(
        fig,
        file=output,
        include_plotlyjs="cdn" if args.cdn else True,
        full_html=True,
        config={
            "responsive": True,
            "displaylogo": False,
            "scrollZoom": True,
            "modeBarButtonsToRemove": ["lasso2d", "select2d"],
        },
    )

    print(f"Filtered loops clustered: {filtered_count}")
    print(f"Total clusters found: {cluster_count}")
    for index, cluster in enumerate(selected_clusters, start=1):
        stats = cluster["stats"]
        print(
            f"Cluster {index}: loops={stats['loop_count']}, "
            f"steps={stats['total_steps']}, diameter={stats['diameter']:.2f}, "
            f"max_loop={stats['max_loop_length']}"
        )
    closest_distance = selected_clusters[0].get("selection_distance")
    if closest_distance is not None:
        print(f"Closest distance: {closest_distance:.2f} lattice units")
    print(f"Saved two-cluster viewer to {output}")


if __name__ == "__main__":
    main()
