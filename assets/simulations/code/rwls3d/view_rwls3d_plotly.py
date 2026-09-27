import argparse
import pickle
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
import plotly.io as pio


COLORS = [
    "#1f2937",
    "#2f6f68",
    "#2d8f78",
    "#79a33d",
    "#c98a1f",
    "#c05621",
    "#7c2d12",
]


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
        for point in loop:
            maximum = max(maximum, point[0], point[1], point[2])
    return maximum + 1


def loop_measurements(loop):
    length = len(loop) - 1
    points = np.asarray(loop, dtype=float)
    span = points.max(axis=0) - points.min(axis=0)
    diameter = float(np.linalg.norm(span))
    return length, diameter


def select_loops(loops, min_length, min_diameter, max_loops, max_points):
    candidates = []

    for loop in loops:
        length, diameter = loop_measurements(loop)
        if length < min_length or diameter < min_diameter:
            continue
        candidates.append((diameter, length, loop))

    candidates.sort(reverse=True, key=lambda item: (item[0], item[1]))

    selected = []
    total_points = 0
    for diameter, length, loop in candidates:
        if max_loops is not None and len(selected) >= max_loops:
            break
        if max_points is not None and total_points + len(loop) > max_points:
            break
        selected.append((diameter, length, loop))
        total_points += len(loop)

    return selected


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
    for a, b in edges:
        x.extend([vertices[a, 0], vertices[b, 0], None])
        y.extend([vertices[a, 1], vertices[b, 1], None])
        z.extend([vertices[a, 2], vertices[b, 2], None])

    return go.Scatter3d(
        x=x,
        y=y,
        z=z,
        mode="lines",
        line=dict(color="rgba(30, 30, 30, 0.42)", width=3),
        hoverinfo="skip",
        showlegend=False,
        name="cube",
    )


def bucket_index(value, edges):
    return int(np.searchsorted(edges, value, side="right"))


def loop_traces(selected, N, bucket_count):
    if not selected:
        return []

    diameters = np.array([item[0] for item in selected], dtype=float)
    max_diameter = max(float(diameters.max()), 1.0)
    edges = np.quantile(diameters, np.linspace(0, 1, bucket_count + 1)[1:-1])
    scale = max(1, N - 1)
    buckets = [[] for _ in range(bucket_count)]

    for diameter, length, loop in selected:
        buckets[bucket_index(diameter, edges)].append((diameter, length, loop))

    traces = []
    for index, bucket in enumerate(buckets):
        if not bucket:
            continue

        x, y, z = [], [], []
        lengths = []
        bucket_diameters = []

        for diameter, length, loop in bucket:
            bucket_diameters.append(diameter)
            lengths.append(length)
            for px, py, pz in loop:
                x.append(px / scale)
                y.append(py / scale)
                z.append(pz / scale)
            x.append(None)
            y.append(None)
            z.append(None)

        mean_diameter = float(np.mean(bucket_diameters))
        line_width = 2 + 5 * (mean_diameter / max_diameter)
        color = COLORS[min(index, len(COLORS) - 1)]
        name = (
            f"diam {min(bucket_diameters) / scale:.3f}-"
            f"{max(bucket_diameters) / scale:.3f}"
        )

        traces.append(
            go.Scatter3d(
                x=x,
                y=y,
                z=z,
                mode="lines",
                line=dict(color=color, width=line_width),
                opacity=0.70,
                name=name,
                hoverinfo="skip",
                legendgroup=f"bucket-{index}",
            )
        )

    return traces


def make_figure(payload, min_length, min_diameter, max_loops, max_points, bucket_count):
    loops = payload["loopsoup"]
    N = payload.get("N") or infer_N(loops)
    c = payload.get("c")
    selected = select_loops(loops, min_length, min_diameter, max_loops, max_points)

    source_count = len(loops)
    shown_count = len(selected)
    shown_steps = sum(item[1] for item in selected)
    largest = max((item[1] for item in selected), default=0)

    traces = [cube_trace()]
    traces.extend(loop_traces(selected, N, bucket_count))

    title = (
        f"3D Random Walk Loop Soup in a Cube "
        f"(N={N}, c={c}, shown={shown_count:,}/{source_count:,})"
    )

    fig = go.Figure(data=traces)
    fig.update_layout(
        title=dict(
            text=title,
            x=0.02,
            y=0.97,
            xanchor="left",
            font=dict(size=20, color="#1f2937"),
        ),
        paper_bgcolor="#f7f5ef",
        plot_bgcolor="#f7f5ef",
        margin=dict(l=0, r=0, t=58, b=0),
        width=None,
        height=900,
        legend=dict(
            title="Loop diameter, rescaled",
            x=0.02,
            y=0.88,
            bgcolor="rgba(255,255,255,0.78)",
            bordercolor="rgba(31,41,55,0.16)",
            borderwidth=1,
            font=dict(size=12),
        ),
        scene=dict(
            xaxis=dict(
                title="x",
                range=[0, 1],
                showbackground=True,
                backgroundcolor="rgba(255,255,255,0.54)",
                gridcolor="rgba(31,41,55,0.14)",
                zerolinecolor="rgba(31,41,55,0.18)",
            ),
            yaxis=dict(
                title="y",
                range=[0, 1],
                showbackground=True,
                backgroundcolor="rgba(255,255,255,0.54)",
                gridcolor="rgba(31,41,55,0.14)",
                zerolinecolor="rgba(31,41,55,0.18)",
            ),
            zaxis=dict(
                title="z",
                range=[0, 1],
                showbackground=True,
                backgroundcolor="rgba(255,255,255,0.54)",
                gridcolor="rgba(31,41,55,0.14)",
                zerolinecolor="rgba(31,41,55,0.18)",
            ),
            aspectmode="cube",
            camera=dict(
                eye=dict(x=1.55, y=1.45, z=1.18),
                center=dict(x=0, y=0, z=0),
                up=dict(x=0, y=0, z=1),
            ),
        ),
        annotations=[
            dict(
                text=(
                    f"Visible filter: length >= {min_length}, lattice diameter >= {min_diameter}. "
                    f"Shown steps: {shown_steps:,}. Largest shown loop: {largest:,} steps. "
                    "Drag to rotate, scroll to zoom, shift-drag to pan."
                ),
                x=0.02,
                y=0.02,
                xref="paper",
                yref="paper",
                showarrow=False,
                align="left",
                bgcolor="rgba(255,255,255,0.78)",
                bordercolor="rgba(31,41,55,0.16)",
                borderwidth=1,
                borderpad=8,
                font=dict(size=13, color="#374151"),
            )
        ],
    )
    return fig, shown_count, shown_steps


def parse_args():
    parser = argparse.ArgumentParser(description="Create a polished Plotly 3D viewer for a RWLS pickle.")
    parser.add_argument("pickle_path", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--min-length", type=int, default=8)
    parser.add_argument("--min-diameter", type=float, default=0.0)
    parser.add_argument("--max-loops", type=int, default=12000)
    parser.add_argument("--max-points", type=int, default=260000)
    parser.add_argument("--buckets", type=int, default=7)
    parser.add_argument(
        "--cdn",
        action="store_true",
        help="Use Plotly from CDN instead of embedding it. Smaller file, needs internet in browser.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    payload = load_payload(args.pickle_path)
    output = args.output or args.pickle_path.with_name(args.pickle_path.stem + "-plotly.html")

    fig, shown_count, shown_steps = make_figure(
        payload,
        min_length=args.min_length,
        min_diameter=args.min_diameter,
        max_loops=args.max_loops,
        max_points=args.max_points,
        bucket_count=args.buckets,
    )

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

    print(f"Shown loops: {shown_count}")
    print(f"Shown steps: {shown_steps}")
    print(f"Saved Plotly viewer to {output}")


if __name__ == "__main__":
    main()
