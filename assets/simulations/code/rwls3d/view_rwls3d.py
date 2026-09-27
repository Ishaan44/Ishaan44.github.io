import argparse
import json
import pickle
from pathlib import Path

import numpy as np


HTML_TEMPLATE = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>3D RWLS Viewer</title>
  <style>
    html, body {{
      margin: 0;
      width: 100%;
      height: 100%;
      overflow: hidden;
      background: #f6f4ef;
      color: #1d1c19;
      font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    }}

    #stage {{
      display: block;
      width: 100vw;
      height: 100vh;
      cursor: grab;
      touch-action: none;
    }}

    #stage:active {{
      cursor: grabbing;
    }}

    .panel {{
      position: fixed;
      top: 16px;
      left: 16px;
      max-width: min(360px, calc(100vw - 32px));
      padding: 12px 14px;
      background: rgba(255, 255, 255, 0.84);
      border: 1px solid rgba(35, 33, 29, 0.14);
      border-radius: 8px;
      box-shadow: 0 12px 36px rgba(35, 33, 29, 0.12);
      backdrop-filter: blur(10px);
      user-select: none;
    }}

    .title {{
      font-weight: 760;
      font-size: 14px;
      line-height: 1.2;
      margin-bottom: 8px;
    }}

    .stats {{
      display: grid;
      grid-template-columns: auto auto;
      gap: 4px 16px;
      font-size: 12px;
      line-height: 1.35;
      color: #514d46;
    }}

    .stats strong {{
      color: #1d1c19;
      font-weight: 680;
      text-align: right;
    }}

    .controls {{
      display: flex;
      align-items: center;
      gap: 8px;
      margin-top: 12px;
    }}

    button {{
      appearance: none;
      border: 1px solid rgba(35, 33, 29, 0.18);
      background: #ffffff;
      color: #1d1c19;
      border-radius: 6px;
      font: inherit;
      font-size: 12px;
      font-weight: 650;
      height: 30px;
      padding: 0 10px;
      cursor: pointer;
    }}

    button:hover {{
      background: #f0ede6;
    }}

    input[type="range"] {{
      width: 120px;
      accent-color: #2f6f68;
    }}
  </style>
</head>
<body>
  <canvas id="stage"></canvas>
  <div class="panel">
    <div class="title">3D RWLS Cube</div>
    <div class="stats">
      <span>N</span><strong id="stat-n"></strong>
      <span>c</span><strong id="stat-c"></strong>
      <span>shown loops</span><strong id="stat-loops"></strong>
      <span>shown steps</span><strong id="stat-steps"></strong>
      <span>source loops</span><strong id="stat-source"></strong>
    </div>
    <div class="controls">
      <button id="reset" type="button">Reset</button>
      <input id="density" type="range" min="0.08" max="1" step="0.01" value="1" aria-label="density">
    </div>
  </div>
  <script>
    const DATA = __DATA__;

    const canvas = document.getElementById("stage");
    const ctx = canvas.getContext("2d");
    const densityControl = document.getElementById("density");

    let width = 0;
    let height = 0;
    let yaw = 0.68;
    let pitch = -0.48;
    let zoom = 1.62;
    let panX = 0;
    let panY = 0;
    let dragging = false;
    let panning = false;
    let lastX = 0;
    let lastY = 0;

    document.getElementById("stat-n").textContent = DATA.N;
    document.getElementById("stat-c").textContent = DATA.c;
    document.getElementById("stat-loops").textContent = DATA.loops.length.toLocaleString();
    document.getElementById("stat-steps").textContent = DATA.shown_steps.toLocaleString();
    document.getElementById("stat-source").textContent = DATA.source_loop_count.toLocaleString();

    document.getElementById("reset").addEventListener("click", () => {{
      yaw = 0.68;
      pitch = -0.48;
      zoom = 1.62;
      panX = 0;
      panY = 0;
      densityControl.value = "1";
      draw();
    }});

    function resize() {{
      const dpr = Math.max(1, Math.min(2, window.devicePixelRatio || 1));
      width = window.innerWidth;
      height = window.innerHeight;
      canvas.width = Math.floor(width * dpr);
      canvas.height = Math.floor(height * dpr);
      canvas.style.width = width + "px";
      canvas.style.height = height + "px";
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      draw();
    }}

    function rotatePoint(point) {{
      const x = point[0] - 0.5;
      const y = point[1] - 0.5;
      const z = point[2] - 0.5;

      const cy = Math.cos(yaw);
      const sy = Math.sin(yaw);
      const cp = Math.cos(pitch);
      const sp = Math.sin(pitch);

      const x1 = cy * x + sy * z;
      const z1 = -sy * x + cy * z;
      const y1 = cp * y - sp * z1;
      const z2 = sp * y + cp * z1;

      return [x1, y1, z2];
    }}

    function project(point) {{
      const rotated = rotatePoint(point);
      const depth = 2.35 - rotated[2];
      const perspective = 1.15 / Math.max(0.35, depth);
      const scale = Math.min(width, height) * 0.82 * zoom;

      return [
        width * 0.5 + panX + rotated[0] * scale * perspective,
        height * 0.52 + panY + rotated[1] * scale * perspective,
        rotated[2],
      ];
    }}

    function drawCube() {{
      const vertices = [
        [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
        [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
      ];
      const edges = [
        [0, 1], [1, 2], [2, 3], [3, 0],
        [4, 5], [5, 6], [6, 7], [7, 4],
        [0, 4], [1, 5], [2, 6], [3, 7],
      ];

      ctx.save();
      ctx.strokeStyle = "rgba(29, 28, 25, 0.24)";
      ctx.lineWidth = 1;
      ctx.beginPath();
      for (const edge of edges) {{
        const a = project(vertices[edge[0]]);
        const b = project(vertices[edge[1]]);
        ctx.moveTo(a[0], a[1]);
        ctx.lineTo(b[0], b[1]);
      }}
      ctx.stroke();
      ctx.restore();
    }}

    function loopColor(index, length) {{
      const t = Math.min(1, Math.log1p(length) / Math.log1p(DATA.max_loop_length || 1));
      const hue = 184 + 138 * t;
      return `hsla(${{hue}}, 58%, ${{24 + 18 * t}}%, ${{0.34 + 0.32 * t}})`;
    }}

    function drawLoops() {{
      const density = Number(densityControl.value);
      const stride = Math.max(1, Math.ceil(1 / density));

      ctx.save();
      ctx.lineCap = "round";
      ctx.lineJoin = "round";

      for (let i = 0; i < DATA.loops.length; i += stride) {{
        const loop = DATA.loops[i];
        if (loop.length < 2) continue;

        ctx.strokeStyle = loopColor(i, loop.length);
        ctx.lineWidth = Math.min(3.2, Math.max(0.75, 0.55 + Math.log1p(loop.length) * 0.13));
        ctx.beginPath();

        const first = project(loop[0]);
        ctx.moveTo(first[0], first[1]);

        for (let j = 1; j < loop.length; j++) {{
          const p = project(loop[j]);
          ctx.lineTo(p[0], p[1]);
        }}

        ctx.stroke();
      }}

      ctx.restore();
    }}

    function draw() {{
      ctx.clearRect(0, 0, width, height);
      ctx.fillStyle = "#f6f4ef";
      ctx.fillRect(0, 0, width, height);
      drawCube();
      drawLoops();
    }}

    canvas.addEventListener("pointerdown", (event) => {{
      dragging = true;
      panning = event.button === 2 || event.shiftKey;
      lastX = event.clientX;
      lastY = event.clientY;
      canvas.setPointerCapture(event.pointerId);
    }});

    canvas.addEventListener("pointermove", (event) => {{
      if (!dragging) return;

      const dx = event.clientX - lastX;
      const dy = event.clientY - lastY;
      lastX = event.clientX;
      lastY = event.clientY;

      if (panning) {{
        panX += dx;
        panY += dy;
      }} else {{
        yaw += dx * 0.008;
        pitch = Math.max(-1.45, Math.min(1.45, pitch + dy * 0.008));
      }}

      draw();
    }});

    canvas.addEventListener("pointerup", (event) => {{
      dragging = false;
      canvas.releasePointerCapture(event.pointerId);
    }});

    canvas.addEventListener("pointercancel", () => {{
      dragging = false;
    }});

    canvas.addEventListener("contextmenu", (event) => event.preventDefault());

    canvas.addEventListener("wheel", (event) => {{
      event.preventDefault();
      const factor = Math.exp(-event.deltaY * 0.0012);
      zoom = Math.max(0.28, Math.min(12, zoom * factor));
      draw();
    }}, {{ passive: false }});

    densityControl.addEventListener("input", draw);
    window.addEventListener("resize", resize);
    resize();
  </script>
</body>
</html>
"""


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


def loop_diameter(loop):
    points = np.asarray(loop, dtype=float)
    if len(points) == 0:
        return 0.0
    span = points.max(axis=0) - points.min(axis=0)
    return float(np.linalg.norm(span))


def filtered_loops(loops, min_length, min_diameter, max_loops, max_points):
    candidates = []

    for loop in loops:
        length = len(loop) - 1
        if length < min_length:
            continue

        diameter = loop_diameter(loop)
        if diameter < min_diameter:
            continue

        candidates.append((length, diameter, loop))

    candidates.sort(key=lambda item: (item[1], item[0]), reverse=True)

    selected = []
    total_points = 0

    for _, _, loop in candidates:
        if max_loops is not None and len(selected) >= max_loops:
            break
        if max_points is not None and total_points + len(loop) > max_points:
            break
        selected.append(loop)
        total_points += len(loop)

    return selected


def normalize_loop(loop, N):
    scale = max(1, N - 1)
    return [[round(p[0] / scale, 6), round(p[1] / scale, 6), round(p[2] / scale, 6)] for p in loop]


def build_view_data(payload, min_length, min_diameter, max_loops, max_points):
    loops = payload["loopsoup"]
    N = payload.get("N")
    if N is None:
        N = max(max(max(point) for point in loop) for loop in loops if loop) + 1

    selected = filtered_loops(loops, min_length, min_diameter, max_loops, max_points)
    normalized = [normalize_loop(loop, N) for loop in selected]
    lengths = [len(loop) - 1 for loop in selected]

    return {
        "N": N,
        "c": payload.get("c"),
        "seed": payload.get("seed"),
        "source_loop_count": len(loops),
        "loops": normalized,
        "shown_steps": int(sum(lengths)),
        "max_loop_length": int(max(lengths) if lengths else 0),
        "min_length": min_length,
        "min_diameter": min_diameter,
    }


def parse_args():
    parser = argparse.ArgumentParser(description="Export a 3D RWLS pickle to an interactive HTML viewer.")
    parser.add_argument("pickle_path", type=Path, help="Input .pkl generated by RWLS_3D.py.")
    parser.add_argument("--output", type=Path, default=None, help="Output .html path.")
    parser.add_argument("--min-length", type=int, default=2, help="Only show loops with at least this many steps.")
    parser.add_argument("--min-diameter", type=float, default=0.0, help="Only show loops with this lattice diameter or larger.")
    parser.add_argument("--max-loops", type=int, default=3000, help="Maximum number of loops to embed.")
    parser.add_argument("--max-points", type=int, default=160000, help="Maximum number of path points to embed.")
    return parser.parse_args()


def main():
    args = parse_args()
    payload = load_payload(args.pickle_path)
    view_data = build_view_data(
        payload,
        min_length=args.min_length,
        min_diameter=args.min_diameter,
        max_loops=args.max_loops,
        max_points=args.max_points,
    )

    output = args.output
    if output is None:
        output = args.pickle_path.with_suffix(".html")

    output.parent.mkdir(parents=True, exist_ok=True)
    # Unescape the template before inserting JSON so data is never rewritten.
    html = HTML_TEMPLATE.replace("{{", "{").replace("}}", "}")
    html = html.replace("__DATA__", json.dumps(view_data, separators=(",", ":")))
    output.write_text(html, encoding="utf-8")

    print(f"Source loops: {view_data['source_loop_count']}")
    print(f"Shown loops: {len(view_data['loops'])}")
    print(f"Shown steps: {view_data['shown_steps']}")
    print(f"Saved viewer to {output}")


if __name__ == "__main__":
    main()
