import argparse
import json
import random
from argparse import Namespace
from pathlib import Path

import plotly.io as pio
from flask import Flask, Response, jsonify, request, send_file

from RWLS_3D import loop_soup_3d, loop_stats
from sprinkling_gluing import (
    build_figure as make_sprinkling_figure,
    choose_closest_macro_clusters,
    clipped_translated_local_dust,
    local_box_from_clusters,
    make_records,
    sprinkle_until_glued,
)
from view_rwls3d_plotly import make_figure as make_loop_soup_figure
from view_rwls3d_two_clusters import make_figure as make_cluster_figure


APP_ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = APP_ROOT / "output"

app = Flask(__name__)


INDEX_HTML = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>3D RWLS Sampler</title>
  <script src="/plotly.min.js"></script>
  <style>
    :root {
      color-scheme: light;
      --ink: #1f2937;
      --muted: #5f6673;
      --line: rgba(31, 41, 55, 0.16);
      --surface: rgba(255, 255, 255, 0.88);
      --paper: #f7f5ef;
      --accent: #0f766e;
      --accent-strong: #115e59;
      --danger: #b42318;
    }

    * {
      box-sizing: border-box;
    }

    html,
    body {
      margin: 0;
      min-height: 100%;
      background: var(--paper);
      color: var(--ink);
      font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    }

    body {
      overflow: hidden;
    }

    #app {
      display: grid;
      grid-template-columns: 320px minmax(0, 1fr);
      height: 100vh;
    }

    aside {
      border-right: 1px solid var(--line);
      background: var(--surface);
      backdrop-filter: blur(12px);
      padding: 18px;
      overflow-y: auto;
    }

    main {
      position: relative;
      min-width: 0;
      min-height: 0;
    }

    h1 {
      margin: 0 0 16px;
      font-size: 19px;
      line-height: 1.15;
      letter-spacing: 0;
    }

    .field {
      display: grid;
      gap: 7px;
      margin-bottom: 14px;
    }

    label {
      color: var(--muted);
      font-size: 12px;
      font-weight: 700;
    }

    input,
    select {
      width: 100%;
      height: 36px;
      border: 1px solid var(--line);
      border-radius: 6px;
      background: #fff;
      color: var(--ink);
      font: inherit;
      font-size: 14px;
      padding: 0 10px;
      outline: none;
    }

    input:focus,
    select:focus {
      border-color: rgba(15, 118, 110, 0.72);
      box-shadow: 0 0 0 3px rgba(15, 118, 110, 0.13);
    }

    .grid-2 {
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 10px;
    }

    .actions {
      display: grid;
      grid-template-columns: 1fr;
      gap: 10px;
      margin-top: 18px;
    }

    button {
      border: 1px solid transparent;
      border-radius: 6px;
      height: 40px;
      padding: 0 12px;
      background: var(--accent);
      color: #fff;
      font: inherit;
      font-size: 14px;
      font-weight: 760;
      cursor: pointer;
    }

    button:hover {
      background: var(--accent-strong);
    }

    button:disabled {
      cursor: wait;
      opacity: 0.72;
    }

    .secondary {
      background: #fff;
      color: var(--ink);
      border-color: var(--line);
    }

    .secondary:hover {
      background: #f0ede6;
    }

    .status {
      min-height: 44px;
      margin-top: 14px;
      padding-top: 12px;
      border-top: 1px solid var(--line);
      color: var(--muted);
      font-size: 12px;
      line-height: 1.45;
    }

    .status strong {
      color: var(--ink);
    }

    .error {
      color: var(--danger);
      font-weight: 700;
    }

    #plot {
      width: 100%;
      height: 100vh;
    }

    .overlay {
      position: absolute;
      inset: 0;
      display: none;
      place-items: center;
      background: rgba(247, 245, 239, 0.56);
      pointer-events: none;
      z-index: 4;
    }

    .overlay.active {
      display: grid;
    }

    .camera-tools {
      position: absolute;
      right: 16px;
      top: 16px;
      z-index: 3;
      display: flex;
      gap: 8px;
      padding: 8px;
      border: 1px solid var(--line);
      border-radius: 8px;
      background: rgba(255, 255, 255, 0.84);
      box-shadow: 0 14px 40px rgba(31, 41, 55, 0.12);
      backdrop-filter: blur(10px);
    }

    .camera-tools button {
      width: 38px;
      height: 34px;
      padding: 0;
      background: #fff;
      color: var(--ink);
      border-color: var(--line);
      font-size: 18px;
      line-height: 1;
    }

    .camera-tools button:hover {
      background: #f0ede6;
    }

    .camera-tools .wide {
      width: auto;
      padding: 0 10px;
      font-size: 12px;
    }

    .loader {
      padding: 13px 16px;
      border: 1px solid var(--line);
      border-radius: 8px;
      background: rgba(255, 255, 255, 0.92);
      color: var(--ink);
      box-shadow: 0 18px 50px rgba(31, 41, 55, 0.13);
      font-size: 14px;
      font-weight: 700;
    }

    @media (max-width: 820px) {
      body {
        overflow: auto;
      }

      #app {
        grid-template-columns: 1fr;
        height: auto;
      }

      aside {
        border-right: 0;
        border-bottom: 1px solid var(--line);
      }

      #plot {
        height: 72vh;
      }
    }
  </style>
</head>
<body>
  <div id="app">
    <aside>
      <h1>3D RWLS Sampler</h1>

      <div class="field">
        <label for="mode">View</label>
        <select id="mode">
          <option value="sprinkling">Sprinkling</option>
          <option value="clusters">Two Clusters</option>
          <option value="soup">Loop Soup</option>
        </select>
      </div>

      <div class="grid-2">
        <div class="field">
          <label for="N">N</label>
          <input id="N" type="number" min="5" max="240" step="1" value="90">
        </div>
        <div class="field">
          <label for="c">c</label>
          <input id="c" type="number" min="0" max="1" step="0.05" value="0.5">
        </div>
      </div>

      <div class="grid-2">
        <div class="field">
          <label for="seed">Seed</label>
          <input id="seed" type="number" step="1" value="">
        </div>
        <div class="field">
          <label for="minLength">Min Length</label>
          <input id="minLength" type="number" min="2" step="1" value="8">
        </div>
      </div>

      <div class="grid-2">
        <div class="field">
          <label for="macroCutoff">Macro Cutoff L</label>
          <input id="macroCutoff" type="number" min="1" step="1" value="8">
        </div>
        <div class="field">
          <label for="epsilonMin">Epsilon Min</label>
          <input id="epsilonMin" type="number" min="0" step="1" value="2">
        </div>
      </div>

      <div class="grid-2">
        <div class="field">
          <label for="cutoffMetric">Cutoff Metric</label>
          <select id="cutoffMetric">
            <option value="length">Time Length</option>
            <option value="diameter">Diameter</option>
          </select>
        </div>
        <div class="field">
          <label for="minMacroClusterLoops">Min Macro Cluster</label>
          <input id="minMacroClusterLoops" type="number" min="1" step="1" value="8">
        </div>
      </div>

      <div class="grid-2">
        <div class="field">
          <label for="dustC">Dust c</label>
          <input id="dustC" type="number" min="0" max="1" step="0.05" value="1">
        </div>
        <div class="field">
          <label for="dustSeed">Dust Seed</label>
          <input id="dustSeed" type="number" step="1" value="">
        </div>
      </div>

      <div class="grid-2">
        <div class="field">
          <label for="localDust">Local Dust</label>
          <select id="localDust">
            <option value="true">Box Around Clusters</option>
            <option value="false">Whole Cube</option>
          </select>
        </div>
        <div class="field">
          <label for="localPadding">Box Padding</label>
          <input id="localPadding" type="number" min="1" step="1" value="8">
        </div>
      </div>

      <div class="grid-2">
        <div class="field">
          <label for="minDiameter">Min Diameter</label>
          <input id="minDiameter" type="number" min="0" step="1" value="0">
        </div>
        <div class="field">
          <label for="maxPoints">Max Points</label>
          <input id="maxPoints" type="number" min="10000" step="10000" value="220000">
        </div>
      </div>

      <div class="grid-2">
        <div class="field">
          <label for="maxLoops">Max Loops</label>
          <input id="maxLoops" type="number" min="100" step="100" value="5000">
        </div>
        <div class="field">
          <label for="rankBy">Cluster Choice</label>
          <select id="rankBy">
            <option value="closest">Closest Pair</option>
            <option value="steps">Steps</option>
            <option value="diameter">Diameter</option>
            <option value="loops">Loop Count</option>
            <option value="max-loop">Max Loop</option>
          </select>
        </div>
      </div>

      <div class="actions">
        <button id="sample" type="button">New Sample</button>
        <button id="sampleSeeded" class="secondary" type="button">Resample Same Seed</button>
        <button id="lowerEpsilon" class="secondary" type="button">Lower Epsilon</button>
      </div>

      <div id="status" class="status">Ready.</div>
    </aside>

    <main>
      <div id="plot"></div>
      <div class="camera-tools" aria-label="camera controls">
        <button id="zoomIn" type="button" title="Zoom in">+</button>
        <button id="zoomOut" type="button" title="Zoom out">-</button>
        <button id="resetCamera" class="wide" type="button" title="Reset camera">Reset</button>
      </div>
      <div id="overlay" class="overlay">
        <div class="loader">Generating sample...</div>
      </div>
    </main>
  </div>

  <script>
    const controls = {
      mode: document.getElementById("mode"),
      N: document.getElementById("N"),
      c: document.getElementById("c"),
      seed: document.getElementById("seed"),
      minLength: document.getElementById("minLength"),
      macroCutoff: document.getElementById("macroCutoff"),
      epsilonMin: document.getElementById("epsilonMin"),
      cutoffMetric: document.getElementById("cutoffMetric"),
      minMacroClusterLoops: document.getElementById("minMacroClusterLoops"),
      dustC: document.getElementById("dustC"),
      dustSeed: document.getElementById("dustSeed"),
      localDust: document.getElementById("localDust"),
      localPadding: document.getElementById("localPadding"),
      minDiameter: document.getElementById("minDiameter"),
      maxPoints: document.getElementById("maxPoints"),
      maxLoops: document.getElementById("maxLoops"),
      rankBy: document.getElementById("rankBy"),
    };

    const status = document.getElementById("status");
    const overlay = document.getElementById("overlay");
    const sampleButton = document.getElementById("sample");
    const sampleSeededButton = document.getElementById("sampleSeeded");
    const lowerEpsilonButton = document.getElementById("lowerEpsilon");
    const zoomInButton = document.getElementById("zoomIn");
    const zoomOutButton = document.getElementById("zoomOut");
    const resetCameraButton = document.getElementById("resetCamera");
    const plotDiv = document.getElementById("plot");
    const defaultCamera = {
      eye: {x: 1.55, y: 1.42, z: 1.16},
      center: {x: 0, y: 0, z: 0},
      up: {x: 0, y: 0, z: 1},
    };
    let currentCamera = cloneCamera(defaultCamera);
    let relayoutListenerAttached = false;

    function cloneCamera(camera) {
      return JSON.parse(JSON.stringify(camera));
    }

    function integerValue(id) {
      const value = controls[id].value;
      return value === "" ? null : Number.parseInt(value, 10);
    }

    function floatValue(id) {
      return Number.parseFloat(controls[id].value);
    }

    function payload(forceSeed) {
      let seed = integerValue("seed");
      if (!forceSeed || seed === null || Number.isNaN(seed)) {
        seed = Math.floor(Math.random() * 2147483647);
        controls.seed.value = seed;
      }

      let dustSeed = integerValue("dustSeed");
      if (!forceSeed || dustSeed === null || Number.isNaN(dustSeed)) {
        dustSeed = Math.floor(Math.random() * 2147483647);
        controls.dustSeed.value = dustSeed;
      }

      return {
        mode: controls.mode.value,
        N: integerValue("N"),
        c: floatValue("c"),
        seed,
        dustSeed,
        macroCutoff: floatValue("macroCutoff"),
        epsilonMin: floatValue("epsilonMin"),
        cutoffMetric: controls.cutoffMetric.value,
        minMacroClusterLoops: integerValue("minMacroClusterLoops"),
        dustC: floatValue("dustC"),
        localDust: controls.localDust.value === "true",
        localPadding: integerValue("localPadding"),
        minLength: integerValue("minLength"),
        minDiameter: floatValue("minDiameter"),
        maxPoints: integerValue("maxPoints"),
        maxLoops: integerValue("maxLoops"),
        rankBy: controls.rankBy.value,
      };
    }

    function setBusy(isBusy) {
      overlay.classList.toggle("active", isBusy);
      sampleButton.disabled = isBusy;
      sampleSeededButton.disabled = isBusy;
      lowerEpsilonButton.disabled = isBusy;
    }

    function updateStatus(meta) {
      const lines = [
        `<strong>${meta.mode_label}</strong>`,
        `N=${meta.N}, c=${meta.c}, seed=${meta.seed}`,
        `source loops=${meta.source_loops.toLocaleString()}, shown loops=${meta.shown_loops.toLocaleString()}`,
      ];
      if (meta.closest_distance !== undefined && meta.closest_distance !== null) {
        lines.push(`closest distance=${meta.closest_distance_rescaled.toFixed(4)} rescaled (${meta.closest_distance.toFixed(2)} lattice units)`);
      }
      if (meta.gluing_cutoff !== undefined && meta.gluing_cutoff !== null) {
        lines.push(`<strong>CONNECTED</strong> at epsilon*=${meta.gluing_cutoff.toFixed(2)} ${meta.cutoff_metric_label}`);
      } else if (meta.mode_label === "Sprinkling") {
        lines.push(`<strong>NOT CONNECTED</strong> down to epsilon=${meta.epsilon_min} ${meta.cutoff_metric_label}`);
      }
      if (meta.dust_loops !== undefined) {
        lines.push(`dust loops=${meta.dust_loops.toLocaleString()}, activated=${meta.activated_dust_loops.toLocaleString()}, bridge=${meta.bridge_dust_loops.toLocaleString()}`);
      }
      if (meta.local_box) {
        lines.push(`local dust box=${meta.local_box.shape.join("x")} with padding ${meta.local_padding}`);
      }
      lines.push(`elapsed=${meta.elapsed_seconds.toFixed(2)}s`);
      status.innerHTML = lines.join("<br>");
    }

    async function generate(forceSeed = false) {
      const body = payload(forceSeed);
      setBusy(true);
      status.textContent = `Generating N=${body.N} sample...`;

      try {
        const response = await fetch("/api/sample", {
          method: "POST",
          headers: {"Content-Type": "application/json"},
          body: JSON.stringify(body),
        });

        const result = await response.json();
        if (!response.ok || result.error) {
          throw new Error(result.error || "Sampling failed.");
        }

        const fig = JSON.parse(result.figure);
        if (fig.layout && fig.layout.scene && fig.layout.scene.camera) {
          currentCamera = cloneCamera(fig.layout.scene.camera);
        } else {
          currentCamera = cloneCamera(defaultCamera);
        }
        await Plotly.react(plotDiv, fig.data, fig.layout, {
          responsive: true,
          displaylogo: false,
          scrollZoom: true,
          modeBarButtonsToRemove: ["lasso2d", "select2d"],
        });
        if (!relayoutListenerAttached) {
          plotDiv.on("plotly_relayout", (event) => {
            if (event["scene.camera"]) {
              currentCamera = cloneCamera(event["scene.camera"]);
            }
          });
          relayoutListenerAttached = true;
        }
        updateStatus(result.meta);
      } catch (error) {
        status.innerHTML = `<span class="error">${error.message}</span>`;
      } finally {
        setBusy(false);
      }
    }

    function zoomCamera(factor) {
      const camera = currentCamera || cloneCamera(defaultCamera);
      const eye = camera.eye || defaultCamera.eye;
      currentCamera = {
        ...camera,
        eye: {
          x: eye.x * factor,
          y: eye.y * factor,
          z: eye.z * factor,
        },
      };
      Plotly.relayout(plotDiv, {"scene.camera": currentCamera});
    }

    sampleButton.addEventListener("click", () => generate(false));
    sampleSeededButton.addEventListener("click", () => generate(true));
    lowerEpsilonButton.addEventListener("click", () => {
      const current = Number.parseFloat(controls.epsilonMin.value);
      const step = controls.cutoffMetric.value === "diameter" ? 0.5 : 1;
      const next = Math.max(0, current - step);
      controls.epsilonMin.value = Number.isInteger(next) ? String(next) : next.toFixed(1);
      generate(true);
    });
    zoomInButton.addEventListener("click", () => zoomCamera(0.78));
    zoomOutButton.addEventListener("click", () => zoomCamera(1.28));
    resetCameraButton.addEventListener("click", () => {
      currentCamera = cloneCamera(defaultCamera);
      Plotly.relayout(plotDiv, {"scene.camera": currentCamera});
    });
    window.addEventListener("resize", () => Plotly.Plots.resize(plotDiv));

    generate(false);
  </script>
</body>
</html>
"""


def plotly_js_path():
    import plotly

    return Path(plotly.__file__).resolve().parent / "package_data" / "plotly.min.js"


@app.get("/")
def index():
    return Response(INDEX_HTML, mimetype="text/html")


@app.get("/plotly.min.js")
def plotly_min_js():
    return send_file(plotly_js_path(), mimetype="application/javascript")


def clamp(value, low, high):
    return max(low, min(high, value))


def build_payload(data):
    N = clamp(int(data.get("N", 90)), 5, 240)
    c = clamp(float(data.get("c", 0.5)), 0.0, 1.0)
    seed = data.get("seed")
    if seed is None:
        seed = random.randrange(2**31)
    seed = int(seed)

    return {
        "N": N,
        "c": c,
        "seed": seed,
        "min_length": max(2, int(data.get("minLength", 8))),
        "min_diameter": max(0.0, float(data.get("minDiameter", 0.0))),
        "max_points": max(10_000, int(data.get("maxPoints", 220_000))),
        "max_loops": max(100, int(data.get("maxLoops", 5_000))),
        "mode": data.get("mode", "clusters"),
        "rank_by": data.get("rankBy", "steps"),
        "macro_cutoff": max(0.0, float(data.get("macroCutoff", 8.0))),
        "epsilon_min": max(0.0, float(data.get("epsilonMin", 1.0))),
        "cutoff_metric": data.get("cutoffMetric", "length"),
        "min_macro_cluster_loops": max(1, int(data.get("minMacroClusterLoops", 8))),
        "dust_c": clamp(float(data.get("dustC", 1.0)), 0.0, 1.0),
        "dust_seed": int(data.get("dustSeed", random.randrange(2**31))),
        "local_dust": bool(data.get("localDust", True)),
        "local_padding": max(1, int(data.get("localPadding", 8))),
    }


def source_payload(loopsoup, params):
    return {
        "N": params["N"],
        "c": params["c"],
        "seed": params["seed"],
        "dimension": 3,
        "loopsoup": loopsoup,
        "stats": loop_stats(loopsoup),
    }


def make_soup_view(payload, params):
    fig, shown_count, shown_steps = make_loop_soup_figure(
        payload,
        min_length=params["min_length"],
        min_diameter=params["min_diameter"],
        max_loops=params["max_loops"],
        max_points=params["max_points"],
        bucket_count=7,
    )

    return fig, {
        "mode_label": "Loop Soup",
        "shown_loops": shown_count,
        "shown_steps": shown_steps,
    }


def make_cluster_view(payload, params):
    cluster_args = Namespace(
        clusters=2,
        min_length=params["min_length"],
        min_diameter=params["min_diameter"],
        min_cluster_loops=1,
        rank_by=params["rank_by"],
        max_points_per_cluster=params["max_points"],
        line_width=7.0,
        height=900,
        show_boxes=True,
        cdn=False,
    )
    fig, filtered_count, cluster_count, selected_clusters = make_cluster_figure(payload, cluster_args)
    shown_loops = sum(cluster["stats"]["loop_count"] for cluster in selected_clusters)
    shown_steps = sum(cluster["stats"]["total_steps"] for cluster in selected_clusters)
    closest_distance = selected_clusters[0].get("selection_distance")

    fig.update_layout(
        title=dict(
            text=(
                f"{'Two Closest Clusters' if params['rank_by'] == 'closest' else 'Two Intersection Clusters'} "
                f"(N={params['N']}, c={params['c']}, seed={params['seed']})"
            )
        )
    )

    return fig, {
        "mode_label": "Two Clusters",
        "shown_loops": shown_loops,
        "shown_steps": shown_steps,
        "filtered_loops": filtered_count,
        "cluster_count": cluster_count,
        "closest_distance": closest_distance,
        "closest_distance_rescaled": (
            closest_distance / max(1, params["N"] - 1)
            if closest_distance is not None
            else None
        ),
    }


def make_sprinkling_view(macro_payload, params):
    macro_records = make_records(
        macro_payload["loopsoup"],
        min_size=params["macro_cutoff"],
        cutoff_metric=params["cutoff_metric"],
    )

    if len(macro_records) < 2:
        raise ValueError("Not enough macro loops. Lower Macro Cutoff L.")

    selected_clusters, closest_distance, cluster_count, eligible_count = choose_closest_macro_clusters(
        macro_records,
        min_macro_cluster_loops=params["min_macro_cluster_loops"],
    )
    dust_loops, local_box = generate_sprinkling_dust(macro_records, selected_clusters, params)
    dust_records = make_records(
        dust_loops,
        min_size=params["epsilon_min"],
        max_size=params["macro_cutoff"],
        cutoff_metric=params["cutoff_metric"],
    )

    if not dust_records:
        raise ValueError("No dust loops in the epsilon window. Lower Epsilon Min or raise Macro Cutoff L.")

    sprinkle_result = sprinkle_until_glued(
        dust_records,
        macro_records,
        selected_clusters[0],
        selected_clusters[1],
        params["epsilon_min"],
    )

    fig = make_sprinkling_figure(
        {
            "macro_c": params["c"],
            "dust_c": params["dust_c"],
            "macro_seed": params["seed"],
            "dust_seed": params["dust_seed"],
        },
        macro_records,
        selected_clusters,
        sprinkle_result,
        params["N"],
        params["macro_cutoff"],
        params["epsilon_min"],
        params["cutoff_metric"],
        show_boxes=True,
        show_context_dust=True,
        max_context_points=params["max_points"],
        height=900,
    )

    gluing_cutoff = sprinkle_result["gluing_cutoff"]
    cutoff_metric_label = "time steps" if params["cutoff_metric"] == "length" else "lattice diameter"

    return fig, {
        "mode_label": "Sprinkling",
        "shown_loops": sum(cluster["stats"]["loop_count"] for cluster in selected_clusters)
        + len(sprinkle_result["activated_dust_indices"]),
        "shown_steps": sum(cluster["stats"]["total_steps"] for cluster in selected_clusters),
        "filtered_loops": len(macro_records),
        "cluster_count": cluster_count,
        "eligible_cluster_count": eligible_count,
        "closest_distance": closest_distance,
        "closest_distance_rescaled": closest_distance / max(1, params["N"] - 1),
        "gluing_cutoff": gluing_cutoff,
        "gluing_cutoff_rescaled": (
            gluing_cutoff / max(1, params["N"] - 1)
            if gluing_cutoff is not None and params["cutoff_metric"] == "diameter"
            else None
        ),
        "epsilon_min": params["epsilon_min"],
        "cutoff_metric": params["cutoff_metric"],
        "cutoff_metric_label": cutoff_metric_label,
        "min_macro_cluster_loops": params["min_macro_cluster_loops"],
        "dust_loops": len(dust_records),
        "bridge_dust_loops": len(sprinkle_result["bridge_dust_indices"]),
        "activated_dust_loops": len(sprinkle_result["activated_dust_indices"]),
        "local_box": local_box,
        "local_padding": params["local_padding"],
    }


def generate_sprinkling_dust(macro_records, selected_clusters, params):
    if not params["local_dust"]:
        dust_loops = loop_soup_3d(
            params["N"],
            params["dust_c"],
            seed=params["dust_seed"],
            progress_every=0,
        )
        return dust_loops, None

    lower, upper, box_shape = local_box_from_clusters(
        macro_records,
        selected_clusters,
        params["N"],
        params["local_padding"],
    )
    dust_loops = clipped_translated_local_dust(
        box_shape,
        lower,
        params["dust_c"],
        params["dust_seed"],
        progress_every=0,
    )
    return dust_loops, {"lower": lower, "upper": upper, "shape": box_shape}


@app.post("/api/sample")
def api_sample():
    import time

    started = time.time()
    params = build_payload(request.get_json(force=True) or {})

    try:
        loopsoup = loop_soup_3d(params["N"], params["c"], seed=params["seed"], progress_every=0)
        payload = source_payload(loopsoup, params)

        if params["mode"] == "soup":
            fig, view_meta = make_soup_view(payload, params)
        elif params["mode"] == "sprinkling":
            fig, view_meta = make_sprinkling_view(payload, params)
        else:
            fig, view_meta = make_cluster_view(payload, params)

        elapsed = time.time() - started
        meta = {
            "N": params["N"],
            "c": params["c"],
            "seed": params["seed"],
            "source_loops": len(loopsoup),
            "elapsed_seconds": elapsed,
            **view_meta,
        }

        return jsonify({"figure": pio.to_json(fig, validate=False), "meta": meta})
    except Exception as exc:
        return jsonify({"error": str(exc)}), 400


def parse_args():
    parser = argparse.ArgumentParser(description="Local interactive 3D RWLS sampler.")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=5055)
    parser.add_argument("--debug", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    print(f"Open http://{args.host}:{args.port} in your browser")
    app.run(host=args.host, port=args.port, debug=args.debug)


if __name__ == "__main__":
    main()
