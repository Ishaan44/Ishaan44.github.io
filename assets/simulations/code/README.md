# Loop-soup simulation sources

Python sources for the square-lattice 2D experiments, circular visualizations,
and 3D cube experiments shown in the simulation gallery. Saved simulations and
exported HTML viewers are distributed separately. This bundle preserves the
original stochastic logic and parameter defaults. Limited repairs make the 2D
export paths relative to each script and correct escaped template braces in the
basic 3D Canvas viewer generator. Braces are corrected before inserting the
simulation data, preserving the exported coordinates.

## Setup

Use Python 3.10 or later. Run these commands from the extracted bundle directory:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

On Windows, activate the environment with `.venv\Scripts\activate` in Command
Prompt. All commands below assume this environment is active.

## 2D scripts

The 2D exporters take their parameters from constants near the bottom of each
file, not command-line flags. Before a first run, set `N = 50` and, for an
animation, `n_steps = 5` in the selected script. For `RWLS_single_loop.py`, use
`N = 50` and `N_STEPS = 5` in its parameter block. Some unchanged defaults are
much larger (including `N = 10000` in `square/RWLS.py`) and can take substantial
time and memory. Small runs may have no eligible loops or transitions for a
chosen cutoff.

Run any desired exporter; these are independent commands, not a required sequence:

```sh
python Loop-soup/square/RWLS.py
python Loop-soup/square/RWLS_cutoff.py
python Loop-soup/square/RWLS_epsilon_compare.py
python Loop-soup/square/RWLS_rewiring.py
python Loop-soup/square/RWLS_rewiring_compare.py
python Loop-soup/square/RWLS_distance.py
python Loop-soup/square/RWLS_distance_compare.py
python Loop-soup/square/RWLS_single_loop.py
```

For circular visualizations, the same commands are available in `Loop-soup/Disk/`
except `RWLS_distance_compare.py` and `RWLS_single_loop.py`. For example:

```sh
python Loop-soup/Disk/RWLS.py
python Loop-soup/Disk/RWLS_rewiring.py
```

PNG, GIF and pickle exports are saved beside their script. Repeating an export
with the same output filename replaces the previous file. The base disk carpet
is conformally mapped from the square; the other circular views clip paths
generated on the square lattice to a circle.

### Desktop distance viewer

This viewer requires a local graphical desktop and a Python installation with
Tk support. Reduce `N` and `n_steps` in the file before a first run:

```sh
python Loop-soup/square/RWLS_distance_interactive.py
```

Click a panel to relocate the observation point; use the pause/play button to
control animation. Tk is an operating-system/Python GUI dependency, not a pip
package. The file currently selects Matplotlib's `TkAgg` backend.

### Streamlit apps

```sh
python -m streamlit run Loop-soup/square/app.py --server.port 8501
python -m streamlit run Loop-soup/app.py --server.port 8502
```

The first app displays the square, and the second displays circularly clipped
views. Open the local URL printed by Streamlit. Start with grid size 50.

## 3D commands

The 3D scripts have command-line options; run a script with `--help` to inspect
them. This small example generates data and exports three alternative viewers:

```sh
python rwls3d/RWLS_3D.py --N 20 --c 0.5 --seed 11 --output rwls3d/output/example.pkl
python rwls3d/view_rwls3d.py rwls3d/output/example.pkl --output rwls3d/output/basic.html
python rwls3d/view_rwls3d_plotly.py rwls3d/output/example.pkl --output rwls3d/output/plotly.html
python rwls3d/view_rwls3d_two_clusters.py rwls3d/output/example.pkl --output rwls3d/output/clusters.html
```

The explicit `--output` above creates `rwls3d/output/example.pkl`; without that
option the N=20, c=0.5 generator writes `rwls3d/output/loopsoup-3d-N20-c0p5.pkl`.

The cluster viewer requires at least two eligible clusters; that condition is
sample-dependent. Open generated HTML files in a browser. Plotly exporters embed
their plotting library by default; choosing `--cdn` requires an internet
connection when viewing. The basic viewer is a standalone Canvas implementation
with no external browser library.

For a retained-largest-loops run, explicitly choose a modest cube size:

```sh
python rwls3d/RWLS_3D_visible.py --N 20 --c 0.5 --seed 11 --min-length 8 --min-diameter 2 --max-loops 1000 --output rwls3d/output/visible.pkl
```

For an independent small-loop sprinkling experiment:

```sh
python rwls3d/sprinkling_gluing.py --N 20 --c 0.5 --macro-seed 11 --dust-seed 29 --macro-cutoff 8 --epsilon-min 2 --output rwls3d/output/sprinkling.html
```

This experiment requires two eligible macro clusters; it may report that no
gluing occurred within the requested cutoff window. It is a stochastic outcome,
not a guaranteed connection. `--same-soup-dust` and `--local-dust` select the
other implemented experiment modes.

### Local Flask app

```sh
python rwls3d/app.py --host 127.0.0.1 --port 5055
```

Open `http://127.0.0.1:5055` and start with a small cube, such as `N = 20`.
This app computes new simulations locally; GitHub Pages only serves the archived
visuals and source downloads.

## Interpretation and platform notes

- Archived visuals are not guaranteed to have been generated using the current
  source defaults. Most 2D scripts do not record a seed or run manifest.
- In this implementation, `c` is used in the comparison `uniform() < c` when
  labelling/retaining loops. The 3D CLI restricts it to `[0, 1]`; some original
  2D app controls and defaults exceed 1, where this comparison saturates. The
  bundle preserves those controls. Treat them as code parameters, not a
  validated extension of the loop-soup intensity beyond 1.
- The 2D diameter cutoff is the maximum coordinate span divided by `N`.
  The distance observable is the maximum, over loops, of the nearest-vertex
  Euclidean distance to the observation point; vertex coordinates are divided
  by `N - 1`. It is not the distance to the union of all loops.
- The single-loop script selects the largest retained loop, despite older
  comments that describe random selection.
- `RWLS_rewiring_compare.py` in both 2D folders and
  `square/RWLS_distance_compare.py` explicitly use multiprocessing `fork`.
  Run those unchanged scripts on Linux or macOS; Windows needs a Linux
  environment or a multiprocessing adaptation.
- NumPy, SciPy, Matplotlib, Pillow, Streamlit, Plotly and Flask are included in
  the consolidated requirements. Versions are not pinned in the original
  projects. Sources are suitable for inspection without running simulations.
