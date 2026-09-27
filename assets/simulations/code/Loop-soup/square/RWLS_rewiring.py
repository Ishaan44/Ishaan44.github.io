from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

def NextTip(N, tip):
    if tip[0] == 0:
        if tip[1] == 0:
            if np.random.randint(2) == 0:
                newtip = [1, 0]
            else:
                newtip = [0, 1]
        elif tip[1] == N - 1:
            if np.random.randint(2) == 0:
                newtip = [1, N - 1]
            else:
                newtip = [0, N - 2]
        else:
            direction = np.random.randint(3)
            if direction == 0:
                newtip = [1, tip[1]]
            elif direction == 1:
                newtip = [0, tip[1] + 1]
            else:
                newtip = [0, tip[1] - 1]
    elif tip[0] == N - 1:
        if tip[1] == 0:
            if np.random.randint(2) == 0:
                newtip = [N - 2, 0]
            else:
                newtip = [N - 1, 1]
        elif tip[1] == N - 1:
            if np.random.randint(2) == 0:
                newtip = [N - 2, N - 1]
            else:
                newtip = [N - 1, N - 2]
        else:
            direction = np.random.randint(3)
            if direction == 0:
                newtip = [N - 2, tip[1]]
            elif direction == 1:
                newtip = [N - 1, tip[1] + 1]
            else:
                newtip = [N - 1, tip[1] - 1]
    elif tip[1] == 0:
        direction = np.random.randint(3)
        if direction == 0:
            newtip = [tip[0], 1]
        elif direction == 1:
            newtip = [tip[0] + 1, 0]
        else:
            newtip = [tip[0] - 1, 0]
    elif tip[1] == N - 1:
        direction = np.random.randint(3)
        if direction == 0:
            newtip = [tip[0], N - 2]
        elif direction == 1:
            newtip = [tip[0] + 1, N - 1]
        else:
            newtip = [tip[0] - 1, N - 1]
    else:
        direction = np.random.randint(4)
        if direction == 0:
            newtip = [tip[0] + 1, tip[1]]
        elif direction == 1:
            newtip = [tip[0] - 1, tip[1]]
        elif direction == 2:
            newtip = [tip[0], tip[1] + 1]
        else:
            newtip = [tip[0], tip[1] - 1]
    return newtip

def LoopSoup(N, c, epsilon=0.0):
    """
    Generate loop soup, discarding loops with spatial diameter < epsilon on the fly.
    Wilson's algorithm still runs in full — the filter saves memory, not compute time.
    """
    loopsoup = []
    tree = np.zeros((N, N))
    for i in range(N):
        tree[i, 0] = 1
        tree[i, N - 1] = 1
        tree[0, i] = 1
        tree[N - 1, i] = 1
    for i in range(1, N - 1):
        for j in range(1, N - 1):
            if tree[i, j] == 1:
                continue
            tip = [i, j]
            path = [tip]
            while tree[tip[0], tip[1]] == 0:
                newtip = NextTip(N, tip)
                path.append(newtip)
                tip = newtip
            while path:
                root = path[0]
                tree[root[0], root[1]] = 1
                label = True
                rmax = 0
                last = len(path) - 1 - path[::-1].index(root)
                now = 0
                while now < last:
                    next = now + 1 + path[now + 1:].index(root)
                    if (r := np.random.uniform()) > rmax:
                        rmax = r
                        label = (np.random.uniform() < c)
                    if label:
                        loop = path[now:next + 1]
                        if spatial_diameter(loop, N) >= epsilon:
                            loopsoup.append(loop)
                    now = next
                path = path[last + 1:]
    return loopsoup

# ── Rewiring operations ──────────────────────────────────────────────────────

# ── Sparse table helpers for O(1) range min/max ───────────────────────────────

def _build_st_max(arr):
    """Sparse table for range-max. Build O(n log n), query O(1)."""
    n = len(arr)
    st = [arr]
    j = 1
    while (1 << j) <= n:
        half = 1 << (j - 1)
        end = n - (1 << j) + 1
        st.append(np.maximum(st[j-1][:end], st[j-1][half:half+end]))
        j += 1
    return st

def _build_st_min(arr):
    """Sparse table for range-min. Build O(n log n), query O(1)."""
    n = len(arr)
    st = [arr]
    j = 1
    while (1 << j) <= n:
        half = 1 << (j - 1)
        end = n - (1 << j) + 1
        st.append(np.minimum(st[j-1][:end], st[j-1][half:half+end]))
        j += 1
    return st

def _qmax(st, l, r):
    """O(1) range-max query over [l, r] inclusive."""
    k = (r - l + 1).bit_length() - 1
    return int(max(st[k][l], st[k][r - (1 << k) + 1]))

def _qmin(st, l, r):
    """O(1) range-min query over [l, r] inclusive."""
    k = (r - l + 1).bit_length() - 1
    return int(min(st[k][l], st[k][r - (1 << k) + 1]))

# ─────────────────────────────────────────────────────────────────────────────

def visit_times(loop):
    """Map each lattice site -> list of times the loop visits it."""
    visits = defaultdict(list)
    for t, p in enumerate(loop):
        visits[tuple(p)].append(t)
    return dict(visits)

def spatial_diameter(loop, N):
    """Bounding-box diameter in scaled [0,1] coordinates."""
    rows = [p[0] for p in loop]
    cols = [p[1] for p in loop]
    return max(max(rows) - min(rows), max(cols) - min(cols)) / N

def SI_pairs(loop, epsilon, N):
    """
    Self-intersection pairs (s,t) valid for a split.
    All bounding-box queries are O(1):
      - loop1 = loop[s:t+1]:         sparse table range max/min
      - loop2 = loop[t:]+loop[1:s+1]: suffix/prefix accumulated arrays
    """
    rows = np.array([p[0] for p in loop], dtype=np.int32)
    cols = np.array([p[1] for p in loop], dtype=np.int32)

    # Sparse tables for loop1 O(1) range queries
    st_max_r = _build_st_max(rows)
    st_min_r = _build_st_min(rows)
    st_max_c = _build_st_max(cols)
    st_min_c = _build_st_min(cols)

    # Suffix arrays for loop2: suf_max_r[i] = max(rows[i:])
    suf_max_r = np.maximum.accumulate(rows[::-1])[::-1]
    suf_min_r = np.minimum.accumulate(rows[::-1])[::-1]
    suf_max_c = np.maximum.accumulate(cols[::-1])[::-1]
    suf_min_c = np.minimum.accumulate(cols[::-1])[::-1]

    # Prefix arrays from index 1: pre_max_r[i] = max(rows[1:i+2])
    pre_max_r = np.maximum.accumulate(rows[1:])
    pre_min_r = np.minimum.accumulate(rows[1:])
    pre_max_c = np.maximum.accumulate(cols[1:])
    pre_min_c = np.minimum.accumulate(cols[1:])

    eps_lat = epsilon * N  # threshold in lattice units

    pairs = []
    for times in visit_times(loop).values():
        for i in range(len(times)):
            for j in range(i + 1, len(times)):
                s, t = times[i], times[j]

                # loop1 = loop[s:t+1] diameter — O(1) via sparse table
                d1 = max(_qmax(st_max_r, s, t) - _qmin(st_min_r, s, t),
                         _qmax(st_max_c, s, t) - _qmin(st_min_c, s, t))
                if d1 < eps_lat:
                    continue

                # loop2 = loop[t:] + loop[1:s+1] — O(1) via suffix/prefix
                mr2, nr2 = int(suf_max_r[t]), int(suf_min_r[t])
                mc2, nc2 = int(suf_max_c[t]), int(suf_min_c[t])
                if s >= 1:
                    mr2 = max(mr2, int(pre_max_r[s-1]))
                    nr2 = min(nr2, int(pre_min_r[s-1]))
                    mc2 = max(mc2, int(pre_max_c[s-1]))
                    nc2 = min(nc2, int(pre_min_c[s-1]))
                if max(mr2-nr2, mc2-nc2) >= eps_lat:
                    pairs.append((s, t))
    return pairs

def I_pairs(loop1, loop2):
    """
    Intersection pairs (s,t): loop1[s] = loop2[t].
    Drives the merge move.
    """
    vt2 = visit_times(loop2)
    pairs = []
    for s, p in enumerate(loop1):
        key = tuple(p)
        if key in vt2:
            for t in vt2[key]:
                pairs.append((s, t))
    return pairs

def split_loop(loop, s, t):
    """
    Split at self-intersection (s,t) into two loops of lengths t-s and T-(t-s).
    loop[s] = loop[t] required.
    """
    loop1 = loop[s:t + 1]               # s -> t
    loop2 = loop[t:] + loop[1:s + 1]    # t -> (wrap) -> s
    return loop1, loop2

def merge_loops(loop1, loop2, s, t):
    """
    Merge at intersection loop1[s] = loop2[t].
    Travels around loop1 from s, then around loop2 from t.
    Result has time-length T1 + T2.
    """
    part1 = loop1[s:] + loop1[1:s + 1]      # around loop1 starting at s
    part2 = loop2[t + 1:] + loop2[1:t + 1]  # around loop2 starting after t
    return part1 + part2

def rewiring_step(soup, epsilon, N):
    """
    One step of M^epsilon: sample uniformly from all valid split/merge events.

    Splits:  enumerate SI_pairs per loop         — O(T log T) per loop (sparse tables)
    Merges:  global site map over all loops      — O(total_length) instead of O(K^2 * avg_T)
    """
    # ── Splits ───────────────────────────────────────────────────────────────
    split_ops = []   # (loop_idx, s, t)
    for i, loop in enumerate(soup):
        for s, t in SI_pairs(loop, epsilon, N):
            split_ops.append((i, s, t))

    # ── Merges via global site map ────────────────────────────────────────────
    # One pass over all loops to build site -> [(loop_idx, time), ...]
    site_map = defaultdict(list)
    for i, loop in enumerate(soup):
        for t, p in enumerate(loop):
            site_map[tuple(p)].append((i, t))

    # Count merge ops per ordered (i<j) loop pair without storing them all
    merge_count = defaultdict(int)
    for entries in site_map.values():
        if len(entries) < 2:
            continue
        by_loop = defaultdict(list)
        for li, ti in entries:
            by_loop[li].append(ti)
        loop_ids = sorted(by_loop.keys())
        for a in range(len(loop_ids)):
            for b in range(a + 1, len(loop_ids)):
                ia, ib = loop_ids[a], loop_ids[b]
                merge_count[(ia, ib)] += len(by_loop[ia]) * len(by_loop[ib])

    n_splits = len(split_ops)
    n_merges = sum(merge_count.values())
    n_total  = n_splits + n_merges

    if n_total == 0:
        return soup, None, n_splits, n_merges

    # ── Sample uniformly ──────────────────────────────────────────────────────
    r = np.random.randint(n_total)

    if r < n_splits:
        i, s, t = split_ops[r]
        l1, l2 = split_loop(soup[i], s, t)
        new_soup = [l for k, l in enumerate(soup) if k != i] + [l1, l2]
        return new_soup, ('split', i, s, t), n_splits, n_merges

    else:
        r2 = r - n_splits
        cumulative = 0
        for (ia, ib), count in merge_count.items():
            if cumulative + count > r2:
                # Enumerate I_pairs only for this one chosen (ia, ib)
                pairs = I_pairs(soup[ia], soup[ib])
                s, t = pairs[r2 - cumulative]
                new_loop = merge_loops(soup[ia], soup[ib], s, t)
                new_soup = [l for k, l in enumerate(soup) if k != ia and k != ib] + [new_loop]
                return new_soup, ('merge', ia, ib, s, t), n_splits, n_merges
            cumulative += count

    return soup, None, n_splits, n_merges  # unreachable

# ── Visualisation ────────────────────────────────────────────────────────────

def draw_soup(ax, soup, N, title='', highlight=()):
    # highlight: list of (loop_idx, color) pairs
    ax.set_facecolor('white')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.axis('off')

    # Square boundary
    ax.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], '-', color='#2C3E50', linewidth=1.5, zorder=1)

    if not soup:
        ax.set_title(title, fontsize=11)
        return

    hl_indices = {idx for idx, _ in highlight}
    ax.set_title(title, fontsize=11)

    nan = [float('nan')]

    # All non-highlighted loops
    xs_bg, ys_bg = [], []
    for idx, loop in enumerate(soup):
        if idx in hl_indices:
            continue
        xs_bg += [p[1] / (N - 1) for p in loop] + nan
        ys_bg += [p[0] / (N - 1) for p in loop] + nan
    if xs_bg:
        ax.plot(xs_bg, ys_bg, '-', color='#1A1A1A', linewidth=0.45, alpha=0.88, zorder=2)

    # Each highlighted loop in its own color
    for idx, color in highlight:
        loop = soup[idx]
        xs = [p[1] / (N - 1) for p in loop]
        ys = [p[0] / (N - 1) for p in loop]
        ax.plot(xs, ys, '-', color=color, linewidth=2.2, alpha=1.0, zorder=4)

# ── Main ─────────────────────────────────────────────────────────────────────

import matplotlib.animation as animation

N = 1000
c = 0.5
epsilon = 0.001  # spatial diameter cutoff in [0,1]
n_steps = 200      # number of rewiring steps
fps = 1           # frames per second — 2 frames per step so 2s per operation

print("Generating loop soup...")
soup = LoopSoup(N, c, epsilon=epsilon)
print(f"Loops with diameter >= {epsilon}: {len(soup)}")

# Run the chain: 2 frames per step — BEFORE (targets highlighted) and AFTER (results highlighted)
print("Running rewiring chain...")
CA = '#2471A3'   # deep blue
CB = '#C0392B'   # crimson
CR = '#FFD700'  # gold — single-loop state (about-to-split / merge result)

# frame: (step, state, highlight, label, n_si, n_i, len_info)
frames = [(0, [loop[:] for loop in soup], [], '', 0, 0, '')]
current = [loop[:] for loop in soup]
for step in range(n_steps):
    prev = [loop[:] for loop in current]
    current, op, n_si, n_i = rewiring_step(current, epsilon, N)
    if op is None:
        print(f"  Step {step}: no valid ops, chain stuck")
        break

    if op[0] == 'split':
        _, i, s, t   = op
        diam_in      = f'{spatial_diameter(prev[i], N):.4f}'
        diam_a       = f'{spatial_diameter(current[-2], N):.4f}'
        diam_b       = f'{spatial_diameter(current[-1], N):.4f}'
        before_hl    = [(i, CR)]
        after_hl     = [(len(current)-2, CA), (len(current)-1, CB)]
        before_label = f'Split.  Yellow loop (diam {diam_in}) will split'
        after_label  = f'Split done.  Blue (diam {diam_a})  +  Red (diam {diam_b})'
        print(f"  Step {step+1}: split  diam {diam_in} → {diam_a} + {diam_b}  "
              f"(K {len(prev)}→{len(current)})  SI={n_si} I={n_i}")
    else:
        _, ia, ib, s, t = op
        diam_a       = f'{spatial_diameter(prev[ia], N):.4f}'
        diam_b       = f'{spatial_diameter(prev[ib], N):.4f}'
        diam_out     = f'{spatial_diameter(current[-1], N):.4f}'
        before_hl    = [(ia, CA), (ib, CB)]
        after_hl     = [(len(current)-1, CR)]
        before_label = f'Merge.  Blue (diam {diam_a})  +  Red (diam {diam_b})'
        after_label  = f'Merge done.  Yellow (diam {diam_out})'
        print(f"  Step {step+1}: merge  {diam_a} + {diam_b} → {diam_out}  "
              f"(K {len(prev)}→{len(current)})  SI={n_si} I={n_i}")
    if op[0] == 'split':
        len_info = f'diam {diam_in} → {diam_a} + {diam_b}'
    else:
        len_info = f'diam {diam_a} + {diam_b} → {diam_out}'
    frames.append((step,     prev,                           before_hl, before_label, n_si, n_i, len_info))
    frames.append((step + 1, [loop[:] for loop in current], after_hl,  after_label,  n_si, n_i, len_info))

# Build animation
print("Rendering frames...")
plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'cm'})
fig, ax = plt.subplots(figsize=(8, 8.5), facecolor='white')
fig.subplots_adjust(left=0.03, right=0.97, top=0.95, bottom=0.08)

# Info panel at the bottom — created once, updated each frame
info_text = fig.text(0.50, 0.02, '', ha='center', va='bottom', fontsize=10,
                     color='#333333')

def update(frame_idx):
    ax.clear()
    step, state, highlight, label, n_si, n_i, len_info = frames[frame_idx]
    title = f'Step {step}'
    if label:
        title += f'\n{label}'
    draw_soup(ax, state, N, title=title, highlight=highlight)
    info_text.set_text(
        f'$\\varepsilon$ = {epsilon}    c = {c}    N = {N}    '
        f'K = {len(state)}    SI = {n_si}    I = {n_i}'
        + (f'    {len_info}' if len_info else '')
    )

anim = animation.FuncAnimation(fig, update, frames=len(frames), interval=1000 // fps)

outpath = Path(__file__).resolve().parent / "loopsoup-rewiring-square.gif"
print("Saving GIF (this may take a moment)...")
anim.save(outpath, writer='pillow', fps=fps, dpi=120)
plt.close()
print(f"Saved to {outpath}")
