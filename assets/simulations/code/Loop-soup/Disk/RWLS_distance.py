from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from collections import defaultdict

# ── Loop soup (Wilson's algorithm) ────────────────────────────────────────────

def spatial_diameter(loop, N):
    rows = [p[0] for p in loop]
    cols = [p[1] for p in loop]
    return max(max(rows) - min(rows), max(cols) - min(cols)) / N

def NextTip(N, tip):
    if tip[0] == 0:
        if tip[1] == 0:
            newtip = [1, 0] if np.random.randint(2) == 0 else [0, 1]
        elif tip[1] == N - 1:
            newtip = [1, N-1] if np.random.randint(2) == 0 else [0, N-2]
        else:
            d = np.random.randint(3)
            newtip = [1, tip[1]] if d == 0 else ([0, tip[1]+1] if d == 1 else [0, tip[1]-1])
    elif tip[0] == N - 1:
        if tip[1] == 0:
            newtip = [N-2, 0] if np.random.randint(2) == 0 else [N-1, 1]
        elif tip[1] == N - 1:
            newtip = [N-2, N-1] if np.random.randint(2) == 0 else [N-1, N-2]
        else:
            d = np.random.randint(3)
            newtip = [N-2, tip[1]] if d == 0 else ([N-1, tip[1]+1] if d == 1 else [N-1, tip[1]-1])
    elif tip[1] == 0:
        d = np.random.randint(3)
        newtip = [tip[0], 1] if d == 0 else ([tip[0]+1, 0] if d == 1 else [tip[0]-1, 0])
    elif tip[1] == N - 1:
        d = np.random.randint(3)
        newtip = [tip[0], N-2] if d == 0 else ([tip[0]+1, N-1] if d == 1 else [tip[0]-1, N-1])
    else:
        d = np.random.randint(4)
        if d == 0:   newtip = [tip[0]+1, tip[1]]
        elif d == 1: newtip = [tip[0]-1, tip[1]]
        elif d == 2: newtip = [tip[0], tip[1]+1]
        else:        newtip = [tip[0], tip[1]-1]
    return newtip

def LoopSoup(N, c, epsilon=0.0):
    loopsoup = []
    tree = np.zeros((N, N))
    for i in range(N):
        tree[i, 0] = tree[i, N-1] = tree[0, i] = tree[N-1, i] = 1
    for i in range(1, N-1):
        for j in range(1, N-1):
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
                label = True; rmax = 0
                last = len(path) - 1 - path[::-1].index(root)
                now = 0
                while now < last:
                    nxt = now + 1 + path[now+1:].index(root)
                    if (r := np.random.uniform()) > rmax:
                        rmax = r
                        label = (np.random.uniform() < c)
                    if label:
                        loop = path[now:nxt+1]
                        if spatial_diameter(loop, N) >= epsilon:
                            loopsoup.append(loop)
                    now = nxt
                path = path[last+1:]
    return loopsoup

# ── Rewiring ──────────────────────────────────────────────────────────────────

def _build_st_max(arr):
    n = len(arr); st = [arr]; j = 1
    while (1 << j) <= n:
        half = 1 << (j-1); end = n - (1 << j) + 1
        st.append(np.maximum(st[j-1][:end], st[j-1][half:half+end])); j += 1
    return st

def _build_st_min(arr):
    n = len(arr); st = [arr]; j = 1
    while (1 << j) <= n:
        half = 1 << (j-1); end = n - (1 << j) + 1
        st.append(np.minimum(st[j-1][:end], st[j-1][half:half+end])); j += 1
    return st

def _qmax(st, l, r):
    k = (r - l + 1).bit_length() - 1
    return int(max(st[k][l], st[k][r - (1 << k) + 1]))

def _qmin(st, l, r):
    k = (r - l + 1).bit_length() - 1
    return int(min(st[k][l], st[k][r - (1 << k) + 1]))

def visit_times(loop):
    visits = defaultdict(list)
    for t, p in enumerate(loop):
        visits[tuple(p)].append(t)
    return dict(visits)

def SI_pairs(loop, epsilon, N):
    rows = np.array([p[0] for p in loop], dtype=np.int32)
    cols = np.array([p[1] for p in loop], dtype=np.int32)
    st_max_r = _build_st_max(rows); st_min_r = _build_st_min(rows)
    st_max_c = _build_st_max(cols); st_min_c = _build_st_min(cols)
    suf_max_r = np.maximum.accumulate(rows[::-1])[::-1]
    suf_min_r = np.minimum.accumulate(rows[::-1])[::-1]
    suf_max_c = np.maximum.accumulate(cols[::-1])[::-1]
    suf_min_c = np.minimum.accumulate(cols[::-1])[::-1]
    pre_max_r = np.maximum.accumulate(rows[1:])
    pre_min_r = np.minimum.accumulate(rows[1:])
    pre_max_c = np.maximum.accumulate(cols[1:])
    pre_min_c = np.minimum.accumulate(cols[1:])
    eps_lat = epsilon * N
    pairs = []
    for times in visit_times(loop).values():
        for i in range(len(times)):
            for j in range(i+1, len(times)):
                s, t = times[i], times[j]
                d1 = max(_qmax(st_max_r, s, t) - _qmin(st_min_r, s, t),
                         _qmax(st_max_c, s, t) - _qmin(st_min_c, s, t))
                if d1 < eps_lat:
                    continue
                mr2, nr2 = int(suf_max_r[t]), int(suf_min_r[t])
                mc2, nc2 = int(suf_max_c[t]), int(suf_min_c[t])
                if s >= 1:
                    mr2 = max(mr2, int(pre_max_r[s-1])); nr2 = min(nr2, int(pre_min_r[s-1]))
                    mc2 = max(mc2, int(pre_max_c[s-1])); nc2 = min(nc2, int(pre_min_c[s-1]))
                if max(mr2-nr2, mc2-nc2) >= eps_lat:
                    pairs.append((s, t))
    return pairs

def I_pairs(loop1, loop2):
    vt2 = visit_times(loop2)
    pairs = []
    for s, p in enumerate(loop1):
        key = tuple(p)
        if key in vt2:
            for t in vt2[key]:
                pairs.append((s, t))
    return pairs

def split_loop(loop, s, t):
    return loop[s:t+1], loop[t:] + loop[1:s+1]

def merge_loops(loop1, loop2, s, t):
    return loop1[s:] + loop1[1:s+1] + loop2[t+1:] + loop2[1:t+1]

def rewiring_step(soup, epsilon, N):
    split_ops = []
    for i, loop in enumerate(soup):
        for s, t in SI_pairs(loop, epsilon, N):
            split_ops.append((i, s, t))
    site_map = defaultdict(list)
    for i, loop in enumerate(soup):
        for t, p in enumerate(loop):
            site_map[tuple(p)].append((i, t))
    merge_count = defaultdict(int)
    for entries in site_map.values():
        if len(entries) < 2:
            continue
        by_loop = defaultdict(list)
        for li, ti in entries:
            by_loop[li].append(ti)
        loop_ids = sorted(by_loop.keys())
        for a in range(len(loop_ids)):
            for b in range(a+1, len(loop_ids)):
                ia, ib = loop_ids[a], loop_ids[b]
                merge_count[(ia, ib)] += len(by_loop[ia]) * len(by_loop[ib])
    n_splits = len(split_ops)
    n_merges = sum(merge_count.values())
    n_total  = n_splits + n_merges
    if n_total == 0:
        return soup, None
    r = np.random.randint(n_total)
    if r < n_splits:
        i, s, t = split_ops[r]
        l1, l2 = split_loop(soup[i], s, t)
        return [l for k, l in enumerate(soup) if k != i] + [l1, l2], ('split', i, s, t)
    else:
        r2 = r - n_splits; cumulative = 0
        for (ia, ib), count in merge_count.items():
            if cumulative + count > r2:
                pairs = I_pairs(soup[ia], soup[ib])
                s, t = pairs[r2 - cumulative]
                new_loop = merge_loops(soup[ia], soup[ib], s, t)
                return [l for k, l in enumerate(soup) if k != ia and k != ib] + [new_loop], ('merge', ia, ib, s, t)
            cumulative += count
    return soup, None

# ── Colors ────────────────────────────────────────────────────────────────────

CA = '#2471A3'   # blue  — solid   (split result A / merge target A)
CB = '#C0392B'   # red   — dashed  (split result B / merge target B)
CR = '#D4AC0D'   # gold           (split target   / merge result)

# ── Distance statistics ───────────────────────────────────────────────────────

def loop_dist(loop, z, N):
    """
    Distance from loop (as compact set of grid points) to z ∈ [0,1]².
    d(L, z) = min_{p ∈ L} ||p/(N-1) - z||_2
    """
    rows = np.array([p[0] for p in loop], dtype=float) / (N - 1)
    cols = np.array([p[1] for p in loop], dtype=float) / (N - 1)
    return float(np.sqrt((rows - z[0])**2 + (cols - z[1])**2).min())

def compute_D(soup, z, N):
    """
    D = max_{L ∈ soup} d(L, z)
    Returns (D_val, D_idx) — the maximum distance and which loop achieves it.
    """
    if not soup:
        return 0.0, -1
    dists = [loop_dist(loop, z, N) for loop in soup]
    idx = int(np.argmax(dists))
    return dists[idx], idx

# ── Visualisation ─────────────────────────────────────────────────────────────

def draw_soup_with_z(ax, soup, N, z, D_val, D_idx, rewiring_hl=()):
    """
    Draw the soup with:
      - background loops in dark grey
      - most-distant loop (D_idx) in orange glow (unless it's also a rewiring hl)
      - witness segment from z to nearest point on D loop
      - dashed circle of radius D_val around z
      - rewiring highlights (gold/blue/red) drawn on top
      - small black dot at z
    rewiring_hl: list of (loop_idx, color, linestyle) tuples
    """
    ax.set_facecolor('white')
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_aspect('equal'); ax.axis('off')
    ax.add_patch(plt.Circle((0.5, 0.5), 0.5, color='#2C3E50', fill=False, linewidth=1.5))

    clip = plt.Circle((0.5, 0.5), 0.5, transform=ax.transData)
    nan  = [float('nan')]

    hl_indices = {h[0] for h in rewiring_hl}

    # Background loops (skip D_idx and rewiring hl indices — they get special treatment)
    xs_bg, ys_bg = [], []
    for idx, loop in enumerate(soup):
        if idx == D_idx or idx in hl_indices:
            continue
        xs_bg += [p[1] / (N-1) for p in loop] + nan
        ys_bg += [p[0] / (N-1) for p in loop] + nan
    if xs_bg:
        line, = ax.plot(xs_bg, ys_bg, '-', color='#1A1A1A', linewidth=0.45, alpha=0.88, zorder=2)
        line.set_clip_path(clip)

    # Most-distant loop: white glow + bright orange
    # (skip if it also appears in rewiring_hl — the rewiring color takes precedence)
    if D_idx >= 0 and D_idx < len(soup) and D_idx not in hl_indices:
        loop = soup[D_idx]
        xs = [p[1] / (N-1) for p in loop]
        ys = [p[0] / (N-1) for p in loop]
        glow, = ax.plot(xs, ys, '-', color='white', linewidth=9, alpha=0.7, zorder=4)
        glow.set_clip_path(clip)
        line, = ax.plot(xs, ys, '-', color='#E67E22', linewidth=4.0, alpha=1.0, zorder=5)
        line.set_clip_path(clip)

    # Witness segment from z to nearest point on D loop (always shown)
    if D_idx >= 0 and D_idx < len(soup):
        loop = soup[D_idx]
        rows_s = np.array([p[0] for p in loop], dtype=float) / (N - 1)
        cols_s = np.array([p[1] for p in loop], dtype=float) / (N - 1)
        dists_s = np.sqrt((rows_s - z[0])**2 + (cols_s - z[1])**2)
        nearest = int(np.argmin(dists_s))
        nx, ny = cols_s[nearest], rows_s[nearest]
        ax.plot([z[1], nx], [z[0], ny], '-', color='#E67E22',
                linewidth=1.8, alpha=0.9, zorder=6)
        ax.plot(nx, ny, 'o', color='#E67E22', markersize=7, zorder=7,
                markeredgecolor='white', markeredgewidth=1.0)

    # Rewiring highlights drawn on top of everything else
    for h in rewiring_hl:
        idx, color = h[0], h[1]
        ls = h[2] if len(h) > 2 else '-'
        lw = 3.0 if ls == '-' else 2.5
        if idx < 0 or idx >= len(soup):
            continue
        loop = soup[idx]
        xs = [p[1] / (N-1) for p in loop]
        ys = [p[0] / (N-1) for p in loop]
        glow, = ax.plot(xs, ys, '-', color='white', linewidth=9, alpha=0.7, zorder=8)
        glow.set_clip_path(clip)
        line, = ax.plot(xs, ys, ls, color=color, linewidth=lw, alpha=0.95, zorder=9)
        line.set_clip_path(clip)

    # Dashed circle of radius D_val around z
    theta = np.linspace(0, 2 * np.pi, 300)
    cx = z[1] + D_val * np.cos(theta)
    cy = z[0] + D_val * np.sin(theta)
    line, = ax.plot(cx, cy, '--', color='#E67E22', linewidth=1.2, alpha=0.5, zorder=3)
    line.set_clip_path(clip)

    # Mark z — small black dot
    ax.plot(z[1], z[0], 'o', color='black', markersize=4, zorder=11)

# ── Main ──────────────────────────────────────────────────────────────────────

N       = 1000
c       = 0.5
epsilon = 0.1
n_steps = 500
fps     = 2

# z ∈ [0,1]²: marked point in the interior of the disk (row, col)
z = (0.5, 0.5)   # centre — change to e.g. (0.3, 0.4) to explore off-centre

print("Generating loop soup...")
soup = LoopSoup(N, c, epsilon=epsilon)
print(f"Loops with diameter ≥ {epsilon}: {len(soup)}")

# ── Run chain, recording D and rewiring highlights at every step ───────────────

current  = [loop[:] for loop in soup]
D_cur, D_idx_cur = compute_D(current, z, N)

# frame tuple: (step, soup_state, D_val, D_idx, rewiring_hl, phase, d_count)
#   step        — rewiring step number (0 = initial)
#   soup_state  — snapshot of the soup
#   D_val       — D value for this state
#   D_idx       — index of most-distant loop in soup_state
#   rewiring_hl — list of (idx, color, linestyle) for highlighted loops
#   phase       — '' | 'before split' | 'after split' | 'before merge' | 'after merge'
#   d_count     — how many entries of D_series / step_series to display on the D(t) plot

frames      = [(0, current, D_cur, D_idx_cur, [], '', 1)]
D_series    = [D_cur]
step_series = [0]
op_types    = ['']   # op_types[k] = '' | 'split' | 'merge', parallel to D_series

print("Running rewiring chain...")
for step in range(n_steps):
    prev_state   = [loop[:] for loop in current]
    D_prev       = D_cur
    D_idx_prev   = D_idx_cur

    new_soup, op = rewiring_step(current, epsilon, N)
    if op is None:
        print(f"  Step {step}: no valid ops, chain stuck")
        break

    D_new, D_idx_new = compute_D(new_soup, z, N)
    n_new = len(new_soup)

    if op[0] == 'split':
        _, i, s, t = op
        before_hl = [(i, CR, '-')]
        after_hl  = [(n_new - 2, CA, '-'), (n_new - 1, CB, '--')]
        op_type   = 'split'
    else:
        _, ia, ib, s, t = op
        before_hl = [(ia, CA, '-'), (ib, CB, '--')]
        after_hl  = [(n_new - 1, CR, '-')]
        op_type   = 'merge'

    d_before = len(D_series)     # D values visible before the op (on D curve)
    D_series.append(D_new)
    step_series.append(step + 1)
    op_types.append(op_type)
    d_after = len(D_series)      # D values visible after the op

    # BEFORE frame: prev state highlighted to show what will change; D curve up to now
    frames.append((step, prev_state, D_prev, D_idx_prev,
                   before_hl, f'before {op_type}', d_before))
    # AFTER frame: new state with result highlighted; D curve extended by one point
    frames.append((step + 1, [loop[:] for loop in new_soup], D_new, D_idx_new,
                   after_hl, f'after {op_type}', d_after))

    current      = new_soup
    D_cur        = D_new
    D_idx_cur    = D_idx_new
    print(f"  Step {step+1}: {op_type}, K={len(current)}, D={D_new:.4f}", flush=True)

# ── Build animation ───────────────────────────────────────────────────────────

print("Rendering animation...")
plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'cm'})

fig, (ax_soup, ax_D) = plt.subplots(
    1, 2, figsize=(14, 7), facecolor='white',
    gridspec_kw={'width_ratios': [1, 1]}
)
fig.suptitle(
    r'Loop Soup Rewiring  $M^\varepsilon$'
    f',   $c = {c}$,   $N = {N}$,   $\\varepsilon = {epsilon}$'
    f',   $z = ({z[0]}, {z[1]})$',
    fontsize=13, fontweight='bold', y=0.98
)
fig.subplots_adjust(top=0.90, wspace=0.35)

D_max_global = max(D_series) * 1.1
D_min_global = 0.0

def update(frame_idx):
    step, state, D_val, D_idx, rewiring_hl, phase, d_count = frames[frame_idx]

    # ── Left panel: soup ──────────────────────────────────────────────────────
    ax_soup.clear()
    draw_soup_with_z(ax_soup, state, N, z, D_val, D_idx, rewiring_hl=rewiring_hl)

    if 'before' in phase:
        op_name   = phase.replace('before ', '')
        phase_str = f'  [about to {op_name}]'
    elif 'after' in phase:
        op_name   = phase.replace('after ', '')
        phase_str = f'  [{op_name}]'
    else:
        phase_str = ''

    ax_soup.set_title(
        f'step {step}{phase_str}   ·   $K = {len(state)}$ loops'
        f'\n$D = {D_val:.4f}$'
        r'   (orange = most distant,  dashed circle = radius $D$)',
        fontsize=10
    )

    # ── Right panel: D(t) curve ───────────────────────────────────────────────
    ax_D.clear()
    ax_D.set_facecolor('white')

    xs = step_series[:d_count]
    ys = D_series[:d_count]
    ax_D.plot(xs, ys, '-', color='#2471A3', linewidth=1.8, zorder=3)
    ax_D.plot(xs[-1], ys[-1], 'o', color='#2471A3', markersize=6, zorder=4)

    # Colour background intervals by operation type
    for k in range(1, d_count):
        col = '#AED6F1' if op_types[k] == 'split' else '#F1948A'
        ax_D.axvspan(step_series[k-1], step_series[k], alpha=0.18,
                     color=col, linewidth=0)

    ax_D.set_xlim(0, max(step_series[-1], 1))
    ax_D.set_ylim(D_min_global, D_max_global)
    ax_D.set_xlabel('rewiring step', fontsize=11)
    ax_D.set_ylabel('$D(t)$', fontsize=13)
    ax_D.set_title(
        r'$D(t) = \max_{L \in \mathrm{soup}} \, d(L,\, z)$' + '\n'
        r'{\footnotesize blue bg = split,  red bg = merge}',
        fontsize=11
    )
    ax_D.grid(True, alpha=0.25, linewidth=0.6)
    for sp in ax_D.spines.values():
        sp.set_linewidth(0.6)

anim = animation.FuncAnimation(fig, update, frames=len(frames), interval=1000 // fps)

outpath = Path(__file__).resolve().parent / 'loopsoup-distance.gif'
print("Saving GIF...")
anim.save(outpath, writer='pillow', fps=fps, dpi=110)
plt.close()
print(f"Saved → {outpath}")
