from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from collections import defaultdict
import multiprocessing

# ── Random walk on N×N grid ───────────────────────────────────────────────────

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

def spatial_diameter(loop, N):
    rows = [p[0] for p in loop]
    cols = [p[1] for p in loop]
    return max(max(rows) - min(rows), max(cols) - min(cols)) / N

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

# ── Rewiring ───────────────────────────────────────────────────────────────────

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

# ── Distance ───────────────────────────────────────────────────────────────────

def loop_dist(loop, z, N):
    rows = np.array([p[0] for p in loop], dtype=float) / (N - 1)
    cols = np.array([p[1] for p in loop], dtype=float) / (N - 1)
    return float(np.sqrt((rows - z[0])**2 + (cols - z[1])**2).min())

def compute_D(soup, z, N):
    if not soup:
        return 0.0, -1
    dists = [loop_dist(loop, z, N) for loop in soup]
    idx = int(np.argmax(dists))
    return dists[idx], idx

# ── Colors ─────────────────────────────────────────────────────────────────────

CA = '#2471A3'
CB = '#C0392B'
CR = '#FFD700'

# ── Chain runner — top-level for multiprocessing fork ─────────────────────────

def run_chain(args):
    init_soup, eps, N, n_steps, seed, z = args
    np.random.seed(seed)
    current = [loop[:] for loop in init_soup]
    D_cur, D_idx_cur = compute_D(current, z, N)
    # frame: (step, state, D_val, D_idx, rewiring_hl, op_type)
    # Only record frames where D changes
    frames = []
    d_changes = 0
    for step in range(n_steps):
        D_prev = D_cur
        new_soup, op = rewiring_step(current, eps, N)
        if op is None:
            print(f"  ε={eps} step {step}: no valid ops, chain stuck", flush=True)
            break
        D_new, D_idx_new = compute_D(new_soup, z, N)
        n_new = len(new_soup)
        if op[0] == 'split':
            _, i, s, t = op
            after_hl = [(n_new-2, CA, '-'), (n_new-1, CB, '--')]
            op_type  = 'split'
        else:
            _, ia, ib, s, t = op
            after_hl = [(n_new-1, CR, '-')]
            op_type  = 'merge'
        if abs(D_new - D_prev) > 1e-10:
            d_changes += 1
            frames.append((step + 1, [loop[:] for loop in new_soup],
                           D_new, D_idx_new, after_hl, op_type))
            print(f"  ε={eps} step {step+1}: D changed {D_prev:.4f} -> {D_new:.4f} ({op_type})", flush=True)
        current   = new_soup
        D_cur     = D_new
        D_idx_cur = D_idx_new
    print(f"  ε={eps}: {d_changes} D-changes out of {n_steps} steps", flush=True)
    return frames

# ── Visualisation ──────────────────────────────────────────────────────────────

def draw_soup_with_z(ax, soup, N, z, D_val, D_idx, rewiring_hl=()):
    ax.set_facecolor('white')
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_aspect('equal'); ax.axis('off')
    ax.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], '-', color='#2C3E50', linewidth=1.5, zorder=1)

    nan = [float('nan')]
    hl_indices = {h[0] for h in rewiring_hl}

    # Background loops
    xs_bg, ys_bg = [], []
    for idx, loop in enumerate(soup):
        if idx == D_idx or idx in hl_indices:
            continue
        xs_bg += [p[1] / (N-1) for p in loop] + nan
        ys_bg += [p[0] / (N-1) for p in loop] + nan
    if xs_bg:
        ax.plot(xs_bg, ys_bg, '-', color='#1A1A1A', linewidth=0.45, alpha=0.88, zorder=2)

    # Most-distant loop in green (unless also a rewiring highlight)
    if 0 <= D_idx < len(soup) and D_idx not in hl_indices:
        loop = soup[D_idx]
        xs = [p[1] / (N-1) for p in loop]
        ys = [p[0] / (N-1) for p in loop]
        ax.plot(xs, ys, '-', color='white',    linewidth=9,   alpha=0.7,  zorder=4)
        ax.plot(xs, ys, '-', color='#27AE60',  linewidth=4.0, alpha=1.0,  zorder=5)

    # Witness segment from z to nearest point on D loop
    if 0 <= D_idx < len(soup):
        loop = soup[D_idx]
        rows_s = np.array([p[0] for p in loop], dtype=float) / (N - 1)
        cols_s = np.array([p[1] for p in loop], dtype=float) / (N - 1)
        nearest = int(np.argmin(np.sqrt((rows_s - z[0])**2 + (cols_s - z[1])**2)))
        nx, ny = cols_s[nearest], rows_s[nearest]
        ax.plot([z[1], nx], [z[0], ny], '-', color='#27AE60', linewidth=1.8, alpha=0.9, zorder=6)
        ax.plot(nx, ny, 'o', color='#27AE60', markersize=7, zorder=7,
                markeredgecolor='white', markeredgewidth=1.0)

    # Rewiring highlights on top
    for h in rewiring_hl:
        idx, color = h[0], h[1]
        ls = h[2] if len(h) > 2 else '-'
        lw = 3.0 if ls == '-' else 2.5
        if idx < 0 or idx >= len(soup):
            continue
        loop = soup[idx]
        xs = [p[1] / (N-1) for p in loop]
        ys = [p[0] / (N-1) for p in loop]
        ax.plot(xs, ys, '-',  color='white', linewidth=9,  alpha=0.7,  zorder=8)
        ax.plot(xs, ys,  ls, color=color,   linewidth=lw, alpha=0.95, zorder=9)

    # Dashed circle of radius D_val around z
    theta = np.linspace(0, 2 * np.pi, 300)
    ax.plot(z[1] + D_val * np.cos(theta), z[0] + D_val * np.sin(theta),
            '--', color='#27AE60', linewidth=1.2, alpha=0.5, zorder=3)

    # Mark z
    ax.plot(z[1], z[0], 'o', color='black', markersize=4, zorder=11)

# ── Main ────────────────────────────────────────────────────────────────────────

N        = 1000
c        = 0.5
epsilons = [0.15, 0.08, 0.03, 0.01]
n_steps  = 50
fps      = 2
z        = (0.5, 0.5)

print("Generating loop soup...")
raw_soup = LoopSoup(N, c)
print(f"Full soup: {len(raw_soup)} loops")

diams = [spatial_diameter(loop, N) for loop in raw_soup]

init_soups = []
for eps in epsilons:
    filtered = [loop for loop, d in zip(raw_soup, diams) if d >= eps]
    print(f"  epsilon={eps}: {len(filtered)} loops")
    init_soups.append(filtered)

print("Running rewiring chains (parallel)...")
seeds = np.random.randint(0, 2**31, size=len(epsilons)).tolist()
args_list = [(init_soups[ei], eps, N, n_steps, seeds[ei], z)
             for ei, eps in enumerate(epsilons)]

ctx = multiprocessing.get_context('fork')
with ctx.Pool(len(epsilons)) as pool:
    all_frames = pool.map(run_chain, args_list)
print("All chains done.")

print("Rendering frames...")
plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'cm'})
fig, axes = plt.subplots(2, 2, figsize=(13, 13), facecolor='white')
fig.subplots_adjust(hspace=0.12, wspace=0.06, left=0.03, right=0.97, top=0.97, bottom=0.04)
fig.text(0.50, 0.01, f'c = {c}    N = {N}    z = {z}',
         ha='center', va='bottom', fontsize=11, color='#333333')

def update(frame_idx):
    for ax, eps, frames in zip(axes.flat, epsilons, all_frames):
        ax.clear()
        if not frames:
            ax.set_facecolor('white'); ax.axis('off')
            ax.set_title(f'$\\varepsilon$ = {eps}    no D-changes in {n_steps} steps', fontsize=10)
            continue
        fidx = min(frame_idx, len(frames) - 1)
        step, state, D_val, D_idx, rewiring_hl, op_type = frames[fidx]
        title = (f'$\\varepsilon$ = {eps}    Step {step}    D = {D_val:.3f}'
                 f'\nD changed via {op_type}  ({fidx+1}/{len(frames)} events)')
        draw_soup_with_z(ax, state, N, z, D_val, D_idx, rewiring_hl=rewiring_hl)
        ax.set_title(title, fontsize=10)

n_frames = max((len(f) for f in all_frames if f), default=1)
anim = animation.FuncAnimation(fig, update, frames=n_frames, interval=1000 // fps)

outpath = Path(__file__).resolve().parent / 'loopsoup-distance-compare-square.gif'
print("Saving GIF...")
anim.save(outpath, writer='pillow', fps=fps, dpi=110)
plt.close()
print(f"Saved → {outpath}")
