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

def LoopSoup(N, c):
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
                        loopsoup.append(path[now:next + 1])
                    now = next
                path = path[last + 1:]
    return loopsoup

def spatial_diameter(loop, N):
    return max(int(loop[:, 0].max()) - int(loop[:, 0].min()),
               int(loop[:, 1].max()) - int(loop[:, 1].min())) / N

# ── Sparse table helpers for O(1) range min/max ───────────────────────────────

def _build_st_max(arr):
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
    k = (r - l + 1).bit_length() - 1
    return int(max(st[k][l], st[k][r - (1 << k) + 1]))

def _qmin(st, l, r):
    k = (r - l + 1).bit_length() - 1
    return int(min(st[k][l], st[k][r - (1 << k) + 1]))

# ── Rewiring operations ───────────────────────────────────────────────────────

def visit_times(loop):
    # loop is numpy (T,2) int16; keys are plain Python int tuples for hashing
    visits = defaultdict(list)
    for t in range(len(loop)):
        visits[(int(loop[t, 0]), int(loop[t, 1]))].append(t)
    return dict(visits)

def SI_pairs(loop, epsilon, N):
    # numpy views — no list comprehension overhead
    rows = loop[:, 0].astype(np.int32)
    cols = loop[:, 1].astype(np.int32)

    st_max_r = _build_st_max(rows)
    st_min_r = _build_st_min(rows)
    st_max_c = _build_st_max(cols)
    st_min_c = _build_st_min(cols)

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
            for j in range(i + 1, len(times)):
                s, t = times[i], times[j]
                d1 = max(_qmax(st_max_r, s, t) - _qmin(st_min_r, s, t),
                         _qmax(st_max_c, s, t) - _qmin(st_min_c, s, t))
                if d1 < eps_lat:
                    continue
                mr2, nr2 = int(suf_max_r[t]), int(suf_min_r[t])
                mc2, nc2 = int(suf_max_c[t]), int(suf_min_c[t])
                if s >= 1:
                    mr2 = max(mr2, int(pre_max_r[s-1]))
                    nr2 = min(nr2, int(pre_min_r[s-1]))
                    mc2 = max(mc2, int(pre_max_c[s-1]))
                    nc2 = min(nc2, int(pre_min_c[s-1]))
                if max(mr2 - nr2, mc2 - nc2) >= eps_lat:
                    pairs.append((s, t))
    return pairs

def I_pairs(loop1, loop2):
    vt2 = visit_times(loop2)
    pairs = []
    for s in range(len(loop1)):
        key = (int(loop1[s, 0]), int(loop1[s, 1]))
        if key in vt2:
            for t in vt2[key]:
                pairs.append((s, t))
    return pairs

def split_loop(loop, s, t):
    return loop[s:t+1], np.concatenate([loop[t:], loop[1:s+1]])

def merge_loops(loop1, loop2, s, t):
    return np.concatenate([loop1[s:], loop1[1:s+1], loop2[t+1:], loop2[1:t+1]])

def rewiring_step(soup, epsilon, N):
    split_ops = []
    for i, loop in enumerate(soup):
        for s, t in SI_pairs(loop, epsilon, N):
            split_ops.append((i, s, t))

    site_map = defaultdict(list)
    for i, loop in enumerate(soup):
        for t in range(len(loop)):
            site_map[(int(loop[t, 0]), int(loop[t, 1]))].append((i, t))

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
                pairs = I_pairs(soup[ia], soup[ib])
                s, t = pairs[r2 - cumulative]
                new_loop = merge_loops(soup[ia], soup[ib], s, t)
                new_soup = [l for k, l in enumerate(soup) if k != ia and k != ib] + [new_loop]
                return new_soup, ('merge', ia, ib, s, t), n_splits, n_merges
            cumulative += count

    return soup, None, n_splits, n_merges

# ── Chain runner — top-level so multiprocessing (fork) can call it ────────────

def run_chain(args):
    init_soup, eps, N, n_steps, seed = args
    np.random.seed(seed)
    CA = '#2471A3'   # deep blue
    CB = '#C0392B'   # crimson
    CR =  '#FFD700'   # gold — single-loop state (about-to-split / merge result)
    cum_splits = 0
    cum_merges = 0
    # frame: (step, state, highlight, label, n_si, n_i, cum_splits, cum_merges)
    frames = [(0, list(init_soup), [], '', 0, 0, 0, 0)]
    current = list(init_soup)
    for step in range(n_steps):
        prev = list(current)
        current, op, n_si, n_i = rewiring_step(current, eps, N)
        if op is None:
            print(f"  ε={eps} step {step}: no valid ops, chain stuck", flush=True)
            break
        if op[0] == 'split':
            cum_splits += 1
            _, i, _, _      = op
            before_hl       = [(i, CR)]
            after_hl        = [(len(current)-2, CA), (len(current)-1, CB)]
            before_label    = 'Split.  Yellow loop will split'
            after_label     = 'Splitting done.  Result: Blue loop  +  Red loop'
        else:
            cum_merges += 1
            _, ia, ib, _, _ = op
            before_hl       = [(ia, CA), (ib, CB)]
            after_hl        = [(len(current)-1, CR)]
            before_label    = 'Merge.  Blue loop  +  Red loop about to merge'
            after_label     = 'Merging done.  Result: Yellow loop'
        total = cum_splits + cum_merges
        print(f"  ε={eps} step {total}: {op[0]}, {len(current)} loops  "
              f"splits={cum_splits}/{total}  merges={cum_merges}/{total}", flush=True)
        frames.append((step,     prev,          before_hl, before_label, n_si, n_i, cum_splits, cum_merges))
        frames.append((step + 1, list(current), after_hl,  after_label,  n_si, n_i, cum_splits, cum_merges))
    return frames

# ── Visualisation ─────────────────────────────────────────────────────────────

def draw_soup(ax, soup, N, title='', highlight=()):
    ax.set_facecolor('white')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.axis('off')

    # Square boundary
    ax.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], '-', color='#2C3E50', linewidth=1.5, zorder=1)

    if not soup:
        ax.set_title(title, fontsize=10)
        return

    hl_indices = {idx for idx, _ in highlight}
    ax.set_title(title, fontsize=10)
    nan_row = np.array([[np.nan, np.nan]])

    # All non-highlighted loops in one plot call via nan-row separator
    parts = []
    for idx, loop in enumerate(soup):
        if idx in hl_indices:
            continue
        parts.append(loop)
        parts.append(nan_row)
    if parts:
        all_pts = np.vstack(parts).astype(float)
        ax.plot(all_pts[:, 1] / (N - 1), all_pts[:, 0] / (N - 1),
                '-', color='#1A1A1A', linewidth=0.45, alpha=0.88, zorder=2)

    for idx, color in highlight:
        loop = soup[idx].astype(float)
        ax.plot(loop[:, 1] / (N - 1), loop[:, 0] / (N - 1),
                '-', color=color, linewidth=2.2, alpha=1.0, zorder=4)

# ── Main ──────────────────────────────────────────────────────────────────────

N        = 1000
c        = 1/2
epsilons = [0.05, 0.01, 0.005, 0.001]
n_steps  = 20
fps      = 2

print("Generating loop soup (this takes a moment)...")
raw_soup = LoopSoup(N, c)
print(f"Full soup: {len(raw_soup)} loops")

# Convert to numpy int16 arrays once — avoids list comprehensions in SI_pairs
full_soup = [np.array(loop, dtype=np.int16) for loop in raw_soup]

all_diams = [spatial_diameter(loop, N) for loop in full_soup]

init_soups = []
for eps in epsilons:
    filtered = [loop for loop, d in zip(full_soup, all_diams) if d >= eps]
    print(f"  epsilon={eps}: {len(filtered)} loops")
    init_soups.append(filtered)

# Run 4 chains in parallel using fork (no serialisation overhead for inputs)
print("Running rewiring chains (parallel)...")
seeds = np.random.randint(0, 2**31, size=len(epsilons)).tolist()
args_list = [(init_soups[ei], eps, N, n_steps, seeds[ei])
             for ei, eps in enumerate(epsilons)]

ctx = multiprocessing.get_context('fork')
with ctx.Pool(len(epsilons)) as pool:
    all_frames = pool.map(run_chain, args_list)
print("All chains done.")

# Build 2×2 animation
print("Rendering frames...")
plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'cm'})
fig, axes = plt.subplots(2, 2, figsize=(13, 13.5), facecolor='white')
fig.subplots_adjust(hspace=0.18, wspace=0.06, left=0.03, right=0.97, top=0.97, bottom=0.06)

# Bottom info strip — static, set once
fig.text(0.50, 0.01, f'c = {c}    N = {N}', ha='center', va='bottom',
         fontsize=11, color='#333333')

def update(frame_idx):
    for ax, eps, frames in zip(axes.flat, epsilons, all_frames):
        ax.clear()
        fidx = min(frame_idx, len(frames) - 1)
        step, state, highlight, label, n_si, n_i, cum_splits, cum_merges = frames[fidx]
        title = f'$\\varepsilon$ = {eps}    Step {step}'
        if label:
            title += f'\n{label}'
        draw_soup(ax, state, N, title=title, highlight=highlight)

n_frames = max(len(f) for f in all_frames)
anim = animation.FuncAnimation(fig, update, frames=n_frames, interval=1000 // fps)

outpath = Path(__file__).resolve().parent / "loopsoup-rewiring-compare-square.gif"
print("Saving GIF (may take a moment)...")
anim.save(outpath, writer='pillow', fps=fps, dpi=110)
plt.close()
print(f"Saved to {outpath}")
