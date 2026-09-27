"""
RWLS_single_loop.py
-------------------
M^epsilon rewiring chain started from a SINGLE loop.

Initial state: run the full Wilson loop soup, keep only loops with spatial
diameter >= MIN_DIAM, then pick one uniformly at random.
The chain's transition kernel is identical to the full M^epsilon chain.
"""

from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from collections import defaultdict

# ── Parameters ────────────────────────────────────────────────────────────────
N        = 2000    # grid size
C        = 0.5     # loop-labelling probability (same as full chain)
MIN_DIAM = 0.1     # keep loops with spatial diameter >= this
EPSILON  = 0.0005   # rewiring diameter cutoff
N_STEPS  = 120
FPS      = 1
SEED     = 0


# ── Wilson loop soup ──────────────────────────────────────────────────────────

def NextTip(N, tip):
    if tip[0] == 0:
        if tip[1] == 0:
            newtip = [1,0] if np.random.randint(2)==0 else [0,1]
        elif tip[1] == N-1:
            newtip = [1,N-1] if np.random.randint(2)==0 else [0,N-2]
        else:
            d = np.random.randint(3)
            newtip = [[1,tip[1]],[0,tip[1]+1],[0,tip[1]-1]][d]
    elif tip[0] == N-1:
        if tip[1] == 0:
            newtip = [N-2,0] if np.random.randint(2)==0 else [N-1,1]
        elif tip[1] == N-1:
            newtip = [N-2,N-1] if np.random.randint(2)==0 else [N-1,N-2]
        else:
            d = np.random.randint(3)
            newtip = [[N-2,tip[1]],[N-1,tip[1]+1],[N-1,tip[1]-1]][d]
    elif tip[1] == 0:
        d = np.random.randint(3)
        newtip = [[tip[0],1],[tip[0]+1,0],[tip[0]-1,0]][d]
    elif tip[1] == N-1:
        d = np.random.randint(3)
        newtip = [[tip[0],N-2],[tip[0]+1,N-1],[tip[0]-1,N-1]][d]
    else:
        d = np.random.randint(4)
        newtip = [[tip[0]+1,tip[1]],[tip[0]-1,tip[1]],[tip[0],tip[1]+1],[tip[0],tip[1]-1]][d]
    return newtip

def spatial_diameter(loop, N):
    rows = [p[0] for p in loop]
    cols = [p[1] for p in loop]
    return max(max(rows)-min(rows), max(cols)-min(cols)) / N

def LoopSoup(N, c, min_diam=0.0):
    np.random.seed(SEED)
    loopsoup = []
    tree = np.zeros((N, N))
    for i in range(N):
        tree[i,0] = tree[i,N-1] = tree[0,i] = tree[N-1,i] = 1
    for i in range(1, N-1):
        for j in range(1, N-1):
            if tree[i,j] == 1:
                continue
            tip = [i, j]
            path = [tip]
            while tree[tip[0]][tip[1]] == 0:
                tip = NextTip(N, tip)
                path.append(tip)
            while path:
                root = path[0]
                tree[root[0]][root[1]] = 1
                label = True;  rmax = 0
                last = len(path)-1 - path[::-1].index(root)
                now = 0
                while now < last:
                    nxt = now + 1 + path[now+1:].index(root)
                    if (r := np.random.uniform()) > rmax:
                        rmax = r
                        label = (np.random.uniform() < c)
                    if label:
                        loop = path[now:nxt+1]
                        if spatial_diameter(loop, N) >= min_diam:
                            loopsoup.append([tuple(p) for p in loop])
                    now = nxt
                path = path[last+1:]
    return loopsoup

def make_initial_loop(N):
    """Run full Wilson soup, filter by MIN_DIAM, pick one uniformly at random."""
    print("  Running Wilson loop soup...")
    soup = LoopSoup(N, C, min_diam=MIN_DIAM)
    print(f"  {len(soup)} loops with diameter >= {MIN_DIAM}")
    loop = max(soup, key=lambda l: spatial_diameter(l, N))
    print(f"  Picked largest loop: length={len(loop)}, diameter={spatial_diameter(loop, N):.4f}")
    return loop


# ── Sparse-table helpers (O(n log n) build, O(1) range max/min) ───────────────

def _build_st_max(arr):
    n, st, j = len(arr), [arr], 1
    while (1 << j) <= n:
        half = 1 << (j-1)
        end  = n - (1 << j) + 1
        st.append(np.maximum(st[j-1][:end], st[j-1][half:half+end]))
        j += 1
    return st

def _build_st_min(arr):
    n, st, j = len(arr), [arr], 1
    while (1 << j) <= n:
        half = 1 << (j-1)
        end  = n - (1 << j) + 1
        st.append(np.minimum(st[j-1][:end], st[j-1][half:half+end]))
        j += 1
    return st

def _qmax(st, l, r):
    k = (r - l + 1).bit_length() - 1
    return int(max(st[k][l], st[k][r - (1 << k) + 1]))

def _qmin(st, l, r):
    k = (r - l + 1).bit_length() - 1
    return int(min(st[k][l], st[k][r - (1 << k) + 1]))


# ── Core rewiring helpers ─────────────────────────────────────────────────────

def visit_times(loop):
    visits = defaultdict(list)
    for t, p in enumerate(loop):
        visits[tuple(p)].append(t)
    return dict(visits)

def SI_pairs(loop, epsilon, N):
    """
    Self-intersection pairs (s, t) valid for a split.
    Both resulting sub-loops must have spatial diameter >= epsilon.
    """
    rows = np.array([p[0] for p in loop], dtype=np.int32)
    cols = np.array([p[1] for p in loop], dtype=np.int32)

    st_max_r = _build_st_max(rows);  st_min_r = _build_st_min(rows)
    st_max_c = _build_st_max(cols);  st_min_c = _build_st_min(cols)

    suf_max_r = np.maximum.accumulate(rows[::-1])[::-1]
    suf_min_r = np.minimum.accumulate(rows[::-1])[::-1]
    suf_max_c = np.maximum.accumulate(cols[::-1])[::-1]
    suf_min_c = np.minimum.accumulate(cols[::-1])[::-1]

    pre_max_r = np.maximum.accumulate(rows[1:])
    pre_min_r = np.minimum.accumulate(rows[1:])
    pre_max_c = np.maximum.accumulate(cols[1:])
    pre_min_c = np.minimum.accumulate(cols[1:])

    eps_lat = epsilon * N
    pairs   = []

    for times in visit_times(loop).values():
        for i in range(len(times)):
            for j in range(i + 1, len(times)):
                s, t = times[i], times[j]

                d1 = max(_qmax(st_max_r, s, t) - _qmin(st_min_r, s, t),
                         _qmax(st_max_c, s, t) - _qmin(st_min_c, s, t))
                if d1 < eps_lat:
                    continue

                mr2 = int(suf_max_r[t]);  nr2 = int(suf_min_r[t])
                mc2 = int(suf_max_c[t]);  nc2 = int(suf_min_c[t])
                if s >= 1:
                    mr2 = max(mr2, int(pre_max_r[s-1]));  nr2 = min(nr2, int(pre_min_r[s-1]))
                    mc2 = max(mc2, int(pre_max_c[s-1]));  nc2 = min(nc2, int(pre_min_c[s-1]))
                if max(mr2 - nr2, mc2 - nc2) >= eps_lat:
                    pairs.append((s, t))
    return pairs

def I_pairs(loop1, loop2):
    vt2   = visit_times(loop2)
    pairs = []
    for s, p in enumerate(loop1):
        for t in vt2.get(tuple(p), []):
            pairs.append((s, t))
    return pairs

def split_loop(loop, s, t):
    return loop[s:t+1], loop[t:] + loop[1:s+1]

def merge_loops(loop1, loop2, s, t):
    part1 = loop1[s:] + loop1[1:s+1]
    part2 = loop2[t+1:] + loop2[1:t+1]
    return part1 + part2

def rewiring_step(soup, epsilon, N):
    """One step of M^epsilon — identical transition kernel to the full chain."""
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
        return soup, None, 0, 0

    r = np.random.randint(n_total)

    if r < n_splits:
        i, s, t = split_ops[r]
        l1, l2  = split_loop(soup[i], s, t)
        new_soup = [l for k, l in enumerate(soup) if k != i] + [l1, l2]
        return new_soup, ('split', i, s, t), n_splits, n_merges
    else:
        r2 = r - n_splits
        cumulative = 0
        for (ia, ib), count in merge_count.items():
            if cumulative + count > r2:
                pairs    = I_pairs(soup[ia], soup[ib])
                s, t     = pairs[r2 - cumulative]
                new_loop = merge_loops(soup[ia], soup[ib], s, t)
                new_soup = [l for k, l in enumerate(soup) if k != ia and k != ib] + [new_loop]
                return new_soup, ('merge', ia, ib, s, t), n_splits, n_merges
            cumulative += count

    return soup, None, n_splits, n_merges


# ── Visualisation ─────────────────────────────────────────────────────────────

def draw_soup(ax, soup, N, title='', highlight=()):
    ax.set_facecolor('white')
    ax.set_xlim(0, 1);  ax.set_ylim(0, 1)
    ax.set_aspect('equal');  ax.axis('off')
    ax.plot([0,1,1,0,0],[0,0,1,1,0],'-', color='#2C3E50', linewidth=1.5, zorder=1)
    if not soup:
        ax.set_title(title, fontsize=11);  return

    hl_indices = {idx for idx, _ in highlight}
    ax.set_title(title, fontsize=11)
    nan = [float('nan')]

    xs_bg, ys_bg = [], []
    for idx, loop in enumerate(soup):
        if idx in hl_indices:
            continue
        xs_bg += [p[1]/(N-1) for p in loop] + nan
        ys_bg += [p[0]/(N-1) for p in loop] + nan
    if xs_bg:
        ax.plot(xs_bg, ys_bg, '-', color='#1A1A1A', linewidth=1.2, alpha=0.88, zorder=2)

    for idx, color in highlight:
        loop = soup[idx]
        ax.plot([p[1]/(N-1) for p in loop], [p[0]/(N-1) for p in loop],
                '-', color=color, linewidth=2.5, alpha=1.0, zorder=4)


# ── Main ──────────────────────────────────────────────────────────────────────

CA = '#2471A3'   # blue
CB = '#C0392B'   # red
CR = '#E67E22'   # orange

print(f"Generating initial loop (Wilson soup filtered at MIN_DIAM={MIN_DIAM}, picking one at random)...")
initial = make_initial_loop(N)
n_si    = len(SI_pairs(initial, EPSILON, N))
print(f"  Length = {len(initial)} steps, "
      f"diameter = {spatial_diameter(initial, N):.4f}, "
      f"valid splits = {n_si}")

soup    = [initial]
frames  = [(0, [loop[:] for loop in soup], [], '', 0, 0, '')]
current = [loop[:] for loop in soup]

print(f"Running {N_STEPS}-step rewiring chain...")
for step in range(N_STEPS):
    prev    = [loop[:] for loop in current]
    current, op, n_si, n_i = rewiring_step(current, EPSILON, N)

    if op is None:
        print(f"  Step {step}: no valid ops — chain stuck (K={len(prev)})")
        break

    if op[0] == 'split':
        _, i, s, t = op
        din  = f'{spatial_diameter(prev[i], N):.4f}'
        da   = f'{spatial_diameter(current[-2], N):.4f}'
        db   = f'{spatial_diameter(current[-1], N):.4f}'
        bhl  = [(i, CR)];  ahl = [(len(current)-2, CA), (len(current)-1, CB)]
        blbl = f'Split  — orange (diam {din}) splits'
        albl = f'Split done  — blue (diam {da}) + red (diam {db})'
        linfo = f'diam {din} → {da}+{db}'
        print(f"  Step {step+1}: split  {din} → {da}+{db}  "
              f"K {len(prev)}→{len(current)}  SI={n_si} I={n_i}")
    else:
        _, ia, ib, s, t = op
        da   = f'{spatial_diameter(prev[ia], N):.4f}'
        db   = f'{spatial_diameter(prev[ib], N):.4f}'
        dout = f'{spatial_diameter(current[-1], N):.4f}'
        bhl  = [(ia, CA), (ib, CB)];  ahl = [(len(current)-1, CR)]
        blbl = f'Merge  — blue (diam {da}) + red (diam {db})'
        albl = f'Merge done  — orange (diam {dout})'
        linfo = f'diam {da}+{db} → {dout}'
        print(f"  Step {step+1}: merge  {da}+{db} → {dout}  "
              f"K {len(prev)}→{len(current)}  SI={n_si} I={n_i}")

    frames.append((step,     prev,                     bhl, blbl, n_si, n_i, linfo))
    frames.append((step+1,   [l[:] for l in current],  ahl, albl, n_si, n_i, linfo))

print(f"Rendering {len(frames)} frames...")
plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'cm'})
fig, ax = plt.subplots(figsize=(7, 7.5), facecolor='white')
fig.subplots_adjust(left=0.03, right=0.97, top=0.93, bottom=0.08)
info_text = fig.text(0.50, 0.02, '', ha='center', va='bottom', fontsize=10, color='#333333')

def update(frame_idx):
    ax.clear()
    step, state, highlight, label, n_si, n_i, linfo = frames[frame_idx]
    title = f'Step {step}' + (f'\n{label}' if label else '')
    draw_soup(ax, state, N, title=title, highlight=highlight)
    info_text.set_text(
        f'$\\varepsilon$ = {EPSILON}    N = {N}    K = {len(state)}    '
        f'SI = {n_si}    I = {n_i}'
        + (f'    {linfo}' if linfo else '')
    )

anim = animation.FuncAnimation(fig, update, frames=len(frames), interval=1000//FPS)

outpath = Path(__file__).resolve().parent / "loopsoup-single-loop.gif"
print("Saving GIF...")
anim.save(outpath, writer='pillow', fps=FPS, dpi=120)
plt.close()
print(f"Saved → {outpath}")
