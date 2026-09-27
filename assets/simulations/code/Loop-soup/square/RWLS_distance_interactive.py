import matplotlib
matplotlib.use('TkAgg')   # change to 'MacOSX' or 'Qt5Agg' if needed
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as mwidgets
import matplotlib.animation as animation
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from collections import defaultdict

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

# ── Pre-render soup backgrounds ───────────────────────────────────────────────

def render_soup_bg(soup, N, size=500):
    fig = Figure(figsize=(size / 100, size / 100), dpi=100, facecolor='white')
    canvas = FigureCanvasAgg(fig)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_facecolor('white')
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_aspect('equal'); ax.axis('off')
    ax.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], '-', color='#2C3E50', linewidth=1.5)
    if soup:
        nan = [float('nan')]
        xs, ys = [], []
        for loop in soup:
            xs += [p[1] / (N - 1) for p in loop] + nan
            ys += [p[0] / (N - 1) for p in loop] + nan
        ax.plot(xs, ys, '-', color='#1A1A1A', linewidth=0.45, alpha=0.88)
    canvas.draw()
    w, h = canvas.get_width_height()
    buf = np.frombuffer(canvas.buffer_rgba(), dtype=np.uint8).copy()
    return buf.reshape(h, w, 4)

# ── Distance overlay (drawn live using current z) ─────────────────────────────

def draw_distance_overlay(ax, soup, N, z):
    D_val, D_idx = compute_D(soup, z, N)
    if D_idx >= 0:
        loop = soup[D_idx]
        xs = [p[1] / (N-1) for p in loop]
        ys = [p[0] / (N-1) for p in loop]
        ax.plot(xs, ys, '-', color='white',   linewidth=9,   alpha=0.7,  zorder=4)
        ax.plot(xs, ys, '-', color='#27AE60', linewidth=4.0, alpha=1.0,  zorder=5)
        rows_s = np.array([p[0] for p in loop], dtype=float) / (N - 1)
        cols_s = np.array([p[1] for p in loop], dtype=float) / (N - 1)
        nearest = int(np.argmin(np.sqrt((rows_s - z[0])**2 + (cols_s - z[1])**2)))
        nx, ny = cols_s[nearest], rows_s[nearest]
        ax.plot([z[1], nx], [z[0], ny], '-', color='#27AE60', linewidth=1.8, alpha=0.9, zorder=6)
        ax.plot(nx, ny, 'o', color='#27AE60', markersize=7, zorder=7,
                markeredgecolor='white', markeredgewidth=1.0)
        theta = np.linspace(0, 2 * np.pi, 300)
        ax.plot(z[1] + D_val * np.cos(theta), z[0] + D_val * np.sin(theta),
                '--', color='#27AE60', linewidth=1.2, alpha=0.5, zorder=3)
    ax.plot(z[1], z[0], 'o', color='black', markersize=6, zorder=11)
    return D_val

# ── Main ────────────────────────────────────────────────────────────────────────

N        = 1000
c        = 0.5
epsilons = [0.15, 0.08, 0.03, 0.01]
n_steps  = 30
fps      = 1   # frames per second for the animation

print("Generating loop soup...")
raw_soup = LoopSoup(N, c)
print(f"Full soup: {len(raw_soup)} loops")
diams = [spatial_diameter(loop, N) for loop in raw_soup]

# Pre-compute chain states and ops (same structure as RWLS_distance_compare)
# frame tuple: (step, soup, op_str)
all_frames = []
for eps in epsilons:
    s = [loop for loop, d in zip(raw_soup, diams) if d >= eps]
    print(f"  ε={eps}: {len(s)} loops — running {n_steps} steps...")
    frames = [(0, s, 'initial state')]
    for step in range(n_steps):
        s, op = rewiring_step(s, eps, N)
        if op is None:
            print(f"    chain stuck at step {step}")
            break
        if op[0] == 'split':
            op_str = f'split loop {op[1]}  →  2 loops  (K: {len(frames[-1][1])} → {len(s)})'
        else:
            op_str = f'merge loops {op[1]} + {op[2]}  →  1 loop  (K: {len(frames[-1][1])} → {len(s)})'
        frames.append((step + 1, s, op_str))
        print(f"    step {step+1}/{n_steps}", end='\r', flush=True)
    print()
    all_frames.append(frames)

# Pre-render every soup state as a pixel array
print("Pre-rendering backgrounds...")
all_images = []
total = sum(len(f) for f in all_frames)
done  = 0
for frames in all_frames:
    imgs = []
    for step, soup, _ in frames:
        imgs.append(render_soup_bg(soup, N))
        done += 1
        print(f"  {done}/{total}", end='\r', flush=True)
    all_images.append(imgs)
print(f"\nReady. Opening window...")

n_total_frames = max(len(f) for f in all_frames)
state = {'z': (0.5, 0.5), 'paused': False, 'frame': 0}

# ── Figure layout ─────────────────────────────────────────────────────────────

plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'cm'})
fig, axes = plt.subplots(2, 2, figsize=(12, 13), facecolor='white')
fig.subplots_adjust(hspace=0.18, wspace=0.06, left=0.03, right=0.97, top=0.95, bottom=0.10)

# Play/Pause button
ax_btn  = fig.add_axes([0.44, 0.02, 0.12, 0.03])
btn     = mwidgets.Button(ax_btn, 'Pause', color='#D6EAF8', hovercolor='#AED6F1')

info_text = fig.text(0.50, 0.065, '', ha='center', va='bottom', fontsize=10, color='#333333')

# ── Animation update ──────────────────────────────────────────────────────────

def draw_frame(frame_idx):
    z = state['z']
    for ax, eps, frames, images in zip(axes.flat, epsilons, all_frames, all_images):
        ax.clear()
        ax.set_facecolor('white')
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_aspect('equal'); ax.axis('off')

        fidx     = min(frame_idx, len(frames) - 1)
        _, soup, op_str = frames[fidx]
        img      = images[fidx]

        ax.imshow(img, extent=[0, 1, 0, 1], origin='upper', aspect='auto', zorder=0)
        D_val = draw_distance_overlay(ax, soup, N, z)
        ax.set_title(
            f'$\\varepsilon$ = {eps}    K = {len(soup)}    D = {D_val:.3f}\n{op_str}',
            fontsize=9
        )

    info_text.set_text(
        f'z = ({z[0]:.3f}, {z[1]:.3f})    Step {frame_idx} / {n_total_frames - 1}'
        f'    Click any panel to move z'
    )

def update(frame_idx):
    state['frame'] = frame_idx
    if state['paused']:
        return
    draw_frame(frame_idx)

def on_pause(_):
    state['paused'] = not state['paused']
    btn.label.set_text('Play' if state['paused'] else 'Pause')

def on_click(event):
    if event.inaxes not in list(axes.flat):
        return
    if event.xdata is None or event.ydata is None:
        return
    state['z'] = (float(np.clip(event.ydata, 0, 1)),
                  float(np.clip(event.xdata, 0, 1)))
    draw_frame(state['frame'])
    fig.canvas.draw_idle()

btn.on_clicked(on_pause)
fig.canvas.mpl_connect('button_press_event', on_click)

anim = animation.FuncAnimation(
    fig, update,
    frames=n_total_frames,
    interval=1000 // fps,
    repeat=True
)

plt.show()
