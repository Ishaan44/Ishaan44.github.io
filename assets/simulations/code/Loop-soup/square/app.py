import matplotlib
matplotlib.use('Agg')

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

# ── Physics / math functions ──────────────────────────────────────────────────

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
                label = True
                rmax = 0
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

def _build_st_max(arr):
    n = len(arr); st = [arr]; j = 1
    while (1 << j) <= n:
        half = 1 << (j-1); end = n - (1 << j) + 1
        st.append(np.maximum(st[j-1][:end], st[j-1][half:half+end]))
        j += 1
    return st

def _build_st_min(arr):
    n = len(arr); st = [arr]; j = 1
    while (1 << j) <= n:
        half = 1 << (j-1); end = n - (1 << j) + 1
        st.append(np.minimum(st[j-1][:end], st[j-1][half:half+end]))
        j += 1
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
        return soup, None, n_splits, n_merges
    r = np.random.randint(n_total)
    if r < n_splits:
        i, s, t = split_ops[r]
        l1, l2 = split_loop(soup[i], s, t)
        return [l for k, l in enumerate(soup) if k != i] + [l1, l2], ('split', i, s, t), n_splits, n_merges
    else:
        r2 = r - n_splits; cumulative = 0
        for (ia, ib), count in merge_count.items():
            if cumulative + count > r2:
                pairs = I_pairs(soup[ia], soup[ib])
                s, t = pairs[r2 - cumulative]
                new_loop = merge_loops(soup[ia], soup[ib], s, t)
                return [l for k, l in enumerate(soup) if k != ia and k != ib] + [new_loop], ('merge', ia, ib, s, t), n_splits, n_merges
            cumulative += count
    return soup, None, n_splits, n_merges

# ── Drawing ───────────────────────────────────────────────────────────────────
# highlight: list of (loop_idx, color, linestyle)  — linestyle '-' or '--'

def draw_soup(ax, soup, N, title='', highlight=()):
    ax.set_facecolor('white')
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_aspect('equal'); ax.axis('off')

    # Square boundary
    ax.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], '-', color='#2C3E50', linewidth=1.5, zorder=1)

    if not soup:
        ax.set_title(title, fontsize=10)
        return
    hl_indices = {h[0] for h in highlight}
    ax.set_title(title, fontsize=10)
    nan = [float('nan')]
    xs_bg, ys_bg = [], []
    for idx, loop in enumerate(soup):
        if idx in hl_indices:
            continue
        xs_bg += [p[1] / (N-1) for p in loop] + nan
        ys_bg += [p[0] / (N-1) for p in loop] + nan
    if xs_bg:
        ax.plot(xs_bg, ys_bg, '-', color='#1A1A1A', linewidth=0.45, alpha=0.88, zorder=2)
    for h in highlight:
        idx, color = h[0], h[1]
        ls = h[2] if len(h) > 2 else '-'
        loop = soup[idx]
        xs = [p[1] / (N-1) for p in loop]
        ys = [p[0] / (N-1) for p in loop]
        lw = 2.5 if ls == '-' else 2.0
        ax.plot(xs, ys, ls, color=color, linewidth=lw, alpha=0.95, zorder=4)

def fig_to_rgb(fig):
    fig.canvas.draw()
    return np.asarray(fig.canvas.buffer_rgba())[:, :, :3]

def render_frame(soup, N, highlight, title, figsize=(7, 7), dpi=150):
    plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'cm'})
    fig, ax = plt.subplots(figsize=figsize, facecolor='white', dpi=dpi)
    draw_soup(ax, soup, N, title=title, highlight=highlight)
    img = fig_to_rgb(fig)
    plt.close(fig)
    return img

def render_grid(panels, N, ncols=2, figsize_per=(5, 5), dpi=120):
    """panels: list of (soup, highlight, title)"""
    plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'cm'})
    k = len(panels)
    ncols = min(ncols, k)
    nrows = (k + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(figsize_per[0]*ncols, figsize_per[1]*nrows),
                              facecolor='white', dpi=dpi)
    axs = np.array(axes).flatten() if k > 1 else [axes]
    for i, (soup, highlight, title) in enumerate(panels):
        draw_soup(axs[i], soup, N, title=title, highlight=highlight)
    for i in range(k, len(axs)):
        axs[i].set_visible(False)
    fig.tight_layout(pad=1.5)
    img = fig_to_rgb(fig)
    plt.close(fig)
    return img

# ── Chain runner ──────────────────────────────────────────────────────────────

def run_chain(soup, eps, N, n_steps):
    """Returns list of (state, highlight, step, n_si, n_i).
       highlight entries are (idx, color, linestyle).
    """
    CA = '#2471A3'   # blue  — solid
    CB = '#C0392B'   # red   — dashed (so both visible when overlapping)
    CR = '#D4AC0D'   # gold
    frames = [([lp[:] for lp in soup], [], 0, 0, 0)]
    current = [lp[:] for lp in soup]
    for step in range(n_steps):
        prev = [lp[:] for lp in current]
        current, op, n_si, n_i = rewiring_step(current, eps, N)
        if op is None:
            break
        if op[0] == 'split':
            _, i, _, _ = op
            frames.append((prev,                       [(i, CR, '-')],                                           step,   n_si, n_i))
            frames.append(([lp[:] for lp in current], [(len(current)-2, CA, '-'), (len(current)-1, CB, '--')],  step+1, n_si, n_i))
        else:
            _, ia, ib, _, _ = op
            frames.append((prev,                       [(ia, CA, '-'), (ib, CB, '--')],  step,   n_si, n_i))
            frames.append(([lp[:] for lp in current], [(len(current)-1, CR, '-')],       step+1, n_si, n_i))
    return frames

# ── Navigation widget ─────────────────────────────────────────────────────────

def nav_buttons(key, n_frames):
    """Renders ← / step info / → and returns current index."""
    fidx = st.session_state.get(key, 0)
    c1, c2, c3 = st.columns([1, 3, 1])
    with c1:
        if st.button("← Prev", key=f'{key}_prev', disabled=(fidx == 0)):
            st.session_state[key] = fidx - 1
            st.rerun()
    with c2:
        st.markdown(
            f"<div style='text-align:center; padding-top:6px'>"
            f"frame <b>{fidx+1}</b> / {n_frames}</div>",
            unsafe_allow_html=True
        )
    with c3:
        if st.button("Next →", key=f'{key}_next', disabled=(fidx == n_frames - 1)):
            st.session_state[key] = fidx + 1
            st.rerun()
    return st.session_state.get(key, 0)

# ── Streamlit app ─────────────────────────────────────────────────────────────

st.set_page_config(page_title="Loop Soup (Square)", layout="wide", page_icon="□")
plt.rcParams.update({'font.family': 'serif', 'mathtext.fontset': 'cm'})

st.markdown(r"# Loop Soup on the Square  $M^\varepsilon$")

# ── Shared sidebar ────────────────────────────────────────────────────────────

with st.sidebar:
    st.header("Soup parameters")
    N = st.slider("Grid size N", 50, 500, 200, 50)
    c = st.slider("Intensity c", 0.5, 4.0, 1.0, 0.5)
    gen_btn = st.button("Generate soup", type="primary", use_container_width=True)

if gen_btn:
    with st.spinner(f"Generating loop soup  (N={N}, c={c})…"):
        soup = LoopSoup(N, c, epsilon=0.0)
    st.session_state['soup'] = soup
    st.session_state['N']    = N
    st.session_state['c']    = c
    st.session_state.pop('chain_imgs', None)
    st.session_state.pop('cmp_imgs',   None)
    st.sidebar.success(f"{len(soup)} loops generated.")

# ── Tabs ──────────────────────────────────────────────────────────────────────

tab_soup, tab_chain, tab_cmp = st.tabs(["Loop Soup", "Rewiring Chain", "Compare ε"])

# ── Tab 1: Loop Soup ──────────────────────────────────────────────────────────

with tab_soup:
    if 'soup' not in st.session_state:
        st.info("Click **Generate soup** in the sidebar to begin.")
    else:
        soup = st.session_state['soup']
        N_s  = st.session_state['N']
        c_s  = st.session_state['c']
        st.write(f"**{len(soup)} loops**   ·   N = {N_s}   ·   c = {c_s}")
        title = f'Loop Soup,   $c = {c_s}$,   $N = {N_s}$\n$K = {len(soup)}$ loops'
        img = render_frame(soup, N_s, [], title, figsize=(7, 7), dpi=150)
        st.image(img, use_container_width=False)

# ── Tab 2: Rewiring Chain ─────────────────────────────────────────────────────

with tab_chain:
    if 'soup' not in st.session_state:
        st.info("Click **Generate soup** in the sidebar first.")
    else:
        col_ctrl, col_viz = st.columns([1, 3])
        with col_ctrl:
            eps_c   = st.slider("ε (cutoff)", 0.001, 0.20, 0.05, 0.001,
                                 format="%.3f", key='eps_chain')
            steps_c = st.slider("Rewiring steps", 5, 200, 20, 5, key='steps_chain')
            run_btn = st.button("▶ Run chain", type="primary", use_container_width=True)

        if run_btn:
            soup = st.session_state['soup']
            N_s  = st.session_state['N']
            c_s  = st.session_state['c']
            with st.spinner("Running rewiring chain…"):
                raw = run_chain(soup, eps_c, N_s, steps_c)
            prog = st.progress(0, text="Rendering frames…")
            imgs = []
            for k, (state, highlight, step, n_si, n_i) in enumerate(raw):
                title = (f'$M^{{\\varepsilon}}$,  $\\varepsilon={eps_c}$,'
                         f'  $c={c_s}$,  $N={N_s}$\n'
                         f'step {step}  ·  $K={len(state)}$'
                         f'  ·  SI={n_si}  ·  I={n_i}')
                imgs.append(render_frame(state, N_s, highlight, title))
                prog.progress((k+1)/len(raw), text=f"Rendering {k+1}/{len(raw)}…")
            prog.empty()
            st.session_state['chain_imgs'] = imgs
            st.session_state['chain_fidx'] = 0

        with col_viz:
            if 'chain_imgs' in st.session_state:
                imgs = st.session_state['chain_imgs']
                fidx = nav_buttons('chain_fidx', len(imgs))
                st.image(imgs[fidx], use_container_width=False)
            elif not run_btn:
                st.info("Set ε and steps, then click **▶ Run chain**.")

# ── Tab 3: Compare ε ─────────────────────────────────────────────────────────

with tab_cmp:
    if 'soup' not in st.session_state:
        st.info("Click **Generate soup** in the sidebar first.")
    else:
        col_ctrl, col_viz = st.columns([1, 3])
        with col_ctrl:
            eps_options = [0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]
            eps_sel = st.multiselect(
                "ε values to compare",
                options=eps_options,
                default=[0.05, 0.01, 0.005, 0.001],
                format_func=str,
                key='eps_cmp'
            )
            steps_cmp = st.slider("Rewiring steps", 5, 200, 20, 5, key='steps_cmp')
            cmp_btn   = st.button("▶ Compare", type="primary", use_container_width=True)

        if cmp_btn and eps_sel:
            soup = st.session_state['soup']
            N_s  = st.session_state['N']
            c_s  = st.session_state['c']
            chain_frames = {}
            for eps_v in eps_sel:
                with st.spinner(f"Running chain  ε = {eps_v}…"):
                    chain_frames[eps_v] = run_chain(soup, eps_v, N_s, steps_cmp)
            max_len = max(len(f) for f in chain_frames.values())
            prog = st.progress(0, text="Rendering comparison frames…")
            imgs = []
            for t in range(max_len):
                panels = []
                for eps_v in eps_sel:
                    fv = chain_frames[eps_v]
                    state, highlight, step, n_si, n_i = fv[min(t, len(fv)-1)]
                    title = (f'$\\varepsilon = {eps_v}$  ·  step {step}'
                             f'  ·  $K = {len(state)}$\n'
                             f'SI = {n_si}  ·  I = {n_i}')
                    panels.append((state, highlight, title))
                ncols = 2 if len(eps_sel) > 1 else 1
                imgs.append(render_grid(panels, N_s, ncols=ncols))
                prog.progress((t+1)/max_len, text=f"Rendering frame {t+1}/{max_len}…")
            prog.empty()
            st.session_state['cmp_imgs'] = imgs
            st.session_state['cmp_fidx'] = 0

        with col_viz:
            if 'cmp_imgs' in st.session_state:
                imgs = st.session_state['cmp_imgs']
                fidx = nav_buttons('cmp_fidx', len(imgs))
                st.image(imgs[fidx], use_container_width=False)
            elif not cmp_btn:
                st.info("Select ε values and click **▶ Compare**.")
