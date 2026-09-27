from pathlib import Path

import matplotlib
matplotlib.use('Agg')  # non-interactive backend — avoids macOS/VS Code display issues
import numpy as np
import matplotlib.pyplot as plt
import scipy
import math
import heapq
import pickle
from collections import deque
from scipy.special import hyp2f1
from scipy.ndimage import map_coordinates

plt.rcParams['figure.dpi'] = 1200

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
    tree[0, 0] = 1
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
                        loopsoup.append(loop)
                    now = next
                path = path[last + 1:]
    return loopsoup

def Carpet(N, loopsoup):
    config = np.ones((2 * N - 1, 2 * N - 1))
    for loop in loopsoup:
        for i in range(len(loop) - 1):
            config[loop[i][0] + loop[i + 1][0], loop[i][1] + loop[i + 1][1]] = 0
    structure = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=bool)
    carpet, n = scipy.ndimage.label(config, structure=structure)
    for i in range(2 * N - 1):
        for j in range(2 * N - 1):
            if carpet[i, j] > 1:
                carpet[i, j] = 0
    return carpet[1:2 * N - 2:2, 1:2 * N - 2:2]

# Precomputed numerically: max |Re f(e^{i theta})| = max |Im f(e^{i theta})| ~= 0.9270373379
SQUARE_RADIUS = 0.9270373379435793

def disk_to_square_hyp(z, shrink=1 - 1e-8):
    """
    Conformal map from unit disk to square (-1,1)^2:
    f(z) = exp(i*pi/4) * z * 2F1(1/2, 1/4; 5/4; z^4),
    rescaled so the image is exactly the square.
    z: complex ndarray with |z| <= 1
    shrink: pull points slightly inside to avoid the branch point at z^4 = 1
    """
    a, b, c = 0.5, 0.25, 1.25
    z = z * shrink
    w = np.exp(1j * np.pi / 4) * z * hyp2f1(a, b, c, z**4)
    return w / SQUARE_RADIUS

def square_image_to_disk_hyp2f1(img):
    """
    Conformally map a square image (N,N) to a disk (N,N) using hyp2f1.
    - Input: 2D integer array img of shape (N, N)
    - Output: 2D integer array out of shape (N, N)
      with 0 outside the unit circle in the output.
    """
    img = np.asarray(img)
    if img.ndim != 2:
        raise ValueError("img must be a 2D array")
    N, M = img.shape
    if N != M:
        raise ValueError("img must be square")

    y_idx, x_idx = np.indices((N, N))
    c = (N - 1) / 2.0
    x = (x_idx - c) / c
    y = (y_idx - c) / c

    z = x + 1j * y
    r = np.abs(z)
    mask = r <= 1.0

    w_sq = np.zeros_like(z, dtype=np.complex128)
    w_sq[mask] = disk_to_square_hyp(z[mask])

    u = w_sq.real
    v = w_sq.imag

    src_x = (u + 1.0) * 0.5 * (N - 1)
    src_y = (v + 1.0) * 0.5 * (N - 1)

    src_x = np.clip(src_x, 0, N - 1)
    src_y = np.clip(src_y, 0, N - 1)

    coords = np.vstack([src_y[mask].ravel(), src_x[mask].ravel()])
    samples = map_coordinates(img, coords, order=0, mode='nearest')

    out = np.full(img.shape, 0, dtype=float)
    out[mask] = samples

    return out

def loop_diameter(loop, N):
    rows = [p[0] for p in loop]
    cols = [p[1] for p in loop]
    # bounding box diameter in scaled coordinates (so epsilon lives in [0,1])
    return max(max(rows) - min(rows), max(cols) - min(cols)) / N

N = 1000
c = 0.5
epsilon = 0.01  # keep loops whose spatial diameter is >= epsilon * domain size

print("Running LoopSoup...")
loopsoup = LoopSoup(N, c)
print(f"Done. {len(loopsoup)} loops found.")

# Apply cutoff: keep only loops with spatial diameter >= epsilon
loopsoup_cutoff = [loop for loop in loopsoup if loop_diameter(loop, N) >= epsilon]
print(f"After cutoff (epsilon={epsilon}): {len(loopsoup_cutoff)} loops kept.")

with open(Path(__file__).resolve().parent / f"loopsoup-{N}-diam-{epsilon}.pkl", "wb") as f:
    pickle.dump(loopsoup_cutoff, f)

print("Drawing loops...")
fig, ax = plt.subplots(figsize=(8, 8), facecolor='white')
ax.set_facecolor('white')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_aspect('equal')
ax.axis('off')

# Draw a white filled circle as background, black circle as boundary
boundary = plt.Circle((0.5, 0.5), 0.5, color='black', fill=False, linewidth=1)
ax.add_patch(boundary)

# Draw all loops in black
nan = [float('nan')]
xs_all, ys_all = [], []
for loop in loopsoup_cutoff:
    xs_all += [p[1] / (N - 1) for p in loop] + nan
    ys_all += [p[0] / (N - 1) for p in loop] + nan
if xs_all:
    line, = ax.plot(xs_all, ys_all, 'k-', linewidth=0.5, alpha=0.8)

# Clip all lines to the disk
clip_circle = plt.Circle((0.5, 0.5), 0.5, transform=ax.transData)
for line in ax.lines:
    line.set_clip_path(clip_circle)

outpath = Path(__file__).resolve().parent / f"loopsoup-{N}-diam-{epsilon}.png"
plt.savefig(outpath, bbox_inches='tight', pad_inches=0.05, dpi=300)
print(f"Saved to {outpath}")
