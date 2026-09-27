from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import scipy
import pickle
from scipy.special import hyp2f1
from scipy.ndimage import map_coordinates

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
                        loop = path[now:next + 1]
                        loopsoup.append(loop)
                    now = next
                path = path[last + 1:]
    return loopsoup

def loop_diameter(loop, N):
    rows = [p[0] for p in loop]
    cols = [p[1] for p in loop]
    return max(max(rows) - min(rows), max(cols) - min(cols)) / N

N = 1000
c = 0.5
epsilons = [0.15, 0.08, 0.03, 0.005]  # large -> small: see more loops as epsilon decreases

print("Running LoopSoup...")
loopsoup = LoopSoup(N, c)
print(f"Done. {len(loopsoup)} loops found.")

# Precompute diameters once so we don't recompute per panel
diameters = [loop_diameter(loop, N) for loop in loopsoup]

print("Drawing...")
fig, axes = plt.subplots(2, 2, figsize=(12, 12), facecolor='white')

nan = [float('nan')]
for ax, eps in zip(axes.flat, epsilons):
    loopsoup_cutoff = [loop for loop, d in zip(loopsoup, diameters) if d >= eps]
    print(f"  epsilon={eps}: {len(loopsoup_cutoff)} loops")

    ax.set_facecolor('white')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(f'ε = {eps}', fontsize=14)

    boundary = plt.Circle((0.5, 0.5), 0.5, color='black', fill=False, linewidth=1)
    ax.add_patch(boundary)

    xs_all, ys_all = [], []
    for loop in loopsoup_cutoff:
        xs_all += [p[1] / (N - 1) for p in loop] + nan
        ys_all += [p[0] / (N - 1) for p in loop] + nan
    if xs_all:
        line, = ax.plot(xs_all, ys_all, 'k-', linewidth=0.4, alpha=0.8)

    clip_circle = plt.Circle((0.5, 0.5), 0.5, transform=ax.transData)
    for line in ax.lines:
        line.set_clip_path(clip_circle)

plt.tight_layout()
outpath = Path(__file__).resolve().parent / "loopsoup-epsilon-compare.png"
plt.savefig(outpath, bbox_inches='tight', pad_inches=0.1, dpi=300)
print(f"Saved to {outpath}")
