from pathlib import Path

import matplotlib
matplotlib.use('Agg')  # non-interactive backend — avoids macOS/VS Code display issues
import numpy as np
import matplotlib.pyplot as plt

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

def loop_diameter(loop, N):
    rows = [p[0] for p in loop]
    cols = [p[1] for p in loop]
    # bounding box diameter in scaled coordinates (so epsilon lives in [0,1])
    return max(max(rows) - min(rows), max(cols) - min(cols)) / N

N = 1000
c = 0.5
epsilon = 0.001  # keep loops whose spatial diameter is >= epsilon * domain size

print("Running LoopSoup...")
loopsoup = LoopSoup(N, c)
print(f"Done. {len(loopsoup)} loops found.")

# Apply cutoff: keep only loops with spatial diameter >= epsilon
loopsoup_cutoff = [loop for loop in loopsoup if loop_diameter(loop, N) >= epsilon]
print(f"After cutoff (epsilon={epsilon}): {len(loopsoup_cutoff)} loops kept.")

print("Drawing loops...")
fig, ax = plt.subplots(figsize=(8, 8), facecolor='white')
ax.set_facecolor('white')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_aspect('equal')
ax.axis('off')

# Square boundary
ax.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], '-', color='black', linewidth=1)

# Draw all loops in black
nan = [float('nan')]
xs_all, ys_all = [], []
for loop in loopsoup_cutoff:
    xs_all += [p[1] / (N - 1) for p in loop] + nan
    ys_all += [p[0] / (N - 1) for p in loop] + nan
if xs_all:
    ax.plot(xs_all, ys_all, 'k-', linewidth=0.5, alpha=0.8)

outpath = Path(__file__).resolve().parent / f"loopsoup-{N}-diam-{epsilon}-square.png"
plt.savefig(outpath, bbox_inches='tight', pad_inches=0.05, dpi=300)
print(f"Saved to {outpath}")
