from pathlib import Path

import matplotlib
matplotlib.use('Agg')  # non-interactive backend — avoids macOS/VS Code display issues
import numpy as np
import matplotlib.pyplot as plt
import scipy
import pickle

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
                    if (r := np.random.uniform()) > rmax:  # fixed walrus operator precedence
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

N = 10000 # use small N to test; change back to 600 once confirmed working
c = .5
print("Running LoopSoup...")
loopsoup = LoopSoup(N, c)
print(f"Done. {len(loopsoup)} loops found.")
with open(Path(__file__).resolve().parent / f"loopsoup-{N}-.5.pkl", "wb") as f:
    pickle.dump(loopsoup, f)
print("Computing carpet...")
carpet = Carpet(N, loopsoup)
print("Saving image...")
plt.imshow(carpet, cmap='binary', interpolation='nearest')
plt.axis('off')
outpath = Path(__file__).resolve().parent / f"loopsoup-{N}-square.png"
plt.savefig(outpath, bbox_inches='tight', pad_inches=0)
print(f"Saved to {outpath}")
