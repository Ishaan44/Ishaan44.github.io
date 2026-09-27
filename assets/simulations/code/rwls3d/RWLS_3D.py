import argparse
import pickle
import time
from pathlib import Path

import numpy as np


def next_tip_3d(N, tip, rng):
    x, y, z = tip
    neighbors = []

    if x > 0:
        neighbors.append((x - 1, y, z))
    if x < N - 1:
        neighbors.append((x + 1, y, z))
    if y > 0:
        neighbors.append((x, y - 1, z))
    if y < N - 1:
        neighbors.append((x, y + 1, z))
    if z > 0:
        neighbors.append((x, y, z - 1))
    if z < N - 1:
        neighbors.append((x, y, z + 1))

    return neighbors[rng.integers(len(neighbors))]


def wired_boundary_cube(N):
    tree = np.zeros((N, N, N), dtype=bool)

    tree[0, :, :] = True
    tree[N - 1, :, :] = True
    tree[:, 0, :] = True
    tree[:, N - 1, :] = True
    tree[:, :, 0] = True
    tree[:, :, N - 1] = True

    return tree


def loop_soup_3d(N, c, seed=None, progress_every=10):
    """
    Random-walk loop soup in a wired N x N x N cube.

    This is a direct 3D analogue of the existing 2D RWLS.py:
    - the cube boundary is wired into the initial tree;
    - each unwired interior vertex launches a random walk until it hits the tree;
    - loops are extracted from repeated vertices in the path;
    - each loop is kept according to the RWLS label rule controlled by c.
    """
    rng = np.random.default_rng(seed)
    loopsoup = []
    tree = wired_boundary_cube(N)
    start_time = time.time()

    processed_path_vertices = 0

    for x in range(1, N - 1):
        if progress_every and (x == 1 or x % progress_every == 0):
            elapsed = time.time() - start_time
            print(
                f"x-slice {x}/{N - 2} | "
                f"loops={len(loopsoup)} | "
                f"processed path vertices={processed_path_vertices} | "
                f"elapsed={elapsed:.1f}s",
                flush=True,
            )

        for y in range(1, N - 1):
            for z in range(1, N - 1):
                if tree[x, y, z]:
                    continue

                tip = (x, y, z)
                path = [tip]

                while not tree[tip]:
                    tip = next_tip_3d(N, tip, rng)
                    path.append(tip)

                while path:
                    root = path[0]
                    tree[root] = True
                    processed_path_vertices += 1

                    label = True
                    rmax = 0.0
                    last = len(path) - 1 - path[::-1].index(root)
                    now = 0

                    while now < last:
                        next_index = now + 1 + path[now + 1 :].index(root)

                        r = rng.random()
                        if r > rmax:
                            rmax = r
                            label = rng.random() < c

                        if label:
                            loopsoup.append(path[now : next_index + 1])

                        now = next_index

                    path = path[last + 1 :]

    return loopsoup


def loop_stats(loopsoup):
    lengths = np.array([len(loop) - 1 for loop in loopsoup], dtype=int)

    if len(lengths) == 0:
        return {
            "count": 0,
            "min_length": 0,
            "max_length": 0,
            "mean_length": 0.0,
            "total_steps": 0,
        }

    return {
        "count": int(len(lengths)),
        "min_length": int(lengths.min()),
        "max_length": int(lengths.max()),
        "mean_length": float(lengths.mean()),
        "total_steps": int(lengths.sum()),
    }


def save_loops(output_path, loopsoup, N, c, seed):
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    payload = {
        "N": N,
        "c": c,
        "seed": seed,
        "dimension": 3,
        "loopsoup": loopsoup,
        "stats": loop_stats(loopsoup),
    }

    with output_path.open("wb") as f:
        pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)


def parse_args():
    parser = argparse.ArgumentParser(description="Simulate a 3D random-walk loop soup in a cube.")
    parser.add_argument("--N", type=int, default=30, help="Cube side length.")
    parser.add_argument("--c", type=float, default=0.5, help="Loop soup intensity/retention parameter.")
    parser.add_argument("--seed", type=int, default=None, help="Random seed.")
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output pickle path. Defaults to rwls3d/output/loopsoup-3d-N{N}-c{c}.pkl.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=5,
        help="Print progress every this many x-slices. Use 0 to disable.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    if args.N < 3:
        raise ValueError("N must be at least 3 so the cube has an interior.")
    if not 0 <= args.c <= 1:
        raise ValueError("c must be between 0 and 1.")

    output = args.output
    if output is None:
        c_label = str(args.c).replace(".", "p")
        output = Path(__file__).resolve().parent / "output" / f"loopsoup-3d-N{args.N}-c{c_label}.pkl"

    print(f"Running 3D RWLS in a {args.N} x {args.N} x {args.N} cube...")
    print(f"c={args.c}, seed={args.seed}")

    start_time = time.time()
    loopsoup = loop_soup_3d(args.N, args.c, seed=args.seed, progress_every=args.progress_every)
    elapsed = time.time() - start_time

    stats = loop_stats(loopsoup)
    print("Done.")
    print(f"Loops found: {stats['count']}")
    print(f"Loop length: min={stats['min_length']}, mean={stats['mean_length']:.2f}, max={stats['max_length']}")
    print(f"Total loop steps: {stats['total_steps']}")
    print(f"Elapsed: {elapsed:.1f}s")

    save_loops(output, loopsoup, args.N, args.c, args.seed)
    print(f"Saved to {output}")


if __name__ == "__main__":
    main()
