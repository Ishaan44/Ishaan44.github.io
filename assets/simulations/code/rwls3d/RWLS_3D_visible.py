import argparse
import heapq
import pickle
import time
from pathlib import Path

import numpy as np

from RWLS_3D import next_tip_3d, wired_boundary_cube


def loop_diameter(loop):
    min_x = max_x = loop[0][0]
    min_y = max_y = loop[0][1]
    min_z = max_z = loop[0][2]

    for x, y, z in loop[1:]:
        min_x = min(min_x, x)
        max_x = max(max_x, x)
        min_y = min(min_y, y)
        max_y = max(max_y, y)
        min_z = min(min_z, z)
        max_z = max(max_z, z)

    dx = max_x - min_x
    dy = max_y - min_y
    dz = max_z - min_z
    return float((dx * dx + dy * dy + dz * dz) ** 0.5)


def maybe_store_top_loop(heap, loop, diameter, max_loops, counter):
    length = len(loop) - 1
    score = (diameter, length, counter)
    item = (score, loop)

    if len(heap) < max_loops:
        heapq.heappush(heap, item)
    elif score > heap[0][0]:
        heapq.heapreplace(heap, item)


def visible_loop_soup_3d(
    N,
    c,
    seed=None,
    min_length=20,
    min_diameter=10.0,
    max_loops=10000,
    progress_every=25,
):
    """
    Exact RWLS traversal, but stores only the largest visible loops.

    This is intended for large N, especially N=1000, where saving every
    tiny loop is neither useful nor renderable.
    """
    rng = np.random.default_rng(seed)
    tree = wired_boundary_cube(N)
    start_time = time.time()

    top_loops = []
    processed_path_vertices = 0
    retained_loop_count = 0
    candidate_loop_count = 0
    stored_counter = 0

    for x in range(1, N - 1):
        if progress_every and (x == 1 or x % progress_every == 0):
            elapsed = time.time() - start_time
            print(
                f"x-slice {x}/{N - 2} | "
                f"retained loops={retained_loop_count} | "
                f"visible candidates={candidate_loop_count} | "
                f"stored top loops={len(top_loops)} | "
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
                            retained_loop_count += 1
                            length = next_index - now

                            if length >= min_length:
                                loop = path[now : next_index + 1]
                                diameter = loop_diameter(loop)

                                if diameter >= min_diameter:
                                    candidate_loop_count += 1
                                    maybe_store_top_loop(top_loops, loop, diameter, max_loops, stored_counter)
                                    stored_counter += 1

                        now = next_index

                    path = path[last + 1 :]

    sorted_top = sorted(top_loops, reverse=True)
    loopsoup = [loop for _, loop in sorted_top]

    stats = {
        "retained_loop_count": retained_loop_count,
        "visible_candidate_loop_count": candidate_loop_count,
        "stored_top_loop_count": len(loopsoup),
        "processed_path_vertices": processed_path_vertices,
        "min_length": min_length,
        "min_diameter": min_diameter,
        "max_loops": max_loops,
        "elapsed": time.time() - start_time,
    }

    return loopsoup, stats


def loop_stats(loopsoup):
    if not loopsoup:
        return {"count": 0, "min_length": 0, "max_length": 0, "mean_length": 0.0}

    lengths = np.array([len(loop) - 1 for loop in loopsoup], dtype=int)
    return {
        "count": int(len(lengths)),
        "min_length": int(lengths.min()),
        "max_length": int(lengths.max()),
        "mean_length": float(lengths.mean()),
        "total_steps": int(lengths.sum()),
    }


def save_payload(path, N, c, seed, loopsoup, source_stats):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    payload = {
        "N": N,
        "c": c,
        "seed": seed,
        "dimension": 3,
        "loopsoup": loopsoup,
        "stats": loop_stats(loopsoup),
        "source_stats": source_stats,
        "note": "This file stores only the largest visible loops from an exact RWLS traversal.",
    }

    with path.open("wb") as f:
        pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)


def parse_args():
    parser = argparse.ArgumentParser(description="Generate only the largest visible loops from 3D RWLS.")
    parser.add_argument("--N", type=int, default=1000)
    parser.add_argument("--c", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--min-length", type=int, default=20)
    parser.add_argument("--min-diameter", type=float, default=20.0)
    parser.add_argument("--max-loops", type=int, default=10000)
    parser.add_argument("--progress-every", type=int, default=25)
    parser.add_argument("--output", type=Path, default=None)
    return parser.parse_args()


def main():
    args = parse_args()

    if args.output is None:
        c_label = str(args.c).replace(".", "p")
        args.output = (
            Path(__file__).resolve().parent
            / "output"
            / f"loopsoup-3d-visible-N{args.N}-c{c_label}.pkl"
        )

    print(f"Running visible/top-loop 3D RWLS in a {args.N} x {args.N} x {args.N} cube...")
    print(
        f"c={args.c}, seed={args.seed}, "
        f"min_length={args.min_length}, min_diameter={args.min_diameter}, max_loops={args.max_loops}"
    )

    loopsoup, source_stats = visible_loop_soup_3d(
        args.N,
        args.c,
        seed=args.seed,
        min_length=args.min_length,
        min_diameter=args.min_diameter,
        max_loops=args.max_loops,
        progress_every=args.progress_every,
    )

    stats = loop_stats(loopsoup)
    save_payload(args.output, args.N, args.c, args.seed, loopsoup, source_stats)

    print("Done.")
    print(f"Retained loops seen: {source_stats['retained_loop_count']}")
    print(f"Visible candidates seen: {source_stats['visible_candidate_loop_count']}")
    print(f"Stored top loops: {stats['count']}")
    print(f"Stored loop length: min={stats['min_length']}, mean={stats['mean_length']:.2f}, max={stats['max_length']}")
    print(f"Elapsed: {source_stats['elapsed']:.1f}s")
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    main()
