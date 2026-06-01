#!/usr/bin/env python3
"""Construct a classical explicit planar comparator.

This is not the Sawin/OpenAI algebraic construction. It is the elementary
multi-direction lattice comparator used in the report to show a genuine finite
coordinate point set with more unit-distance edges than the ordinary square grid
of the same dimensions.

No third-party Python packages are required.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Tuple

Point = Tuple[float, float]


def build_points(m: int, d: int) -> List[Point]:
    scale = math.sqrt(d)
    return [(i / scale, j / scale) for i in range(-m, m + 1) for j in range(-m, m + 1)]


def count_unit_edges(m: int, d: int) -> Dict[str, int]:
    # For d=5, integer differences with squared length 5 are (±1,±2),(±2,±1).
    vectors = [(1, 2), (1, -2), (-1, 2), (-1, -2), (2, 1), (2, -1), (-2, 1), (-2, -1)]
    vertices = [(i, j) for i in range(-m, m + 1) for j in range(-m, m + 1)]
    vertex_set = set(vertices)
    edges = set()
    for i, j in vertices:
        for di, dj in vectors:
            other = (i + di, j + dj)
            if other in vertex_set:
                edges.add(tuple(sorted(((i, j), other))))
    square_grid_edges = 2 * (2 * m + 1) * (2 * m)
    return {
        "m": m,
        "d": d,
        "n_points": len(vertices),
        "unit_edges_multi_direction": len(edges),
        "unit_edges_square_grid_same_side_length": square_grid_edges,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m", type=int, default=10, help="integer grid half-width")
    parser.add_argument("--d", type=int, default=5, help="squared length used for scaling; currently d=5")
    parser.add_argument("--json-out", default="explicit_planar_comparator.json")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.d != 5:
        raise SystemExit("This simple comparator currently implements d=5 only.")
    summary = count_unit_edges(args.m, args.d)
    print(json.dumps(summary, indent=2, sort_keys=True))
    Path(args.json_out).write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")


if __name__ == "__main__":
    main()
