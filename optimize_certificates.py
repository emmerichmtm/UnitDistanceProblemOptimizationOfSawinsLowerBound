#!/usr/bin/env python3
"""Greedy parameter search for Sawin-style unit-distance certificates.

The script is intentionally lightweight and reproducible. It first verifies the
built-in Sawin and optimized examples, then performs a deterministic greedy
search using the same default set T as Sawin's published example.

No third-party Python packages are required.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

from verify_certificates import (
    Certificate,
    OPTIMIZED_CERTIFICATE,
    SAWIN_CERTIFICATE,
    admissibility_witness,
    delta_components,
    e_value,
    legendre_symbol,
    product,
    quadratic_status_in_q_sqrt_d,
    verify_certificate,
)


@dataclass(frozen=True)
class SearchResult:
    certificate: Certificate
    delta: Decimal
    numerator: Decimal
    denominator: Decimal
    budget_lhs: int
    budget_rhs: int


def primes_up_to(n: int) -> List[int]:
    """Return all primes <= n."""
    if n < 2:
        return []
    sieve = [True] * (n + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            for j in range(i * i, n + 1, i):
                sieve[j] = False
    return [i for i, flag in enumerate(sieve) if flag]


def pick_k(p: int, c: float) -> int:
    """First-order heuristic for k(p)."""
    return max(1, int(math.floor(1.0 / (2.0 * c * math.log(p)) - 1.0)))


def budget_capacity(T: List[int]) -> int:
    """Available budget for S_Q when no selected prime splits in Q."""
    return ((len(T) - 1) ** 2) // 4 - len(T) - 1


def candidate_prime_cost(p: int, T: List[int]) -> int:
    """Budget cost of a candidate prime p."""
    D = product(T)
    return 2 if quadratic_status_in_q_sqrt_d(p, D) == "split" else 1


def is_admissible_candidate(p: int, T: List[int]) -> bool:
    """Whether p passes the implemented admissibility witness condition."""
    return not admissibility_witness(p, T).startswith("NO")


def greedy_select_S_Q(T: List[int], pmax: int, c: float, allow_split: bool = False) -> Tuple[List[int], Dict[int, int]]:
    """Greedily select S_Q using benefit/cost scores."""
    capacity = budget_capacity(T)
    scored: List[Tuple[float, int, int]] = []
    for p in primes_up_to(pmax):
        if p in T:
            continue
        if not is_admissible_candidate(p, T):
            continue
        cost = candidate_prime_cost(p, T)
        if cost == 2 and not allow_split:
            continue
        k_p = pick_k(p, c)
        benefit = (1.0 / (4.0 * e_value(p, T))) * math.log(k_p + 1.0)
        score = benefit / cost
        scored.append((score, p, cost))
    scored.sort(reverse=True)
    chosen: List[int] = []
    k: Dict[int, int] = {}
    used = 0
    for _, p, cost in scored:
        if used + cost <= capacity:
            chosen.append(p)
            k[p] = pick_k(p, c)
            used += cost
    chosen.sort()
    return chosen, k


def evaluate_certificate(cert: Certificate) -> SearchResult:
    """Evaluate one certificate into a SearchResult."""
    components = delta_components(cert)
    verification = verify_certificate(cert)
    return SearchResult(
        certificate=cert,
        delta=components["delta"],
        numerator=components["numerator"],
        denominator=components["denominator"],
        budget_lhs=int(verification["budget_lhs"]),
        budget_rhs=int(verification["budget_rhs"]),
    )


def search(T: List[int], pmax: int, c: float, r_min: float, r_max: float, r_steps: int) -> SearchResult:
    """Run a deterministic greedy search over S_Q and a grid search over R."""
    S_Q, k = greedy_select_S_Q(T=T, pmax=pmax, c=c)
    best: SearchResult | None = None
    for step in range(r_steps + 1):
        if r_steps == 0:
            R_float = r_min
        else:
            R_float = r_min + (r_max - r_min) * step / r_steps
        cert = Certificate(
            key="greedy_search",
            name="Greedy search candidate",
            T=T,
            S_Q=S_Q,
            k=k,
            R=Decimal(str(R_float)),
        )
        result = evaluate_certificate(cert)
        if best is None or result.delta > best.delta:
            best = result
    assert best is not None
    return best


def result_to_dict(result: SearchResult) -> Dict[str, object]:
    """Serialize SearchResult."""
    cert = result.certificate
    return {
        "name": cert.name,
        "key": cert.key,
        "T": cert.T,
        "S_Q": cert.S_Q,
        "k": {str(p): cert.k[p] for p in cert.S_Q},
        "R": str(cert.R),
        "delta": str(result.delta),
        "numerator": str(result.numerator),
        "denominator": str(result.denominator),
        "budget_lhs": result.budget_lhs,
        "budget_rhs": result.budget_rhs,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pmax", type=int, default=300, help="largest candidate prime considered")
    parser.add_argument("--c", type=float, default=0.015, help="target exponent used in the k(p) heuristic")
    parser.add_argument("--r-min", type=float, default=40.0, help="minimum R for grid search")
    parser.add_argument("--r-max", type=float, default=90.0, help="maximum R for grid search")
    parser.add_argument("--r-steps", type=int, default=500, help="number of R grid intervals")
    parser.add_argument("--json-out", default="optimization_results.json", help="path for JSON output")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    built_in = [evaluate_certificate(SAWIN_CERTIFICATE), evaluate_certificate(OPTIMIZED_CERTIFICATE)]
    greedy = search(
        T=SAWIN_CERTIFICATE.T,
        pmax=args.pmax,
        c=args.c,
        r_min=args.r_min,
        r_max=args.r_max,
        r_steps=args.r_steps,
    )
    payload = {
        "built_in_certificates": [result_to_dict(result) for result in built_in],
        "greedy_search": result_to_dict(greedy),
        "parameters": vars(args),
    }
    print("Built-in certificate deltas:")
    for result in built_in:
        print(f"  {result.certificate.name}: {result.delta}")
    print("Greedy search candidate:")
    print(f"  |S_Q| = {len(greedy.certificate.S_Q)}")
    print(f"  R = {greedy.certificate.R}")
    print(f"  delta = {greedy.delta}")
    Path(args.json_out).write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


if __name__ == "__main__":
    main()
