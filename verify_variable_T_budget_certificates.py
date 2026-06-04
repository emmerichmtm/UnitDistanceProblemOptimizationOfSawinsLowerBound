#!/usr/bin/env python3
"""Verify suitable-budget variable-T addendum certificates.

This script is dependency-free and is designed to live in the flat repository
root together with verify_certificates.py. It generates the deterministic
variable-T certificates used in the addendum, verifies them with the same
finite arithmetic pipeline as the main report, and writes JSON/TXT/PGFPlots
outputs.

Important:
- Only certificates that pass verify_certificate(...) are written to the
  verified results table/plot.
- Infeasible attempted budgets, such as the first-61-odd-prime T, are written
  only to the rejected output.
"""

from __future__ import annotations

import heapq
import json
import math
from decimal import Decimal
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from verify_certificates import Certificate, verify_certificate, delta_components, format_text_report

PRIME_LIMIT = 300_000
SUITABLE_BUDGETS = [17, 27, 39, 57, 58, 59, 62, 66, 67, 81, 91, 101, 121]
REJECTED_BUDGETS_TO_CHECK = [61]


def primes_upto(n: int) -> List[int]:
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, int(n**0.5) + 1):
        if sieve[p]:
            start = p * p
            sieve[start : n + 1 : p] = b"\x00" * (((n - start) // p) + 1)
    return [i for i, value in enumerate(sieve) if value]


PRIMES = primes_upto(PRIME_LIMIT)


@lru_cache(None)
def legendre_symbol(a: int, p: int) -> int:
    a %= p
    if a == 0:
        return 0
    residue = pow(a, (p - 1) // 2, p)
    return -1 if residue == p - 1 else residue


def first_odd_primes(count: int) -> List[int]:
    return [p for p in PRIMES if p > 2][:count]


def status_in_q_sqrt_product_t(p: int, T: List[int]) -> str:
    if p == 2:
        return "ramified"
    d_mod = 1
    for q in T:
        d_mod = (d_mod * (q % p)) % p
        if d_mod == 0:
            return "ramified"
    symbol = legendre_symbol(d_mod, p)
    return "split" if symbol == 1 else "inert"


def admissible(p: int, T: List[int], T_set: set[int]) -> bool:
    if p == 2 or p in T_set or p % 4 == 1:
        return True
    return any(q != p and legendre_symbol(q, p) == -1 for q in T)


def parity_feasible_T(T: List[int]) -> bool:
    return sum(1 for q in T if q % 4 == 3) % 2 == 1


def available_s_q_budget(T: List[int]) -> int:
    return ((len(T) - 1) ** 2) // 4 - len(T) - 1


def select_no_split_admissible_S_Q(T: List[int], budget: int) -> List[int]:
    T_set = set(T)
    selected: List[int] = []
    for p in PRIMES:
        if admissible(p, T, T_set) and status_in_q_sqrt_product_t(p, T) != "split":
            selected.append(p)
            if len(selected) >= budget:
                return selected
    raise RuntimeError(f"prime table too small: selected {len(selected)} of {budget}")


def fast_delta_numerator_denominator(T: List[int], S_Q: List[int], k: List[int], R: float) -> tuple[float, float]:
    T_set = set(T)
    ln_prod_T = sum(math.log(q) for q in T)
    numerator = (
        math.log(1.0 - 1.0 / R)
        + 0.5 * math.log(2.0 * math.pi / math.e)
        - 0.125 * (math.log(4.0) + ln_prod_T)
        - 0.5 * math.log(0.5 * (math.log(4.0) + ln_prod_T))
    )
    log_product_term = 0.0
    for p, k_p in zip(S_Q, k):
        e = 2 if p == 2 or p in T_set else 1
        numerator += (1.0 / (4.0 * e)) * math.log(k_p + 1.0)
        log_product_term += (k_p / (2.0 * e)) * math.log(p)
    denominator = math.log(2.0 * R) + log_product_term
    return numerator, denominator


def fast_delta(T: List[int], S_Q: List[int], k: List[int], R: float) -> float:
    numerator, denominator = fast_delta_numerator_denominator(T, S_Q, k, R)
    return numerator / denominator


def optimize_k_for_fixed_R(T: List[int], S_Q: List[int], R: float, max_steps: int = 200_000) -> List[int]:
    """Greedy integer multiplicity optimization by marginal improvement.

    It starts from k(p)=1 for every selected prime and repeatedly applies the
    best integer increment whose marginal ratio still improves the current
    quotient numerator/denominator.
    """
    T_set = set(T)
    k = [1] * len(S_Q)
    numerator, denominator = fast_delta_numerator_denominator(T, S_Q, k, R)
    heap: list[tuple[float, int, float, float]] = []
    for i, p in enumerate(S_Q):
        e = 2 if p == 2 or p in T_set else 1
        a = 1.0 / (4.0 * e)
        b = math.log(p) / (2.0 * e)
        next_ratio = a * (math.log(3.0) - math.log(2.0)) / b
        heapq.heappush(heap, (-next_ratio, i, a, b))
    steps = 0
    while heap and steps < max_steps:
        ratio, i, a, b = heapq.heappop(heap)
        ratio = -ratio
        current = numerator / denominator
        if ratio <= current + 1e-15:
            break
        old_k = k[i]
        numerator += a * (math.log(old_k + 2.0) - math.log(old_k + 1.0))
        denominator += b
        k[i] += 1
        new_k = k[i]
        new_ratio = a * (math.log(new_k + 2.0) - math.log(new_k + 1.0)) / b
        heapq.heappush(heap, (-new_ratio, i, a, b))
        steps += 1
    return k


def optimize_certificate_for_budget(N: int) -> Optional[Certificate]:
    T = first_odd_primes(N)
    if not parity_feasible_T(T):
        return None
    budget = available_s_q_budget(T)
    if budget <= 0:
        return None
    S_Q = select_no_split_admissible_S_Q(T, budget)

    # Coarse then local fine search over R.
    best: Optional[tuple[float, float, List[int]]] = None
    for R in [r / 2.0 for r in range(20, 121)]:
        k = optimize_k_for_fixed_R(T, S_Q, R)
        d = fast_delta(T, S_Q, k, R)
        if best is None or d > best[0]:
            best = (d, R, k)
    assert best is not None
    _, R0, _ = best
    for R in [round(R0 - 1.0 + 0.02 * i, 2) for i in range(101) if R0 - 1.0 + 0.02 * i > 1.01]:
        k = optimize_k_for_fixed_R(T, S_Q, R)
        d = fast_delta(T, S_Q, k, R)
        if d > best[0]:
            best = (d, R, k)

    _, R, k_values = best
    return Certificate(
        key=f"variable_T_budget_{N}",
        name=f"Verified variable-T budget certificate N={N}",
        T=T,
        S_Q=S_Q,
        k={p: int(k) for p, k in zip(S_Q, k_values)},
        R=Decimal(str(R)),
    )


def attempted_infeasible_certificate(N: int) -> Certificate:
    """Build the same deterministic candidate even if T fails parity, for rejection auditing."""
    T = first_odd_primes(N)
    budget = available_s_q_budget(T)
    S_Q = select_no_split_admissible_S_Q(T, budget)
    best: Optional[tuple[float, float, List[int]]] = None
    for R in [round(32.8 + 0.01 * i, 2) for i in range(81)]:
        k = optimize_k_for_fixed_R(T, S_Q, R)
        d = fast_delta(T, S_Q, k, R)
        if best is None or d > best[0]:
            best = (d, R, k)
    assert best is not None
    _, R, k_values = best
    return Certificate(
        key=f"rejected_variable_T_budget_{N}",
        name=f"Rejected variable-T budget candidate N={N}",
        T=T,
        S_Q=S_Q,
        k={p: int(k) for p, k in zip(S_Q, k_values)},
        R=Decimal(str(R)),
    )


def plot_coordinates(results: List[Dict[str, object]]) -> str:
    return "\n  ".join(f"({len(r['T'])},{r['delta']})" for r in results)


def write_plot(results: List[Dict[str, object]], path: Path) -> None:
    coords = "\n  ".join(f"({len(r['T'])},{Decimal(str(r['delta'])):.18f})" for r in results)
    text = rf"""% Plot input generated by verify_variable_T_budget_certificates.py.
\begin{{tikzpicture}}
\begin{{axis}}[
  width=0.93\linewidth,
  height=0.56\linewidth,
  xlabel={{ramified-prime budget $N=\#T$}},
  ylabel={{verified exponent gain $\delta$}},
  xmin=10, xmax=126,
  ymin=0.012, ymax=0.0335,
  xtick={{13,17,27,39,57,62,67,81,101,121}},
  ytick={{0.015,0.020,0.025,0.030}},
  scaled y ticks=false,
  yticklabel style={{/pgf/number format/fixed,/pgf/number format/precision=3}},
  grid=both,
  legend columns=3,
  legend style={{font=\scriptsize, at={{(0.5,-0.18)}}, anchor=north, fill=white, fill opacity=0.92, draw=black!20}},
  clip=false,
]
\addplot+[mark=*, mark size=1.8pt, thick] coordinates {{
  {coords}
}};
\addlegendentry{{verified variable-$T$ reruns}}

\addplot+[only marks, mark=square*, mark size=3.2pt] coordinates {{(13,0.0141144286784982)}};
\addlegendentry{{Sawin published reference}}

\addplot+[only marks, mark=triangle*, mark size=3.3pt] coordinates {{(13,0.0152628688170072)}};
\addlegendentry{{Emmerich v1 best fixed-$T$}}

\draw[dashed, black!40] (axis cs:67,0.012) -- (axis cs:67,0.03118533444317659);
\node[font=\scriptsize, anchor=west, align=left] at (axis cs:69,0.0320)
  {{$N=67$\\$\delta=0.0311853344\ldots$}};
\node[font=\scriptsize, anchor=west, align=left] at (axis cs:82,0.0286)
  {{declining side:\\root-discriminant cost dominates}};
\node[font=\scriptsize, anchor=west, align=left] at (axis cs:22,0.0190)
  {{rising side:\\larger budget for $S_Q$}};
\end{{axis}}
\end{{tikzpicture}}
"""
    path.write_text(text, encoding="utf-8")


def main() -> None:
    verified_results: List[Dict[str, object]] = []
    rejected_results: List[Dict[str, object]] = []

    for N in SUITABLE_BUDGETS:
        cert = optimize_certificate_for_budget(N)
        if cert is None:
            continue
        result = verify_certificate(cert)
        if not result["passed"]:
            raise RuntimeError(f"Budget N={N} was expected to pass but failed: {result['checks']}")
        verified_results.append(result)

    for N in REJECTED_BUDGETS_TO_CHECK:
        cert = attempted_infeasible_certificate(N)
        result = verify_certificate(cert)
        rejected_results.append(result)

    Path("verified_variable_T_budget_certificates.json").write_text(json.dumps(verified_results, indent=2, sort_keys=True), encoding="utf-8")
    Path("verified_variable_T_budget_certificates.txt").write_text(format_text_report(verified_results) + "\n", encoding="utf-8")
    Path("rejected_variable_T_budget_candidates.json").write_text(json.dumps(rejected_results, indent=2, sort_keys=True), encoding="utf-8")
    Path("rejected_variable_T_budget_candidates.txt").write_text(format_text_report(rejected_results) + "\n", encoding="utf-8")
    write_plot(verified_results, Path("delta_curve_T_budget_verified_pgfplots.tex"))

    best = max(verified_results, key=lambda r: Decimal(str(r["delta"])))
    print(f"Verified {len(verified_results)} variable-T certificates.")
    print(f"Best verified: N={len(best['T'])}, R={best['R']}, delta={best['delta']}")
    for rejected in rejected_results:
        failed = [name for name, ok in rejected["checks"].items() if not ok]
        print(f"Rejected N={len(rejected['T'])}: failed {failed}")


if __name__ == "__main__":
    main()
