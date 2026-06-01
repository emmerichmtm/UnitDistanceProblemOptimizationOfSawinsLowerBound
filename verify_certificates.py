#!/usr/bin/env python3
"""Verify explicit unit-distance lower-bound certificates.

This script verifies two finite certificates associated with Sawin's explicit
unit-distance lower-bound criterion:

1. ``sawin``: the published example from the proof of Theorem 1 in
   W. Sawin, "An explicit lower bound for the unit distance problem",
   arXiv:2605.20579.
2. ``optimized``: the optimized candidate reported in the accompanying note.

The script checks finite arithmetic side conditions and evaluates the exponent
formula. It does not construct planar coordinates.

No third-party Python packages are required.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
from pathlib import Path
from typing import Dict, Iterable, List, Optional

getcontext().prec = 90

# A high-precision decimal expansion of pi sufficient for the displayed checks.
PI = Decimal(
    "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899"
)


@dataclass(frozen=True)
class Certificate:
    """Finite data for one Sawin-style certificate."""

    key: str
    name: str
    T: List[int]
    S_Q: List[int]
    k: Dict[int, int]
    R: Decimal
    claimed_delta: Optional[Decimal] = None


SAWIN_CERTIFICATE = Certificate(
    key="sawin",
    name="Sawin published example",
    T=[3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43],
    S_Q=[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 47, 71, 79, 97, 101, 107, 109, 139, 151, 163, 167, 179],
    k={
        2: 50,
        3: 31,
        5: 21,
        7: 17,
        11: 14,
        13: 13,
        17: 12,
        19: 11,
        23: 10,
        29: 10,
        47: 8,
        71: 7,
        79: 7,
        97: 7,
        101: 7,
        107: 7,
        109: 7,
        139: 6,
        151: 6,
        163: 6,
        167: 6,
        179: 6,
    },
    R=Decimal("72"),
    claimed_delta=Decimal("0.014114428678498239"),
)


OPTIMIZED_CERTIFICATE = Certificate(
    key="optimized",
    name="Optimized candidate",
    T=[3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43],
    S_Q=[2, 3, 47, 71, 79, 97, 101, 107, 109, 139, 151, 163, 167, 179, 191, 211, 223, 239, 241, 251, 257, 263],
    k={
        2: 47,
        3: 29,
        47: 7,
        71: 6,
        79: 6,
        97: 6,
        101: 6,
        107: 6,
        109: 6,
        139: 5,
        151: 5,
        163: 5,
        167: 5,
        179: 5,
        191: 5,
        211: 5,
        223: 5,
        239: 5,
        241: 5,
        251: 5,
        257: 5,
        263: 4,
    },
    R=Decimal("66.72240803"),
    claimed_delta=Decimal("0.01517180563721325"),
)


CERTIFICATES = {
    SAWIN_CERTIFICATE.key: SAWIN_CERTIFICATE,
    OPTIMIZED_CERTIFICATE.key: OPTIMIZED_CERTIFICATE,
}


def is_prime(n: int) -> bool:
    """Return True iff n is prime."""
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def product(values: Iterable[int]) -> int:
    """Integer product."""
    result = 1
    for value in values:
        result *= value
    return result


def legendre_symbol(a: int, p: int) -> int:
    """Legendre symbol (a/p), for odd prime p."""
    if p == 2:
        raise ValueError("Legendre symbol is only implemented here for odd primes")
    a %= p
    if a == 0:
        return 0
    residue = pow(a, (p - 1) // 2, p)
    return -1 if residue == p - 1 else residue


def quadratic_status_in_q_sqrt_d(p: int, D: int) -> str:
    """Return ramified/split/inert status of p in Q(sqrt(D))."""
    if p == 2:
        # In the present certificates D is odd and D == 3 mod 4, so 2 is ramified
        # in the quadratic field discriminant.
        return "ramified"
    if D % p == 0:
        return "ramified"
    symbol = legendre_symbol(D, p)
    if symbol == 1:
        return "split"
    if symbol == -1:
        return "inert"
    return "ramified"


def admissibility_witness(p: int, T: List[int]) -> str:
    """Return a simple admissibility witness used in Sawin's Lemma-12 setting."""
    if p % 4 == 1:
        return "p == 1 mod 4"
    if p == 2:
        return "inert in Q(sqrt(5)) in the Kronecker-symbol sense"
    for q in T:
        if p != q and legendre_symbol(q, p) == -1:
            return f"inert in Q(sqrt({q}))"
    return "NO WITNESS FOUND"


def e_value(p: int, T: List[int]) -> int:
    """The ramification factor e(p) used in the exponent formula."""
    return 2 if p == 2 or p in T else 1


def delta_components(cert: Certificate) -> Dict[str, Decimal]:
    """Evaluate numerator, denominator, and delta for one certificate."""
    prod_T = Decimal(product(cert.T))
    R = cert.R
    numerator = (
        (Decimal(1) - Decimal(1) / R).ln()
        + Decimal("0.5") * ((Decimal(2) * PI / Decimal(1).exp()).ln())
        + sum(
            (Decimal(1) / Decimal(4 * e_value(p, cert.T)))
            * Decimal(cert.k[p] + 1).ln()
            for p in cert.S_Q
        )
        - Decimal("0.125") * (Decimal(4) * prod_T).ln()
        - Decimal("0.5") * ((Decimal(4) * prod_T).sqrt().ln()).ln()
    )
    log_product_term = sum(
        Decimal(cert.k[p]) / Decimal(2 * e_value(p, cert.T)) * Decimal(p).ln()
        for p in cert.S_Q
    )
    denominator = (Decimal(2) * R * log_product_term.exp() + Decimal(1)).ln()
    return {
        "numerator": numerator,
        "denominator": denominator,
        "delta": numerator / denominator,
    }


def verify_certificate(cert: Certificate) -> Dict[str, object]:
    """Verify finite checks for a certificate and return a JSON-serializable result."""
    D = product(cert.T)
    statuses = {str(p): quadratic_status_in_q_sqrt_d(p, D) for p in cert.S_Q}
    split_count = sum(1 for status in statuses.values() if status == "split")
    budget_lhs = len(cert.T) + len(cert.S_Q) + split_count + 1
    budget_rhs = (len(cert.T) - 1) ** 2 // 4
    witnesses = {str(p): admissibility_witness(p, cert.T) for p in cert.S_Q}
    components = delta_components(cert)
    checks = {
        "T_all_prime": all(is_prime(q) for q in cert.T),
        "S_Q_all_prime": all(is_prime(p) for p in cert.S_Q),
        "T_all_odd": all(q % 2 == 1 for q in cert.T),
        "T_has_odd_number_3_mod_4": sum(1 for q in cert.T if q % 4 == 3) % 2 == 1,
        "k_positive_for_all_S_Q": all(cert.k.get(p, 0) >= 1 for p in cert.S_Q),
        "R_greater_than_1": cert.R > Decimal(1),
        "no_split_primes_in_S_Q": split_count == 0,
        "budget_satisfied": budget_lhs <= budget_rhs,
        "admissibility_witnesses_found": all(not value.startswith("NO") for value in witnesses.values()),
    }
    return {
        "key": cert.key,
        "name": cert.name,
        "T": cert.T,
        "S_Q": cert.S_Q,
        "k": {str(key): value for key, value in cert.k.items()},
        "R": str(cert.R),
        "D": D,
        "statuses_in_Q_sqrt_D": statuses,
        "admissibility_witnesses": witnesses,
        "split_count": split_count,
        "budget_lhs": budget_lhs,
        "budget_rhs": budget_rhs,
        "numerator": str(components["numerator"]),
        "denominator": str(components["denominator"]),
        "delta": str(components["delta"]),
        "claimed_delta": str(cert.claimed_delta) if cert.claimed_delta is not None else None,
        "checks": checks,
        "passed": all(checks.values()),
    }


def format_text_report(results: List[Dict[str, object]]) -> str:
    """Human-readable report for terminal and TXT output."""
    blocks: List[str] = []
    for result in results:
        blocks.append("=" * 78)
        blocks.append(str(result["name"]))
        blocks.append(f"D = {result['D']}")
        blocks.append(f"R = {result['R']}")
        blocks.append(f"budget = {result['budget_lhs']} <= {result['budget_rhs']}")
        blocks.append(f"delta = {result['delta']}")
        blocks.append(f"passed = {result['passed']}")
        blocks.append("checks:")
        for name, value in result["checks"].items():  # type: ignore[union-attr]
            blocks.append(f"  - {name}: {value}")
        blocks.append("")
    return "\n".join(blocks)


def select_certificates(case: str) -> List[Certificate]:
    """Resolve the CLI --case option."""
    if case == "both":
        return [SAWIN_CERTIFICATE, OPTIMIZED_CERTIFICATE]
    return [CERTIFICATES[case]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        choices=["both", "sawin", "optimized"],
        default="both",
        help="which certificate(s) to verify",
    )
    parser.add_argument(
        "--json-out",
        default="certificate_verification_results.json",
        help="path for JSON output",
    )
    parser.add_argument(
        "--txt-out",
        default="certificate_verification_results.txt",
        help="path for text output",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    results = [verify_certificate(cert) for cert in select_certificates(args.case)]
    text_report = format_text_report(results)
    print(text_report)
    Path(args.json_out).write_text(json.dumps(results, indent=2, sort_keys=True), encoding="utf-8")
    Path(args.txt_out).write_text(text_report + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
