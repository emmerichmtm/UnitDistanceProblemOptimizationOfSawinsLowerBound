#!/usr/bin/env python3
"""Rudolph-style integer evolutionary search for unit-distance certificates.

This script implements a flat-file, dependency-free evolutionary search over
integer-coded certificate data for Sawin-style lower-bound certificates.

The search uses an integer-native mutation operator inspired by G. Rudolph's
maximum-entropy mutation principle for unbounded integer programming.  In this
implementation, integer steps are sampled from a two-sided geometric law
(a discrete-Laplace distribution)

    P(Z = z) = (1-a)/(1+a) * a**abs(z),  0 < a < 1,

then clipped or repaired where finite bounds are needed.  The real parameter
R is represented rationally as R_num / R_den, so the evolutionary state is fully
integer-coded.

The script compares three baselines:

1. Sawin's published explicit certificate, giving delta ~= 0.0141144.
2. The greedy/optimized certificate from the accompanying note, giving
   delta ~= 0.0151718.
3. The best certificate found by the Rudolph-style integer EA.

No third-party Python packages are required.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

from verify_certificates import (
    Certificate,
    OPTIMIZED_CERTIFICATE,
    SAWIN_CERTIFICATE,
    admissibility_witness,
    delta_components,
    e_value,
    product,
    quadratic_status_in_q_sqrt_d,
    verify_certificate,
)


@dataclass(frozen=True)
class EAConfig:
    """Configuration of the integer evolutionary algorithm."""

    pmax: int = 300
    mu: int = 16
    lambda_: int = 96
    generations: int = 120
    seed: int = 12345
    alpha_index: float = 0.42
    alpha_k: float = 0.35
    alpha_r: float = 0.35
    replace_probability: float = 0.18
    mutate_k_probability: float = 0.35
    mutate_r_probability: float = 0.85
    r_den: int = 1000
    r_min_num: int = 40_000
    r_max_num: int = 90_000
    k_min: int = 1
    k_max: int = 80
    target_size: int = 22


@dataclass
class Individual:
    """Integer-coded individual.

    indices are positions in the admissible prime pool.
    k_values are multiplicities paired with the normalized selected primes.
    r_num encodes R = r_num / r_den.
    """

    indices: List[int]
    k_values: List[int]
    r_num: int
    delta_float: float = float("-inf")
    delta_decimal: str = ""
    valid: bool = False


def primes_up_to(n: int) -> List[int]:
    """Return all primes <= n."""
    if n < 2:
        return []
    sieve = [True] * (n + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(n**0.5) + 1):
        if sieve[i]:
            for j in range(i * i, n + 1, i):
                sieve[j] = False
    return [i for i, flag in enumerate(sieve) if flag]


def candidate_pool(T: List[int], pmax: int, allow_split: bool = False) -> List[int]:
    """Admissible rational primes for S_Q.

    Unlike the greedy search in earlier versions, this pool includes ramified
    primes p in T, matching Sawin's published example.
    """
    D = product(T)
    pool: List[int] = []
    for p in primes_up_to(pmax):
        if admissibility_witness(p, T).startswith("NO"):
            continue
        status = quadratic_status_in_q_sqrt_d(p, D)
        if status == "split" and not allow_split:
            continue
        pool.append(p)
    return pool


def budget_rhs(T: List[int]) -> int:
    return ((len(T) - 1) ** 2) // 4


def two_sided_geometric(rng: random.Random, alpha: float) -> int:
    """Sample an integer from a two-sided geometric distribution.

    P(Z=z) is proportional to alpha**abs(z).  This is a discrete-Laplace law,
    an integer-native mutation distribution with exponential tails.
    """
    if not (0.0 < alpha < 1.0):
        raise ValueError("alpha must be in (0,1)")
    # Magnitude distribution: P(M=0)=(1-a)/(1+a), P(M=m)=2(1-a)/(1+a)a^m.
    p0 = (1.0 - alpha) / (1.0 + alpha)
    if rng.random() < p0:
        return 0
    # Conditional on nonzero, P(M=m | M>0)=(1-alpha)*alpha**(m-1), m>=1.
    u = rng.random()
    magnitude = 1 + int(math.floor(math.log(1.0 - u) / math.log(alpha)))
    return magnitude if rng.random() < 0.5 else -magnitude


def heuristic_k(p: int, c: float = 0.015) -> int:
    """First-order k(p) heuristic used for initialization."""
    return max(1, int(math.floor(1.0 / (2.0 * c * math.log(p)) - 1.0)))


def selected_primes(ind: Individual, pool: Sequence[int], cfg: EAConfig) -> List[int]:
    """Map and repair index list to a sorted unique list of selected primes."""
    n = len(pool)
    unique = sorted({max(0, min(n - 1, idx)) for idx in ind.indices})
    # Repair duplicates by filling with nearby unused indices, then globally.
    if len(unique) < cfg.target_size:
        used = set(unique)
        for idx in list(unique):
            for step in range(1, n):
                for cand in (idx - step, idx + step):
                    if 0 <= cand < n and cand not in used:
                        used.add(cand)
                        unique.append(cand)
                        break
                if len(unique) >= cfg.target_size:
                    break
            if len(unique) >= cfg.target_size:
                break
        if len(unique) < cfg.target_size:
            for cand in range(n):
                if cand not in used:
                    used.add(cand)
                    unique.append(cand)
                    if len(unique) >= cfg.target_size:
                        break
    unique = sorted(unique[: cfg.target_size])
    return [pool[idx] for idx in unique]


def certificate_from_individual(ind: Individual, pool: Sequence[int], cfg: EAConfig, key: str) -> Certificate:
    """Decode an individual into a Certificate."""
    S_Q = selected_primes(ind, pool, cfg)
    # Pair k-values with selected primes after repair/sorting. If needed, extend.
    k_values = list(ind.k_values[: cfg.target_size])
    while len(k_values) < cfg.target_size:
        k_values.append(1)
    k = {p: max(cfg.k_min, min(cfg.k_max, int(k_values[i]))) for i, p in enumerate(S_Q)}
    R = Decimal(ind.r_num) / Decimal(cfg.r_den)
    return Certificate(
        key=key,
        name="Rudolph-style integer EA candidate",
        T=list(SAWIN_CERTIFICATE.T),
        S_Q=S_Q,
        k=k,
        R=R,
    )



def fast_delta_float(cert: Certificate) -> float:
    """Fast double-precision delta used inside the evolutionary loop."""
    prod_T = float(product(cert.T))
    R = float(cert.R)
    numerator = (
        math.log(1.0 - 1.0 / R)
        + 0.5 * math.log(2.0 * math.pi / math.e)
        + sum((1.0 / (4.0 * e_value(p, cert.T))) * math.log(cert.k[p] + 1.0) for p in cert.S_Q)
        - 0.125 * math.log(4.0 * prod_T)
        - 0.5 * math.log(math.log(math.sqrt(4.0 * prod_T)))
    )
    log_product_term = sum(
        cert.k[p] / (2.0 * e_value(p, cert.T)) * math.log(p)
        for p in cert.S_Q
    )
    denominator = math.log(2.0 * R) + log_product_term
    # log(2R exp(log_product_term)+1) is equal to the above up to a negligible
    # correction for the huge values involved in these certificates.
    return numerator / denominator


def fast_passes_side_conditions(cert: Certificate) -> bool:
    """Fast side-condition check for candidates from the prefiltered pool."""
    if len(cert.S_Q) != len(set(cert.S_Q)):
        return False
    if len(cert.S_Q) != 22:
        return False
    if any(cert.k.get(p, 0) < 1 for p in cert.S_Q):
        return False
    if cert.R <= Decimal(1):
        return False
    D = product(cert.T)
    split_count = sum(1 for p in cert.S_Q if quadratic_status_in_q_sqrt_d(p, D) == "split")
    return len(cert.T) + len(cert.S_Q) + split_count + 1 <= ((len(cert.T)-1)**2)//4

def evaluate(ind: Individual, pool: Sequence[int], cfg: EAConfig) -> Individual:
    """Evaluate individual in place and return it."""
    cert = certificate_from_individual(ind, pool, cfg, key="rudolph_ea")
    if not fast_passes_side_conditions(cert):
        ind.delta_float = float("-inf")
        ind.delta_decimal = "-Infinity"
        ind.valid = False
        return ind
    delta = fast_delta_float(cert)
    ind.delta_float = delta
    ind.delta_decimal = f"{delta:.16g}"
    ind.valid = True
    return ind


def individual_from_certificate(cert: Certificate, pool: Sequence[int], cfg: EAConfig) -> Individual:
    """Seed an individual from an existing certificate."""
    index_of = {p: i for i, p in enumerate(pool)}
    indices = [index_of[p] for p in cert.S_Q if p in index_of]
    # If a baseline contains primes outside pool because pmax is too small, repair later.
    k_values = [cert.k[p] for p in cert.S_Q if p in index_of]
    r_num = int((cert.R * Decimal(cfg.r_den)).to_integral_value())
    ind = Individual(indices=indices, k_values=k_values, r_num=r_num)
    return normalize_individual(ind, pool, cfg)


def greedy_seed(pool: Sequence[int], cfg: EAConfig, c: float = 0.015) -> Individual:
    """Construct a deterministic greedy seed from the candidate pool."""
    T = SAWIN_CERTIFICATE.T
    scored: List[Tuple[float, int]] = []
    for idx, p in enumerate(pool):
        k_p = heuristic_k(p, c)
        benefit = (1.0 / (4.0 * e_value(p, T))) * math.log(k_p + 1.0)
        scored.append((benefit, idx))
    scored.sort(reverse=True)
    chosen = sorted(idx for _, idx in scored[: cfg.target_size])
    primes = [pool[idx] for idx in chosen]
    k_values = [heuristic_k(p, c) for p in primes]
    r_num = int(round(66.722 * cfg.r_den))
    return Individual(indices=chosen, k_values=k_values, r_num=r_num)


def normalize_individual(ind: Individual, pool: Sequence[int], cfg: EAConfig) -> Individual:
    """Repair bounds, length, uniqueness, and k/R ranges."""
    n = len(pool)
    indices = [max(0, min(n - 1, int(idx))) for idx in ind.indices]
    if not indices:
        indices = list(range(min(cfg.target_size, n)))
    # Fill to target size before selected_primes repair.
    while len(indices) < cfg.target_size:
        indices.append(indices[-1])
    indices = indices[: cfg.target_size]
    k_values = [max(cfg.k_min, min(cfg.k_max, int(k))) for k in ind.k_values]
    while len(k_values) < cfg.target_size:
        p = pool[indices[len(k_values) % len(indices)]]
        k_values.append(heuristic_k(p))
    k_values = k_values[: cfg.target_size]
    r_num = max(cfg.r_min_num, min(cfg.r_max_num, int(ind.r_num)))
    return Individual(indices=indices, k_values=k_values, r_num=r_num)


def mutate(parent: Individual, pool: Sequence[int], cfg: EAConfig, rng: random.Random) -> Individual:
    """Mutate one parent using integer-native steps."""
    child = Individual(
        indices=list(parent.indices),
        k_values=list(parent.k_values),
        r_num=parent.r_num,
    )
    n = len(pool)
    # Mutate selected prime indices.
    for i in range(cfg.target_size):
        if rng.random() < cfg.replace_probability:
            child.indices[i] = max(
                0,
                min(n - 1, child.indices[i] + two_sided_geometric(rng, cfg.alpha_index)),
            )
    # Occasionally replace one selected index completely; helps escape local neighborhoods.
    if rng.random() < cfg.replace_probability:
        j = rng.randrange(cfg.target_size)
        child.indices[j] = rng.randrange(n)
    # Mutate k-values.
    for i in range(cfg.target_size):
        if rng.random() < cfg.mutate_k_probability:
            child.k_values[i] += two_sided_geometric(rng, cfg.alpha_k)
    # Mutate rational R numerator.
    if rng.random() < cfg.mutate_r_probability:
        # A unit step corresponds to 1/R_den. Multiply for practically useful R motion.
        child.r_num += 25 * two_sided_geometric(rng, cfg.alpha_r)
    return normalize_individual(child, pool, cfg)


def run_ea(cfg: EAConfig) -> Tuple[Individual, List[Dict[str, object]], List[int]]:
    """Run a (mu+lambda) integer EA and return best individual and history."""
    rng = random.Random(cfg.seed)
    pool = candidate_pool(SAWIN_CERTIFICATE.T, cfg.pmax, allow_split=False)
    if len(pool) < cfg.target_size:
        raise ValueError("candidate pool is too small; increase --pmax")
    population: List[Individual] = [
        individual_from_certificate(SAWIN_CERTIFICATE, pool, cfg),
        individual_from_certificate(OPTIMIZED_CERTIFICATE, pool, cfg),
        greedy_seed(pool, cfg),
    ]
    # Add randomized perturbations of the greedy seed.
    while len(population) < cfg.mu:
        population.append(mutate(greedy_seed(pool, cfg), pool, cfg, rng))
    population = [evaluate(ind, pool, cfg) for ind in population]
    population.sort(key=lambda ind: ind.delta_float, reverse=True)
    history: List[Dict[str, object]] = []
    for generation in range(cfg.generations + 1):
        best = population[0]
        history.append(
            {
                "generation": generation,
                "best_delta": best.delta_decimal,
                "best_R": str(Decimal(best.r_num) / Decimal(cfg.r_den)),
                "valid_population": sum(1 for ind in population if ind.valid),
            }
        )
        if generation == cfg.generations:
            break
        offspring: List[Individual] = []
        for _ in range(cfg.lambda_):
            parent = rng.choice(population[: cfg.mu])
            offspring.append(evaluate(mutate(parent, pool, cfg, rng), pool, cfg))
        population = sorted(population + offspring, key=lambda ind: ind.delta_float, reverse=True)[: cfg.mu]
    return population[0], history, list(pool)


def result_summary(best: Individual, pool: Sequence[int], cfg: EAConfig) -> Dict[str, object]:
    """Create a JSON-serializable summary for the best EA candidate."""
    cert = certificate_from_individual(best, pool, cfg, key="rudolph_ea")
    verification = verify_certificate(cert)
    comps = delta_components(cert)
    return {
        "method": "Rudolph-style integer EA with two-sided geometric mutation",
        "R_representation": {
            "R_num": best.r_num,
            "R_den": cfg.r_den,
            "R": str(cert.R),
        },
        "certificate": {
            "T": cert.T,
            "S_Q": cert.S_Q,
            "k": {str(p): cert.k[p] for p in cert.S_Q},
        },
        "verification": verification,
        "delta_components": {k: str(v) for k, v in comps.items()},
    }


def baseline_summary(cert: Certificate) -> Dict[str, object]:
    comps = delta_components(cert)
    verification = verify_certificate(cert)
    return {
        "name": cert.name,
        "key": cert.key,
        "R": str(cert.R),
        "S_Q": cert.S_Q,
        "k": {str(p): cert.k[p] for p in cert.S_Q},
        "delta": str(comps["delta"]),
        "passes_all_checks": verification["passed"],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pmax", type=int, default=300, help="largest candidate prime considered")
    parser.add_argument("--mu", type=int, default=16, help="number of parents retained")
    parser.add_argument("--lambda", dest="lambda_", type=int, default=96, help="offspring per generation")
    parser.add_argument("--generations", type=int, default=120, help="number of generations")
    parser.add_argument("--seed", type=int, default=12345, help="random seed")
    parser.add_argument("--r-den", type=int, default=1000, help="denominator for rational R")
    parser.add_argument("--r-min", type=float, default=40.0, help="minimum R")
    parser.add_argument("--r-max", type=float, default=90.0, help="maximum R")
    parser.add_argument("--json-out", default="rudolph_integer_ea_results.json", help="JSON output path")
    parser.add_argument("--txt-out", default="rudolph_integer_ea_results.txt", help="text output path")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cfg = EAConfig(
        pmax=args.pmax,
        mu=args.mu,
        lambda_=args.lambda_,
        generations=args.generations,
        seed=args.seed,
        r_den=args.r_den,
        r_min_num=int(round(args.r_min * args.r_den)),
        r_max_num=int(round(args.r_max * args.r_den)),
    )
    best, history, pool = run_ea(cfg)
    ea_summary = result_summary(best, pool, cfg)
    sawin = baseline_summary(SAWIN_CERTIFICATE)
    greedy = baseline_summary(OPTIMIZED_CERTIFICATE)
    output = {
        "config": cfg.__dict__,
        "candidate_pool_size": len(pool),
        "candidate_pool": pool,
        "baselines": {
            "sawin_published": sawin,
            "greedy_optimized": greedy,
        },
        "rudolph_integer_ea": ea_summary,
        "history": history,
    }
    Path(args.json_out).write_text(json.dumps(output, indent=2), encoding="utf-8")
    lines = [
        "Rudolph-style integer EA comparison",
        "====================================",
        f"Sawin published delta:   {sawin['delta']}",
        f"Greedy/optimized delta:  {greedy['delta']}",
        f"Integer-EA best delta:   {ea_summary['delta_components']['delta']}",
        f"Integer-EA R:            {ea_summary['R_representation']['R_num']}/{ea_summary['R_representation']['R_den']} = {ea_summary['R_representation']['R']}",
        f"Integer-EA passes checks: {ea_summary['verification']['passed']}",
        "",
        "Best Integer-EA S_Q:",
        str(ea_summary["certificate"]["S_Q"]),
        "",
        "Best Integer-EA k values:",
        str(ea_summary["certificate"]["k"]),
        "",
        "Interpretation:",
        "The integer EA is initialized with Sawin's certificate, the current greedy/optimized certificate,",
        "and a deterministic greedy seed. Any improvement is therefore measured against both baselines.",
    ]
    Path(args.txt_out).write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
