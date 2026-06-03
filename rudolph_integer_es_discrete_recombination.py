#!/usr/bin/env python3
"""Rudolph-style integer ES with discrete recombination.

This script tests a two-parent discrete recombination variant of the
Rudolph-style integer evolutionary search used for Sawin-style unit-distance
certificates.  Object variables remain integer-coded: selected-prime indices,
integer multiplicities k(p), and the rational numerator R_num for R=R_num/R_den.

Discrete recombination is componentwise: each offspring variable is inherited
from one of two selected parents with probability 1/2.  After recombination,
the same integer-native two-sided geometric mutation operator is applied.

No third-party Python packages are required.
"""

from __future__ import annotations

import argparse
import json
import random
from decimal import Decimal
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from verify_certificates import OPTIMIZED_CERTIFICATE, SAWIN_CERTIFICATE, delta_components, verify_certificate
from rudolph_integer_ea import (
    EAConfig,
    Individual,
    baseline_summary,
    candidate_pool,
    certificate_from_individual,
    evaluate,
    greedy_seed,
    individual_from_certificate,
    mutate,
    normalize_individual,
    result_summary,
)


def discrete_recombine(a: Individual, b: Individual, cfg: EAConfig, rng: random.Random) -> Individual:
    """Componentwise discrete recombination of two integer-coded parents."""
    indices: List[int] = []
    k_values: List[int] = []
    for i in range(cfg.target_size):
        ai = a.indices[i % len(a.indices)]
        bi = b.indices[i % len(b.indices)]
        ak = a.k_values[i % len(a.k_values)]
        bk = b.k_values[i % len(b.k_values)]
        indices.append(ai if rng.random() < 0.5 else bi)
        k_values.append(ak if rng.random() < 0.5 else bk)
    r_num = a.r_num if rng.random() < 0.5 else b.r_num
    return Individual(indices=indices, k_values=k_values, r_num=r_num)


def run_recombination_es(cfg: EAConfig) -> Tuple[Individual, List[Dict[str, object]], List[int]]:
    """Run a (mu+lambda) ES with discrete recombination followed by integer mutation."""
    rng = random.Random(cfg.seed)
    pool = candidate_pool(SAWIN_CERTIFICATE.T, cfg.pmax, allow_split=False)
    if len(pool) < cfg.target_size:
        raise ValueError("candidate pool is too small; increase --pmax")

    population: List[Individual] = [
        individual_from_certificate(SAWIN_CERTIFICATE, pool, cfg),
        individual_from_certificate(OPTIMIZED_CERTIFICATE, pool, cfg),
        greedy_seed(pool, cfg),
    ]
    # Also seed from the best simple Rudolph-style integer ES candidate from v21.
    rudolph_v21_SQ = [2, 3, 5, 47, 71, 79, 97, 101, 107, 109, 139, 151, 163, 167, 179, 191, 211, 223, 239, 241, 251, 263]
    rudolph_v21_k = {2:49, 3:29, 5:19, 47:7, 71:7, 79:6, 97:6, 101:6, 107:6, 109:6, 139:6, 151:5, 163:5, 167:6, 179:5, 191:5, 211:5, 223:5, 239:5, 241:5, 251:5, 263:5}
    from verify_certificates import Certificate
    rudolph_v21_cert = Certificate(
        key="rudolph_v21",
        name="Simple Rudolph-style integer ES v21 certificate",
        T=list(SAWIN_CERTIFICATE.T),
        S_Q=rudolph_v21_SQ,
        k=rudolph_v21_k,
        R=Decimal(6672416) / Decimal(100000),
    )
    population.append(individual_from_certificate(rudolph_v21_cert, pool, cfg))

    while len(population) < cfg.mu:
        population.append(mutate(greedy_seed(pool, cfg), pool, cfg, rng))
    population = [evaluate(ind, pool, cfg) for ind in population]
    population.sort(key=lambda ind: ind.delta_float, reverse=True)

    history: List[Dict[str, object]] = []
    for generation in range(cfg.generations + 1):
        best = population[0]
        if generation % 10 == 0 or generation == cfg.generations:
            history.append({
                "generation": generation,
                "best_delta": best.delta_decimal,
                "best_R": str(Decimal(best.r_num) / Decimal(cfg.r_den)),
                "valid_population": sum(1 for ind in population if ind.valid),
            })
        if generation == cfg.generations:
            break
        offspring: List[Individual] = []
        elite = population[: cfg.mu]
        for _ in range(cfg.lambda_):
            p1 = rng.choice(elite)
            p2 = rng.choice(elite)
            child = normalize_individual(discrete_recombine(p1, p2, cfg, rng), pool, cfg)
            child = mutate(child, pool, cfg, rng)
            offspring.append(evaluate(child, pool, cfg))
        population = sorted(population + offspring, key=lambda ind: ind.delta_float, reverse=True)[: cfg.mu]
    return population[0], history, list(pool)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pmax", type=int, default=300)
    parser.add_argument("--mu", type=int, default=24)
    parser.add_argument("--lambda", dest="lambda_", type=int, default=144)
    parser.add_argument("--generations", type=int, default=160)
    parser.add_argument("--seeds", type=int, nargs="*", default=[12345, 54321, 98765, 20260601])
    parser.add_argument("--r-den", type=int, default=100000)
    parser.add_argument("--r-min", type=float, default=40.0)
    parser.add_argument("--r-max", type=float, default=90.0)
    parser.add_argument("--json-out", default="rudolph_integer_es_discrete_recombination_results.json")
    parser.add_argument("--txt-out", default="rudolph_integer_es_discrete_recombination_results.txt")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    sawin = baseline_summary(SAWIN_CERTIFICATE)
    greedy = baseline_summary(OPTIMIZED_CERTIFICATE)
    simple_best_delta = Decimal("0.0152616610684192942425450138731409223894393631967856067120539712741573810118422083384893589")
    simple_best_R = "66.72416"

    runs = []
    best_summary = None
    best_delta = Decimal("-Infinity")
    for seed in args.seeds:
        cfg = EAConfig(
            pmax=args.pmax,
            mu=args.mu,
            lambda_=args.lambda_,
            generations=args.generations,
            seed=seed,
            r_den=args.r_den,
            r_min_num=int(round(args.r_min * args.r_den)),
            r_max_num=int(round(args.r_max * args.r_den)),
        )
        best, history, pool = run_recombination_es(cfg)
        summary = result_summary(best, pool, cfg)
        d = Decimal(summary["delta_components"]["delta"])
        runs.append({"seed": seed, "config": cfg.__dict__, "summary": summary, "history": history})
        if d > best_delta:
            best_delta = d
            best_summary = summary

    assert best_summary is not None
    improvement_over_simple = best_delta - simple_best_delta
    output = {
        "baselines": {
            "sawin_published": sawin,
            "greedy_optimized": greedy,
            "simple_rudolph_integer_es": {
                "delta": str(simple_best_delta),
                "R": simple_best_R,
            },
        },
        "method": "Rudolph-style integer ES with discrete recombination",
        "best_discrete_recombination": best_summary,
        "improvement_over_simple_rudolph_delta": str(improvement_over_simple),
        "improved_over_simple_rudolph": bool(improvement_over_simple > 0),
        "runs": runs,
        "interpretation": (
            "Discrete recombination was tested as a two-parent componentwise recombination operator. "
            "The best result should be compared with the simple Rudolph-style integer ES baseline."
        ),
    }
    Path(args.json_out).write_text(json.dumps(output, indent=2), encoding="utf-8")
    lines = [
        "Rudolph-style integer ES with discrete recombination",
        "====================================================",
        f"Sawin published delta:          {sawin['delta']}",
        f"Greedy/optimized delta:         {greedy['delta']}",
        f"Simple Rudolph integer ES:      {simple_best_delta}",
        f"Discrete recombination best:    {best_delta}",
        f"Improvement over simple ES:      {improvement_over_simple}",
        f"Improved over simple ES?         {improvement_over_simple > 0}",
        "",
        "Best discrete-recombination certificate:",
        f"R = {best_summary['R_representation']['R']}",
        f"S_Q = {best_summary['certificate']['S_Q']}",
        f"k = {best_summary['certificate']['k']}",
        "",
        "Interpretation:",
        "Discrete recombination was tested as an additional operator. If it does not improve",
        "the simple Rudolph-style integer ES baseline, it should be reported as a negative",
        "algorithmic experiment rather than as the main result.",
    ]
    Path(args.txt_out).write_text("\n".join(lines)+"\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
