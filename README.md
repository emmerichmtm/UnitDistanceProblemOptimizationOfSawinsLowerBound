# Optimizing Explicit Unit-Distance Lower-Bound Certificates

This repository contains the arXiv report and flat-file Python code for a reproducible certificate-verification and integer-optimization study inspired by Sawin's explicit lower bound for the unit-distance problem.

The repository is intentionally flat: all scripts and documents live in the root directory. No third-party Python packages are required.

## Main contribution

The code verifies four certificate levels:

1. **Sawin's published example**, used as a validation target for the verification pipeline.
2. **A greedy optimized certificate**, obtained by a deterministic budget heuristic.
3. **A Tailored Integer Evolution Strategy certificate**, using integer-valued mutation in the integer-programming/evolution-strategy tradition.
4. **A Tailored Integer Evolution Strategy certificate with discrete recombination**, the best certificate reported in the paper.

Subject to Sawin's explicit criterion being applied exactly as cited in the report, the best verified certificate supports the cautious statement

```text
u(n) > n^1.0152
```

for arbitrarily large `n`. The sharper decimal value should be checked independently with interval arithmetic and mathematical review.

## Files

- `2606.03419v1.pdf` -- the accompanying arXiv report.
- `verify_certificates.py` -- verifies Sawin's published certificate and the greedy optimized certificate.
- `verify_all_certificates_integer_evolution_strategy.py` -- verifies all four certificate levels, including the integer-evolution-strategy certificates.
- `optimize_certificates.py` -- verifies the built-in examples and runs the deterministic greedy parameter optimization.
- `rudolph_integer_ea.py` -- support module implementing the integer-coded evolutionary algorithm components.
- `rudolph_integer_es_discrete_recombination.py` -- runs the Tailored Integer Evolution Strategy with discrete recombination.

## Quick start

Requires Python 3.11+ for `Decimal.ln()` and `Decimal.exp()`.

Verify the two built-in baseline certificates:

```bash
python verify_certificates.py
```

This writes:

```text
certificate_verification_results.json
certificate_verification_results.txt
```

Run the full four-certificate verification:

```bash
python verify_all_certificates_integer_evolution_strategy.py
```

This verifies:

```text
Sawin published example
Greedy optimized certificate
Tailored Integer Evolution Strategy certificate
Tailored Integer Evolution Strategy with discrete recombination
```

Run the greedy optimization strategy:

```bash
python optimize_certificates.py --pmax 300 --c 0.015 --r-min 40 --r-max 90 --r-steps 500
```

This writes:

```text
optimization_results.json
```

Run the Tailored Integer Evolution Strategy with discrete recombination:

```bash
python rudolph_integer_es_discrete_recombination.py
```

This writes:

```text
rudolph_integer_es_discrete_recombination_results.json
rudolph_integer_es_discrete_recombination_results.txt
```

With the default seeds, the reported run reproduces the best certificate value

```text
delta ≈ 0.0152628688170072
```

## What the verifier checks

For each certificate, the verifier checks:

- primality of the elements of `T` and `S_Q`,
- the parity condition on primes in `T` congruent to `3 mod 4`,
- splitting/ramification/inertness in `Q(sqrt(D))`, where `D = prod(T)`,
- absence of split primes in `S_Q` for the reported examples,
- the Golod--Shafarevich budget inequality,
- admissibility witnesses for primes in `S_Q`,
- positivity of all multiplicities `k(p)`,
- the exponent formula `delta = numerator / denominator`.

## Coordinate realization

The repository does **not** construct actual planar coordinates for the optimized algebraic certificate. The finite tuple `(T, S_Q, k, R)` is sufficient for the implemented certificate verification but does not specify an explicit finite class-field-tower level, integral basis, ideals, embeddings, normalization, or bounded window. Those data would be needed for a literal coordinate plot.

## Reproducibility note

The scripts use deterministic seeds where randomized optimization is involved. The default Tailored Integer Evolution Strategy run performs a lightweight computational experiment intended to be reproducible on standard hardware.
