# Optimizing an Explicit Unit-Distance Lower-Bound Certificate

This repository contains a self-contained mathematical report and flat-file Python code for a reproducible certificate verification and parameter-optimization study inspired by Sawin's explicit lower bound for the unit-distance problem.

The repository is intentionally flat: all scripts and documents live in the root directory. No third-party Python packages are required.

## Main contribution

The code verifies two finite certificates:

1. **Sawin's published example**, used as a validation target for the verification pipeline.
2. **An optimized candidate**, which passes the implemented arithmetic checks and gives

```text
delta = 0.01517180563721325...
```

Subject to Sawin's explicit criterion being applied exactly as cited in the report, this supports the cautious clean statement

```text
u(n) > n^1.015
```

for arbitrarily large `n`. The sharper decimal value should be checked independently with interval arithmetic and mathematical review.

## Files

- `unit_distance_optimization_report_v15.tex` — LaTeX source for the report.
- `unit_distance_optimization_report_v15.pdf` — compiled report.
- `verify_certificates.py` — verifies Sawin's published certificate and the optimized certificate.
- `optimize_certificates.py` — verifies the built-in examples and runs a deterministic greedy parameter search.
- `coordinate_pipeline_attempt.py` — records why the finite certificate does not by itself yield planar coordinates.
- `explicit_planar_comparator.py` — counts a classical explicit planar comparator based on squared length `5`.

## Quick start

Requires Python 3.11+ for `Decimal.ln()` and `Decimal.exp()`.

Verify both built-in certificates:

```bash
python verify_certificates.py
```

This writes:

```text
certificate_verification_results.json
certificate_verification_results.txt
```

Run the greedy search:

```bash
python optimize_certificates.py --pmax 300 --c 0.015 --r-min 40 --r-max 90 --r-steps 500
```

This writes:

```text
optimization_results.json
```

Record the coordinate-pipeline status:

```bash
python coordinate_pipeline_attempt.py
```

This writes:

```text
coordinate_pipeline_attempt.json
coordinate_pipeline_attempt.txt
```

Count the explicit classical planar comparator:

```bash
python explicit_planar_comparator.py --m 10
```

This writes:

```text
explicit_planar_comparator.json
```

## What the verifier checks

For each certificate, the verifier checks:

- primality of the elements of `T` and `S_Q`,
- the parity condition on primes in `T` congruent to `3 mod 4`,
- splitting/ramification/inertness in `Q(sqrt(D))`, where `D = prod(T)`,
- absence of split primes in `S_Q` for the reported examples,
- the Golod-Shafarevich budget inequality,
- admissibility witnesses for primes in `S_Q`,
- positivity of all multiplicities `k(p)`,
- the exponent formula `delta = numerator / denominator`.

## Coordinate realization

The repository does **not** construct actual planar coordinates for the optimized algebraic candidate. The finite tuple `(T, S_Q, k, R)` is sufficient for the implemented certificate verification but does not specify an explicit finite class-field-tower level, integral basis, ideals, embeddings, normalization, or bounded window. Those data would be needed for a literal coordinate plot in `[-10,10]^2`.

The file `explicit_planar_comparator.py` is included only as a genuine elementary coordinate comparator. It is not the optimized Sawin-style construction.
