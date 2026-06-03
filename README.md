# Optimizing an Explicit Unit-Distance Lower-Bound Certificate

This repository contains a self-contained mathematical report and flat-file Python code for a reproducible certificate verification and parameter-optimization study inspired by Sawin's explicit lower bound for the unit-distance problem.

The repository is intentionally flat: all scripts and documents live in the root directory. No third-party Python packages are required.

## Main contribution

The code verifies four finite certificates:

1. **Sawin's published example**, used as a validation target for the verification pipeline.
2. **An optimized candidate**, based on the greedy pipeline
3. **Tailored Integer ES Optimization" 

Subject to Sawin's explicit criterion being applied exactly as cited in the report, this supports the cautious clean statement

```text
u(n) > n^1.015
```

for arbitrarily large `n`. The sharper decimal value should be checked independently with interval arithmetic and mathematical review.

## Files

- `2606.03419v1.pdf` -- the accompanying Arxiv report
-  `verify_certificates.py` -- verifies Sawin's published certificate and the optimized certificate.
- `optimize_certificates.py` -- verifies the built-in examples and runs a deterministic greedy parameter search.
- `coordinate_pipeline_attempt.py` -- records why the finite certificate does not by itself yield planar coordinates.
- `explicit_planar_comparator.py` -- counts a classical explicit planar comparator based on squared length `5`.

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

## Rudolph-style integer evolutionary search

A separate flat-file script implements an integer-coded evolutionary search inspired by Rudolph's integer-programming EA:

```bash
python rudolph_integer_ea.py --generations 300 --mu 24 --lambda 160 --pmax 500 --seed 20260601 --r-den 100000
```

The mutation operator is integer-native, using a two-sided geometric / discrete-Laplace step distribution. The real parameter is represented rationally as `R_num / R_den`, so the full evolutionary state is integer-coded.

The bundled run compares three certificates:

- Sawin published example: `delta ≈ 0.0141144286784982`
- previous greedy/optimized candidate: `delta ≈ 0.0151718056372133`
- Rudolph-style integer EA candidate: `delta ≈ 0.0152616610684193`

See `RUDOLPH_INTEGER_EA.md`, `rudolph_integer_ea_results.txt`, and `rudolph_integer_ea_results.json` for details.


## Version note

This v21 bundle intentionally uses the simpler Rudolph-style integer evolutionary search that produced the certificate with

```text
rho delta = 0.0152616610684193...
```

The later self-adaptive mutation-rate experiment is not used in the report, because it did not provide a meaningful methodological improvement for the paper.

## Discrete recombination variant

The file `rudolph_integer_es_discrete_recombination.py` tests a two-parent discrete-recombination variant of the Rudolph-style integer ES. Run:

```bash
python rudolph_integer_es_discrete_recombination.py
```

The recorded run found a small improvement over the non-recombining integer ES:

- simple Rudolph-style integer ES: `delta ≈ 0.0152616610684193`
- discrete-recombination variant: `delta ≈ 0.0152628688170072`

A self-adaptive mutation-rate variant was also tested but was not retained as the main result, since it appeared to converge prematurely and did not improve the certificate.


## Additional v23 verification

Run the full four-certificate verification with:

```bash
python verify_all_certificates_v23.py
```

This checks Sawin, the greedy certificate, the simple Rudolph-style integer ES certificate, and the discrete-recombination variant.
