# Addendum Python files

These files augment the unit-distance optimization repository with scripts for the enlarged-`T` addendum triggered by Francesco Cordella's personal communication/shared draft-report of 4 June 2026, after publication of v1.

The scripts are intentionally flat and dependency-free.  Place them in the repository root next to `verify_certificates.py`, or use the copy of `verify_certificates.py` included here for a standalone check.

## Files

### Core verification pipeline dependency

- `verify_certificates.py`  
  Original repository verification pipeline. It defines the `Certificate` dataclass, evaluates the exponent formula, and checks primality, parity of `T`, splitting, admissibility witnesses, positivity of multiplicities, and the Golod--Shafarevich budget.

### Canonical Cordella addendum verifier

- `cordella_certificate_data.txt`  
  Canonical compact data for the verified Cordella enlarged-`T` certificate: `#T=67`, `#S_Q=1021`, `R=33.066458`, and `delta=0.031185334443176590595661...`.

- `verify_cordella_enlarged_T_certificate.py`  
  Parses `cordella_certificate_data.txt`, constructs a `Certificate`, and verifies it with `verify_certificate(...)`.

Run:

```bash
python verify_cordella_enlarged_T_certificate.py
```

Expected summary:

```text
passed = True
#T = 67
#S_Q = 1021
split_count = 0
budget = 1089 <= 1089
delta = 0.03118533444317659059566...
```

Outputs:

```text
cordella_enlarged_T_verification.json
cordella_enlarged_T_verification.txt
```

### Verified suitable-budget reruns

- `verify_variable_T_budget_certificates.py`  
  Generates deterministic suitable-budget candidates, verifies each candidate with `verify_certificate(...)`, and writes only passing certificates to the verified result outputs. Non-passing attempted budgets, such as the first-61-odd-prime candidate, are written only to rejected/audit outputs and are not reported as results.

Run:

```bash
python verify_variable_T_budget_certificates.py
```

Outputs:

```text
verified_variable_T_budget_certificates.json
verified_variable_T_budget_certificates.txt
rejected_variable_T_budget_candidates.json
rejected_variable_T_budget_candidates.txt
delta_curve_T_budget_verified_pgfplots.tex
```

The currently verified best rerun is at `#T=67`, with `R=33.06` and `delta=0.031185334443042...`. The canonical Cordella value uses `R=33.066458` and verifies separately with the canonical-data script above.

### Reconstruction and exploratory helper scripts

These are included for provenance and development continuity. They are useful for reproducing how the addendum scripts were developed, but the current report should cite the canonical-data verifier and the verified-only budget script above.

- `addendum_certificate_data.py`  
  Deterministic constructors for enlarged-`T` certificates using the same finite predicates as the verification pipeline.

- `verify_addendum_cordella_certificate.py`  
  Earlier constructor-based Cordella addendum verifier.

- `reproduce_addendum_budget_curve.py`  
  Earlier script for producing the addendum budget curve data.

- `variable_T_budget_es.py`  
  Lightweight exploratory variable-`T` search helper used while developing the addendum plot. Not all exploratory candidates are feasible; always use the verified-only outputs for reported results.

## Recommended commands before committing

```bash
python verify_certificates.py
python verify_cordella_enlarged_T_certificate.py
python verify_variable_T_budget_certificates.py
```

Commit the generated JSON/TXT files if you want the repository to include reproducible verification artifacts.

## Reporting rule

Only certificates with `passed = True` from `verify_certificate(...)` should be reported as results. Rejected candidates may be kept as audit artifacts, but they should not be included in result tables or plots.
