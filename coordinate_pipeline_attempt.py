#!/usr/bin/env python3
"""Document the coordinate-generation status for the optimized certificate.

The optimized finite certificate verifies arithmetic side conditions and an
exponent formula. It does not specify the additional number-field data needed to
enumerate actual planar coordinates. This script records that status in JSON and
TXT form so the repository contains a reproducible scope statement.

No third-party Python packages are required.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from verify_certificates import OPTIMIZED_CERTIFICATE, verify_certificate


MISSING_DATA = [
    "an explicit finite level L of the relevant Golod-Shafarevich tower",
    "a defining polynomial or integral basis for L over Q",
    "the associated CM extension K and relevant fractional ideals",
    "explicit embeddings/projections into the plane",
    "a normalization and bounded window for extracting a finite point set",
    "a certified unit-distance edge enumeration for the extracted coordinates",
]


def build_status() -> dict:
    verification = verify_certificate(OPTIMIZED_CERTIFICATE)
    return {
        "certificate": OPTIMIZED_CERTIFICATE.name,
        "verification_passed": verification["passed"],
        "delta": verification["delta"],
        "coordinate_status": "BLOCKED_BEFORE_COORDINATES",
        "reason": (
            "The finite tuple (T, S_Q, k, R) verifies Sawin-style certificate "
            "conditions but does not specify a computable finite tower level, "
            "basis, ideals, embeddings, normalization, or coordinate window."
        ),
        "missing_data": MISSING_DATA,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-out", default="coordinate_pipeline_attempt.json")
    parser.add_argument("--txt-out", default="coordinate_pipeline_attempt.txt")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    status = build_status()
    text = [
        f"certificate: {status['certificate']}",
        f"verification_passed: {status['verification_passed']}",
        f"delta: {status['delta']}",
        f"coordinate_status: {status['coordinate_status']}",
        "",
        "reason:",
        str(status["reason"]),
        "",
        "missing data:",
    ]
    text.extend(f"  - {item}" for item in MISSING_DATA)
    Path(args.json_out).write_text(json.dumps(status, indent=2, sort_keys=True), encoding="utf-8")
    Path(args.txt_out).write_text("\n".join(text) + "\n", encoding="utf-8")
    print("\n".join(text))


if __name__ == "__main__":
    main()
