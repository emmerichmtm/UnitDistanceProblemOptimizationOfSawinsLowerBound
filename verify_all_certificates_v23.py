#!/usr/bin/env python3
"""Verify all certificate levels reported in the v23 report.

This script is intentionally flat and dependency-free.  It imports the reusable
verification routines from verify_certificates.py and checks four certificates:
Sawin's published example, the greedy optimized certificate, the simple
Rudolph-style integer ES certificate, and the Rudolph-style integer ES with
discrete recombination.
"""
from __future__ import annotations

import json
from decimal import Decimal
from pathlib import Path

from verify_certificates import Certificate, verify_certificate, format_text_report, SAWIN_CERTIFICATE, OPTIMIZED_CERTIFICATE

SIMPLE_RUDOLPH_CERTIFICATE = Certificate(
    key="simple_rudolph_integer_es",
    name="Rudolph-style integer ES certificate",
    T=[3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43],
    S_Q=[2, 3, 5, 47, 71, 79, 97, 101, 107, 109, 139, 151, 163, 167, 179, 191, 211, 223, 239, 241, 251, 263],
    k={2: 49, 3: 29, 5: 19, 47: 7, 71: 7, 79: 6, 97: 6, 101: 6, 107: 6, 109: 6, 139: 6, 151: 5, 163: 5, 167: 6, 179: 5, 191: 5, 211: 5, 223: 5, 239: 5, 241: 5, 251: 5, 263: 5},
    R=Decimal("66.72416"),
    claimed_delta=Decimal("0.015261661068419294"),
)

DISCRETE_RECOMBINATION_CERTIFICATE = Certificate(
    key="rudolph_integer_es_discrete_recombination",
    name="Rudolph-style integer ES with discrete recombination",
    T=[3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43],
    S_Q=[2, 3, 5, 47, 71, 79, 97, 101, 107, 109, 139, 151, 163, 167, 179, 191, 211, 223, 239, 241, 251, 263],
    k={2: 49, 3: 29, 5: 19, 47: 8, 71: 7, 79: 7, 97: 6, 101: 6, 107: 6, 109: 6, 139: 6, 151: 6, 163: 5, 167: 6, 179: 5, 191: 5, 211: 5, 223: 5, 239: 5, 241: 5, 251: 5, 263: 5},
    R=Decimal("66.72416"),
    claimed_delta=Decimal("0.015262868817007192"),
)

CERTIFICATES = [SAWIN_CERTIFICATE, OPTIMIZED_CERTIFICATE, SIMPLE_RUDOLPH_CERTIFICATE, DISCRETE_RECOMBINATION_CERTIFICATE]


def main() -> None:
    out_json = Path("certificate_verification_all_v23.json")
    out_txt = Path("certificate_verification_all_v23.txt")
    results = [verify_certificate(cert) for cert in CERTIFICATES]
    out_json.write_text(json.dumps(results, indent=2), encoding="utf-8")
    out_txt.write_text(format_text_report(results), encoding="utf-8")
    print(format_text_report(results))


if __name__ == "__main__":
    main()
