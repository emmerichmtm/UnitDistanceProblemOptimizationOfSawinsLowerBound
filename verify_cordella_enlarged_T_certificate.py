#!/usr/bin/env python3
"""Verify Cordella's enlarged-T addendum certificate from canonical data.

The script expects `cordella_certificate_data.txt` in the same directory and
imports the existing repository verification pipeline from `verify_certificates.py`.
It parses the canonical compact text certificate, constructs a Certificate, and
runs verify_certificate(...) without third-party dependencies.
"""
from __future__ import annotations

import argparse
import json
import re
from decimal import Decimal
from pathlib import Path

from verify_certificates import Certificate, verify_certificate, format_text_report


def parse_certificate(path: Path) -> Certificate:
    text = path.read_text(encoding="utf-8")
    t_match = re.search(r"T \(67 primes\):\n(.+?)\n\nS_Q", text, re.S)
    if not t_match:
        raise ValueError("Could not find T block")
    T = [int(x) for x in re.findall(r"\d+", t_match.group(1))]
    sq_text = text.split("S_Q with multiplicities", 1)[1]
    pairs = [(int(p), int(k)) for p, k in re.findall(r"(\d+)\s*:\s*(\d+)", sq_text)]
    S_Q = [p for p, _ in pairs]
    k = {p: v for p, v in pairs}
    return Certificate(
        key="cordella_enlarged_T_canonical",
        name="Cordella enlarged-T canonical certificate",
        T=T,
        S_Q=S_Q,
        k=k,
        R=Decimal(33066458) / Decimal(1000000),
        claimed_delta=Decimal("0.031185334443176590595661"),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", default="cordella_certificate_data.txt")
    parser.add_argument("--json-out", default="cordella_enlarged_T_verification.json")
    parser.add_argument("--txt-out", default="cordella_enlarged_T_verification.txt")
    args = parser.parse_args()

    cert = parse_certificate(Path(args.data))
    result = verify_certificate(cert)
    Path(args.json_out).write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    Path(args.txt_out).write_text(format_text_report([result]) + "\n", encoding="utf-8")
    print(format_text_report([result]))
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
