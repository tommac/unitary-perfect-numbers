#!/usr/bin/env python3
"""
Deeply re-verify previous claims about all-Higgs m values.

The earlier all_higgs_m_scan.py used a fast Higgs check with short
timeouts, so its "all-Higgs" classification of large m has false
positives: some m it marked all-Higgs in fact contain a large prime r
such that r-1 has a forbidden factor pattern but the recursive check
ran out of time.

This script takes a list of candidate m, runs an unbounded Higgs check
on each prime factor of 2^m+1, and prints VERIFIED ALL-HIGGS / NON-HIGGS
/ TIMEOUT for each. Output is line-buffered so progress is visible.
"""
from __future__ import annotations
import json
import signal
import sys
from typing import List

from unitary_perfect_research import is_higgs_prime, is_probable_prime


class _T(Exception):
    pass


def _h(signum, frame):
    raise _T()


signal.signal(signal.SIGALRM, _h)


def safe_higgs(p: int, secs: float) -> bool | None:
    signal.setitimer(signal.ITIMER_REAL, secs)
    try:
        return is_higgs_prime(p)
    except _T:
        return None
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


def verify(m: int, cache, timeout_per_prime: float = 180.0):
    primes = sorted(cache.get(m, {}))
    if not primes:
        return "no-cache", None
    for p in primes:
        if p < 10**8:
            if not is_higgs_prime(p):
                return "NON-HIGGS", p
            continue
        r = safe_higgs(p, timeout_per_prime)
        if r is False:
            return "NON-HIGGS", p
        if r is None:
            return "TIMEOUT", p
    return "ALL-HIGGS", None


def main() -> None:
    if hasattr(sys, "set_int_max_str_digits"):
        sys.set_int_max_str_digits(0)
    import os
    cache_path = os.environ.get(
        "FACTOR_CACHE",
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "factor_cache.json"),
    )
    with open(cache_path) as f:
        raw = json.load(f)
    cache = {int(m): {int(p): int(e) for p, e in d.items()} for m, d in raw.items()}

    candidates: List[int]
    if len(sys.argv) > 1:
        candidates = [int(arg) for arg in sys.argv[1:]]
    else:
        candidates = [
            107, 117, 122, 183, 223, 335, 347, 419, 457, 503,
            613, 701, 749, 751, 967, 997, 1033, 1051, 1093,
            1123, 1171, 1202,
        ]

    truly_all_higgs: List[int] = []
    non_higgs: List[int] = []
    unknown: List[int] = []
    for m in candidates:
        status, witness = verify(m, cache)
        print(f"m={m}: {status}" + (f"  ({len(str(witness))}-digit prime)" if witness else ""), flush=True)
        if status == "ALL-HIGGS":
            truly_all_higgs.append(m)
        elif status == "NON-HIGGS":
            non_higgs.append(m)
        else:
            unknown.append(m)

    print()
    print(f"Truly all-Higgs (verified): {truly_all_higgs}")
    print(f"Now known non-Higgs:         {non_higgs}")
    print(f"Still unknown:               {unknown}")


if __name__ == "__main__":
    main()
