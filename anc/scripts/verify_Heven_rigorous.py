#!/usr/bin/env python3
"""
For each odd Higgs-cubefree k in [k_min, k_max], use the factor cache
+ FactorDB to RIGOROUSLY determine whether m = 2k is in H_even.

Output for each candidate:
  IN_HEVEN:  every prime factor of 2^m+1 verified 3-Higgs
  NOT_HEVEN: some prime factor of 2^m+1 verified non-3-Higgs
  UNKNOWN:   one or more prime factors timed out the Higgs check
             or the factorization is incomplete.
"""
from __future__ import annotations

import argparse
import json
import os
import signal
import sys
from typing import List

from unitary_perfect_research import factor, is_higgs_prime, is_probable_prime


class _T(Exception):
    pass


def _h(signum, frame):
    raise _T()


signal.signal(signal.SIGALRM, _h)

# Hardcoded set of large primes already verified to be non-3-Higgs (via
# their q-1 factorization showing either v_2 > 3 or a non-Higgs sub-prime).
# Used to avoid repeated slow recursive descents.
_KNOWN_NON_HIGGS = {
    # 86-digit factor of 2^566+1: q-1 has v_2(q-1)=142 > 3
    15541351137805832567355695254588151253139260287603415802670284661840802443732044480513,
    # 91-digit factor of 2^1202+1: q-1 has v_2(q-1)=6 > 3
    3571063752373727959434120513220011301363863161567383982128636656862221716975343445451820353,
    # 115-digit factor of 2^758+1: q-1 has v_2(q-1)=190 > 3
    1231312693637327475383720003129487931408741852202045208374953444316525475550246778962497737069179202636323133325313,
}


def safe_higgs(p: int, secs: float) -> bool | None:
    if p in _KNOWN_NON_HIGGS:
        return False
    signal.setitimer(signal.ITIMER_REAL, secs)
    try:
        return is_higgs_prime(p)
    except _T:
        return None
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


def is_higgs_cubefree(k, cache_p_higgs=None):
    if k == 1:
        return True
    fac = factor(k)
    if any(e > 3 for e in fac.values()):
        return False
    for p in fac:
        if cache_p_higgs is not None and p in cache_p_higgs:
            if not cache_p_higgs[p]:
                return False
            continue
        if not is_higgs_prime(p):
            return False
    return True


def classify_m(m: int, cache, deep_secs: float = 60.0) -> tuple[str, int | None]:
    """Returns (status, witness). status in {IN_HEVEN, NOT_HEVEN, UNKNOWN_full,
    UNKNOWN_partial}. witness is the bad prime if NOT_HEVEN.

    The cache only stores prime factors (not composite cofactors). To
    detect partial factorizations, we test whether the product of
    listed prime factors (with multiplicity) equals 2^m + 1. If not,
    the factorization is partial and we cannot rule m IN_HEVEN.
    """
    if m not in cache:
        return "UNKNOWN_no_factor", None
    factor_dict = cache[m]
    # First check the hardcoded large non-Higgs witnesses. This avoids
    # spending timeout budget on other large factors when a decisive
    # witness is already recorded.
    for p in factor_dict:
        if p in _KNOWN_NON_HIGGS:
            return "NOT_HEVEN", p
    # Next: hunt for a non-Higgs prime among the known factors. This is
    # rigorous regardless of completeness.
    timeout_seen = False
    for p in sorted(factor_dict):
        if not is_probable_prime(p):
            continue
        if p < 10**8:
            if not is_higgs_prime(p):
                return "NOT_HEVEN", p
        else:
            r = safe_higgs(p, deep_secs)
            if r is False:
                return "NOT_HEVEN", p
            if r is None:
                timeout_seen = True
    # Now check completeness: known primes^exponents == 2^m+1 ?
    known_product = 1
    for p, e in factor_dict.items():
        known_product *= p ** e
    target = 2**m + 1
    if known_product != target:
        # Partial factorization; we cannot conclude IN_HEVEN.
        return "UNKNOWN_partial_cofactor", None
    if timeout_seen:
        return "UNKNOWN_partial", None
    return "IN_HEVEN", None


def main() -> None:
    if hasattr(sys, "set_int_max_str_digits"):
        sys.set_int_max_str_digits(0)
    parser = argparse.ArgumentParser()
    parser.add_argument("--k-min", type=int, default=1)
    parser.add_argument("--k-max", type=int, default=300)
    parser.add_argument("--deep-secs", type=float, default=60.0)
    parser.add_argument(
        "--cache",
        type=str,
        default=os.environ.get(
            "FACTOR_CACHE",
            os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "factor_cache.json"),
        ),
        help="path to factor cache (default: ../factor_cache.json relative to this script)",
    )
    args = parser.parse_args()

    with open(args.cache) as f:
        cache = json.load(f)
    cache = {int(m): {int(p): int(e) for p, e in d.items()} for m, d in cache.items()}

    in_heven = []
    not_heven = []
    unknown_higgs_timeout = []
    unknown_partial_factor = []
    unknown_nofactor = []
    structurally_excluded = 0

    for k in range(args.k_min, args.k_max + 1, 2):
        if k % 2 == 0:
            continue
        # Structural pre-screen
        if not is_higgs_cubefree(k):
            structurally_excluded += 1
            continue
        m = 2 * k
        status, witness = classify_m(m, cache, args.deep_secs)
        if status == "IN_HEVEN":
            in_heven.append(m)
            print(f"  IN_HEVEN: m={m}", flush=True)
        elif status == "NOT_HEVEN":
            not_heven.append((m, witness))
            print(f"  NOT_HEVEN: m={m} (witness {len(str(witness))}-digit prime)", flush=True)
        elif status == "UNKNOWN_partial":
            unknown_higgs_timeout.append(m)
            print(f"  UNKNOWN_higgs_timeout: m={m}", flush=True)
        elif status == "UNKNOWN_partial_cofactor":
            unknown_partial_factor.append(m)
            print(f"  UNKNOWN_partial_factor: m={m}", flush=True)
        else:
            unknown_nofactor.append(m)
            print(f"  UNKNOWN_no_factor: m={m}", flush=True)

    print()
    print(f"=== Summary k in [{args.k_min}, {args.k_max}] (odd, Higgs-cubefree pre-screen) ===")
    print(f"Structurally excluded (k not Higgs-cubefree):  {structurally_excluded}")
    print(f"IN_HEVEN (full factor + all primes Higgs):      {len(in_heven)}")
    print(f"NOT_HEVEN (verified non-Higgs prime in 2^m+1):  {len(not_heven)}")
    print(f"UNKNOWN_higgs_timeout (slow Higgs check):       {len(unknown_higgs_timeout)}")
    print(f"UNKNOWN_partial_factor (composite cofactor):    {len(unknown_partial_factor)}")
    print(f"UNKNOWN_no_factor (m not in cache):             {len(unknown_nofactor)}")
    print()
    unknowns = unknown_higgs_timeout + unknown_partial_factor + unknown_nofactor
    if unknowns:
        print(f"Certified IN_HEVEN values in [{2*args.k_min}, {2*args.k_max}]: {in_heven}")
        print(f"Unresolved m in this run: {unknowns}")
    else:
        print(f"Verified H_even ∩ [{2*args.k_min}, {2*args.k_max}]: {in_heven}")
    if unknown_higgs_timeout:
        print(f"Higgs-timeout m: {unknown_higgs_timeout}")
    if unknown_partial_factor:
        print(f"Partial-cofactor m: {unknown_partial_factor}")
    if unknown_nofactor:
        print(f"No-factor m: {unknown_nofactor}")


if __name__ == "__main__":
    main()
