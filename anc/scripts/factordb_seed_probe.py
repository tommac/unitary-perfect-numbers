#!/usr/bin/env python3
"""
Query FactorDB for factorizations of 2^a + 1 seeds.

This is an optional network-backed helper for extending the seed frontier when
the local Pollard-Rho/Pollard-p-1 stack times out. It prints only whether a
fully factored seed survives the same local filters used by seed_factor_probe.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
import urllib.parse
import urllib.request
from typing import Dict, List, Tuple

from unitary_perfect_research import (
    SearchConfig,
    format_factor,
    is_higgs_prime,
    is_probable_prime,
    seed_divisors_allowed,
    zsigmondy_exponent_allowed,
)


def factordb_query(n: int, timeout: float) -> dict:
    url = "https://factordb.com/api?query=" + urllib.parse.quote(str(n))
    with urllib.request.urlopen(url, timeout=timeout) as response:
        return json.loads(response.read().decode())


def factors_from_response(data: dict) -> Dict[int, int] | None:
    if data.get("status") != "FF":
        return None
    factors: Dict[int, int] = {}
    for p_text, exponent in data.get("factors", []):
        factors[int(p_text)] = int(exponent)
    return dict(sorted(factors.items()))


def partial_rejection_reasons(
    data: dict,
    max_odd_exp: int,
    max_prime: int,
    max_component: int,
    higgs_filter: bool,
) -> List[str]:
    """Return conclusive rejection reasons from known prime factors in a CF response."""
    reasons: List[str] = []
    for p_text, exponent in data.get("factors", []):
        p = int(p_text)
        exponent = int(exponent)
        if not is_probable_prime(p):
            continue
        if p > max_prime:
            reasons.append("partial-prime-cap")
        if exponent > max_odd_exp:
            reasons.append("partial-seed-exponent")
        if p**exponent > max_component:
            reasons.append("partial-component-cap")
        if (
            higgs_filter
            and p <= max_prime
            and p**exponent <= max_component
            and not is_higgs_prime(p)
        ):
            reasons.append("partial-non-higgs")
    return sorted(set(reasons))


def survives(
    factors: Dict[int, int],
    max_odd_exp: int,
    max_prime: int,
    max_component: int,
    higgs_filter: bool,
) -> Tuple[bool, List[str]]:
    reasons: List[str] = []
    if any(p > max_prime for p in factors):
        reasons.append("prime-cap")
    if any(e > max_odd_exp for e in factors.values()):
        reasons.append("seed-exponent")
    if any(p**e > max_component for p, e in factors.items()):
        reasons.append("component-cap")
    if higgs_filter and not reasons and any(not is_higgs_prime(p) for p in factors):
        reasons.append("non-higgs")
    return not reasons, reasons


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("exponents", nargs="*", type=int)
    parser.add_argument("--min-even-exp", type=int)
    parser.add_argument("--max-even-exp", type=int)
    parser.add_argument("--timeout", type=float, default=30.0)
    parser.add_argument("--max-odd-exp", type=int, default=64)
    parser.add_argument("--max-prime", type=int, default=10**12)
    parser.add_argument("--max-component", type=int, default=10**18)
    parser.add_argument("--higgs-filter", action="store_true")
    parser.add_argument("--zsigmondy-exponent-filter", action="store_true")
    parser.add_argument("--seed-divisor-filter", action="store_true")
    parser.add_argument("--quiet-failures", action="store_true")
    parser.add_argument("--progress-interval", type=int, default=25)
    args = parser.parse_args()

    exponents = args.exponents
    if not exponents:
        if args.min_even_exp is None or args.max_even_exp is None:
            parser.error("provide exponents or both --min-even-exp and --max-even-exp")
        exponents = list(range(args.min_even_exp, args.max_even_exp + 1))

    cfg = SearchConfig(
        min_even_exp=min(exponents),
        max_even_exp=max(exponents),
        max_odd_exp=args.max_odd_exp,
        max_bases=1,
        max_prime=args.max_prime,
        max_component=args.max_component,
        max_solutions=1,
        higgs_filter=args.higgs_filter,
        zsigmondy_exponent_filter=args.zsigmondy_exponent_filter,
        seed_divisor_filter=args.seed_divisor_filter,
    )
    survivors = []
    incomplete = []
    reason_counts: Counter[str] = Counter()
    checked = 0
    for a in exponents:
        if not zsigmondy_exponent_allowed(2, a, cfg):
            reason_counts["skipped-zsigmondy"] += 1
            if not args.quiet_failures:
                print(f"{a}: skipped=zsigmondy-exponent-filter", flush=True)
            continue
        if args.seed_divisor_filter and not seed_divisors_allowed(a, cfg):
            reason_counts["skipped-seed-divisor"] += 1
            if not args.quiet_failures:
                print(f"{a}: skipped=seed-divisor-filter", flush=True)
            continue
        try:
            data = factordb_query(2**a + 1, args.timeout)
        except Exception as exc:
            incomplete.append((a, f"error:{exc.__class__.__name__}"))
            print(f"{a}: incomplete=error:{exc.__class__.__name__}", flush=True)
            continue
        factors = factors_from_response(data)
        if factors is None:
            reasons = partial_rejection_reasons(
                data,
                args.max_odd_exp,
                args.max_prime,
                args.max_component,
                args.higgs_filter,
            )
            if reasons:
                for reason in reasons:
                    reason_counts[reason] += 1
                if not args.quiet_failures:
                    print(
                        f"{a}: partial-reject={data.get('status')} "
                        f"reasons={','.join(reasons)}",
                        flush=True,
                    )
            else:
                incomplete.append((a, data.get("status")))
                print(f"{a}: incomplete={data.get('status')}", flush=True)
            continue
        checked += 1
        ok, reasons = survives(
            factors,
            args.max_odd_exp,
            args.max_prime,
            args.max_component,
            args.higgs_filter,
        )
        for reason in reasons or ["survivor"]:
            reason_counts[reason] += 1
        if ok or not args.quiet_failures:
            print(f"{a}: status=FF survives={ok} reasons={','.join(reasons) or '-'}", flush=True)
            print(f"  {format_factor(factors)}", flush=True)
        if ok:
            survivors.append((a, factors))
        elif (
            args.quiet_failures
            and args.progress_interval > 0
            and checked % args.progress_interval == 0
        ):
            print(f"checked_ff: {checked}", flush=True)

    print(f"survivors: {len(survivors)}", flush=True)
    print(f"incomplete: {len(incomplete)}", flush=True)
    for a, status in incomplete:
        print(f"{a}: {status}", flush=True)
    if reason_counts:
        print("reason_counts:", flush=True)
        for reason, count in sorted(reason_counts.items()):
            print(f"  {reason}: {count}", flush=True)


if __name__ == "__main__":
    main()
