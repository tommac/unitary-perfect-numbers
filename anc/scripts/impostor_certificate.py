#!/usr/bin/env python3
"""
Bounded certificate generator for the five impostor SCC kernels.

NOTE (version 3 of the paper): point 1 below assumes Hypothesis H3, that every
prime factor of a unitary perfect number is 3-Higgs, which is unproven. The
v2-overshoot filter (point 2) is invalid, because v_2(p^e+1) is not monotone in
e, and is withdrawn; the ten candidates it was used for are settled by the exact
search in ../search/. The script is kept as the record of the version 2 run.

For each impostor kernel K with seed congruence a == r mod N and each a
in [r, max_a] within that class, certify that no unitary perfect number
has 2^a || n and K as forced odd components, by exhibiting one of:

    seed-zsigmondy        Zsigmondy/Higgs exponent filter rejects a itself
    seed-non-higgs        some proper divisor m | a, a/m odd has 2^m+1 with
                          a non-Higgs prime (rigorous)
    v2-overshoot          deterministic lower-bound closure with all
                          factored proper-divisor seed contributions
                          exceeds budget a+1 (withdrawn in version 3)

If neither holds, a is reported as "unresolved" and counted.

This is a finitely checkable certificate. Version 2 argued as follows:
    1. Non-Higgs prime in 2^m+1 forces a non-Higgs prime in 2^a+1, which
       contradicts that every prime factor of a UPN is 3-Higgs. (Valid only
       under Hypothesis H3.)
    2. Lower-bound closure is monotone: more seed primes only increase
       the v2 contribution, so overshoot is preserved when more factors
       are added. (False: v_2(p^e+1) can decrease as e grows.)
"""

from __future__ import annotations

from collections import Counter, defaultdict
import argparse
import json
import os
import sys
from math import isqrt
from typing import Dict, List, Tuple

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CACHE_FILE = os.environ.get("FACTOR_CACHE", os.path.join(SCRIPT_DIR, "..", "factor_cache.json"))

from unitary_perfect_research import (
    SearchConfig,
    factor,
    format_factor,
    is_higgs_prime,
    is_probable_prime,
    q_plus_factor_tuple,
    zsigmondy_exponent_allowed,
)

try:
    from factordb_seed_probe import factordb_query, factors_from_response
except Exception:  # pragma: no cover
    factordb_query = None
    factors_from_response = None


IMPOSTORS: Tuple[Tuple[str, Dict[int, int], int, int], ...] = (
    ("3^2 5^3", {3: 2, 5: 3}, 10, 20),
    ("3^4 41", {3: 4, 41: 1}, 9, 18),
    ("5^2 13^2", {5: 2, 13: 2}, 6, 12),
    ("5^4 157^2 313", {5: 4, 157: 2, 313: 1}, 130, 260),
    ("5^4 29 157^2 313", {5: 4, 29: 1, 157: 2, 313: 1}, 26, 52),
)


def proper_seed_divisors(a: int) -> List[int]:
    """m such that m | a, a/m odd, m < a."""
    out = set()
    for d in range(1, isqrt(a) + 1):
        if a % d == 0:
            q = a // d
            if q % 2 == 1 and d < a:
                out.add(d)
            if d % 2 == 1 and q < a:
                out.add(q)
    return sorted(out)


_CACHE_LOAD_FAILED = False


def _load_cache_from_disk() -> Dict[int, Dict[int, int]]:
    global _CACHE_LOAD_FAILED
    if hasattr(sys, "set_int_max_str_digits"):
        sys.set_int_max_str_digits(0)
    if not os.path.exists(CACHE_FILE):
        return {}
    try:
        with open(CACHE_FILE) as f:
            raw = json.load(f)
        return {int(m): {int(p): int(e) for p, e in d.items()} for m, d in raw.items()}
    except Exception as e:
        # SAFETY: if load fails, set a flag so save refuses to write.
        # Otherwise a transient JSON error would let us overwrite the cache
        # with a tiny replacement.
        _CACHE_LOAD_FAILED = True
        print(f"WARNING: cache load failed ({e}); save disabled to prevent overwrite",
              file=sys.stderr)
        return {}


def _save_cache_to_disk(cache: Dict[int, Dict[int, int]]) -> None:
    if _CACHE_LOAD_FAILED:
        print("REFUSING to save: prior load failed; manual recovery required",
              file=sys.stderr)
        return
    try:
        # Write to a temp file and rename atomically.
        tmp = CACHE_FILE + ".tmp"
        with open(tmp, "w") as f:
            json.dump({str(m): {str(p): e for p, e in d.items()} for m, d in cache.items()}, f)
        os.replace(tmp, CACHE_FILE)
    except Exception as e:
        print(f"WARNING: cache save failed: {e}", file=sys.stderr)


def precompute_factorizations(
    needed_ms: List[int],
    per_call_timeout: float = 4.0,
    use_factordb: bool = False,
    factordb_timeout: float = 30.0,
    verbose: bool = True,
) -> Tuple[Dict[int, Dict[int, int]], List[int]]:
    """Factor 2^m+1 for m in needed_ms. Returns (cache, still_skipped).

    First try local factor with per_call_timeout. If that times out and
    use_factordb is True, fall back to FactorDB. Anything still failing is
    appended to still_skipped."""
    import signal
    import sys

    class _T(Exception):
        pass

    def _h(signum, frame):
        raise _T()

    signal.signal(signal.SIGALRM, _h)
    cache: Dict[int, Dict[int, int]] = _load_cache_from_disk()
    if verbose and cache:
        print(f"  loaded {len(cache)} cached factorizations from disk", flush=True)
    skipped_local: List[int] = []
    ms = sorted(set(needed_ms))
    total = len(ms)
    # Only attempt local factoring on m up to a soft cap to avoid
    # Pollard-Rho being unresponsive to SIGALRM on huge n.
    local_cap = 130
    for i, m in enumerate(ms):
        if m in cache:
            continue
        if m > local_cap:
            skipped_local.append(m)
            continue
        signal.setitimer(signal.ITIMER_REAL, per_call_timeout)
        try:
            cache[m] = factor(2**m + 1)
        except _T:
            skipped_local.append(m)
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
        if verbose and (i + 1) % 25 == 0:
            print(f"  ... local factored {i+1}/{total} (cached: {len(cache)}, skipped: {len(skipped_local)})", flush=True)
    if verbose:
        print(f"  local pass done: cached {len(cache)}, {len(skipped_local)} skipped or timed out", flush=True)
    if not use_factordb or factordb_query is None:
        _save_cache_to_disk(cache)
        return cache, skipped_local
    still_skipped: List[int] = []
    skipped_to_check = [m for m in skipped_local if m not in cache]
    for i, m in enumerate(skipped_to_check):
        try:
            data = factordb_query(2**m + 1, factordb_timeout)
            factors_map = factors_from_response(data)
            if factors_map is not None:
                cache[m] = factors_map
            else:
                # Save partial factor data; we will scan the known prime
                # factors for non-Higgs witnesses even when the full
                # factorization is unavailable.
                partial: Dict[int, int] = {}
                for p_text, exp in data.get("factors", []):
                    p = int(p_text)
                    e = int(exp)
                    if is_probable_prime(p):
                        partial[p] = partial.get(p, 0) + e
                if partial:
                    cache[m] = partial
                else:
                    still_skipped.append(m)
        except Exception:
            still_skipped.append(m)
        if verbose and (i + 1) % 10 == 0:
            print(f"  ... factordb {i+1}/{len(skipped_to_check)} (last m={m})", flush=True)
        if (i + 1) % 25 == 0:
            _save_cache_to_disk(cache)
    _save_cache_to_disk(cache)
    return cache, still_skipped


def collect_needed_ms(max_a: int, factor_limit: int) -> List[int]:
    """Collect all m values needed for the certificate.

    For each impostor and each a in its class up to max_a, include all
    proper divisors m | a with a/m odd that satisfy m <= factor_limit,
    plus a itself (regardless of factor_limit so the seed-non-Higgs
    check on m = a always has data when FactorDB can provide it)."""
    needed = set()
    for label, forced, residue, modulus in IMPOSTORS:
        a = residue
        while a <= max_a:
            for m in proper_seed_divisors(a):
                if m <= factor_limit:
                    needed.add(m)
            # Always include a itself; precompute will route to FactorDB
            # if a > local_cap. This makes the m=a non-Higgs check available.
            needed.add(a)
            a += modulus
    return sorted(needed)


_higgs_memo: Dict[int, bool] = {}

# Hardcoded list of large primes verified non-3-Higgs by extended check.
# These shortcut the per-cert Higgs verification when they recur.
_KNOWN_NON_HIGGS = {
    234166918021181298229848923411137491270809,  # 42-digit, factor of 2^1051+1, v_3(p-1) issue
    12366664051937762461361289279871201655390521,  # 44-digit, factor of 2^3099+1
}
for _p in _KNOWN_NON_HIGGS:
    _higgs_memo[_p] = False


def _fast_higgs(p: int, time_limit: float = 1.0) -> bool | None:
    """Returns is_higgs_prime(p) if computable within time_limit, else None."""
    if p in _higgs_memo:
        return _higgs_memo[p]
    import signal

    class _T(Exception):
        pass

    def _h(signum, frame):
        raise _T()

    old = signal.signal(signal.SIGALRM, _h)
    signal.setitimer(signal.ITIMER_REAL, time_limit)
    try:
        result = is_higgs_prime(p)
        _higgs_memo[p] = result
        return result
    except _T:
        return None
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old)


def non_higgs_in(factors: Dict[int, int]) -> int | None:
    # First pass: cheap primes (<=10^6) for fast lookup
    for p in factors:
        if p > 10**6:
            continue
        if not is_probable_prime(p):
            continue
        if not is_higgs_prime(p):
            return p
    # Second pass: medium primes (<=10^15) with short timeout
    for p in factors:
        if p <= 10**6 or p > 10**15:
            continue
        if not is_probable_prime(p):
            continue
        result = _fast_higgs(p, time_limit=2.0)
        if result is False:
            return p
    # Third pass: large primes with capped time budget.
    # Cache via _higgs_memo means repeated checks are free; we spend more
    # time per prime here because the result is informative if found.
    for p in factors:
        if p <= 10**15:
            continue
        if not is_probable_prime(p):
            continue
        digits = len(str(p))
        # Generous up to 60s for 100-digit primes; 90s for larger.
        # Primes that still time out remain "unknown".
        time_limit = min(90.0, max(15.0, 1.0 * digits))
        result = _fast_higgs(p, time_limit=time_limit)
        if result is False:
            return p
    return None


def seed_non_higgs_witness(
    a: int,
    factor_cache: Dict[int, Dict[int, int]],
) -> Tuple[int, int] | None:
    """Return (m, bad_prime) if some factored divisor m | a, a/m odd (including a itself) exposes
    a non-Higgs prime. Otherwise None."""
    for m in proper_seed_divisors(a) + [a]:
        if m not in factor_cache:
            continue
        bad = non_higgs_in(factor_cache[m])
        if bad is not None:
            return (m, bad)
    return None


def seed_lower_bound(
    a: int,
    factor_cache: Dict[int, Dict[int, int]],
) -> Dict[int, int]:
    """Lower bound on factors of 2^a+1.

    If a itself is factored, use it directly (the strongest lower bound,
    in fact it's the full factorization). Otherwise take the max
    valuation across all factored proper divisors m | a with a/m odd."""
    seed_lower: Dict[int, int] = {}
    if a in factor_cache:
        return dict(factor_cache[a])
    for m in proper_seed_divisors(a):
        if m not in factor_cache:
            continue
        for p, e in factor_cache[m].items():
            seed_lower[p] = max(seed_lower.get(p, 0), e)
    return seed_lower


def cascade_overshoot(
    a: int,
    forced: Dict[int, int],
    seed_lower: Dict[int, int],
    max_rounds: int,
    max_bases: int,
    max_exp: int,
    component_cap: int = 10**30,
) -> Tuple[str, Dict[str, int | str]]:
    targets = dict(seed_lower)
    for p, e in forced.items():
        targets[p] = max(targets.get(p, 0), e)
    budget = a + 1
    last_v2 = 0
    last_bases = len(targets)

    for step in range(max_rounds):
        if len(targets) > max_bases:
            return "growth-bases", {"step": step, "bases": len(targets), "last_v2": last_v2}
        if max(targets.values()) > max_exp:
            return "growth-exp", {"step": step, "max_exp": max(targets.values()), "last_v2": last_v2, "last_bases": last_bases}
        needs = defaultdict(int, seed_lower)
        v2 = 0
        for p, e in sorted(targets.items()):
            if p**e > component_cap:
                continue
            for r, valuation in q_plus_factor_tuple(p, e):
                if r == 2:
                    v2 += valuation
                else:
                    needs[r] += valuation
        last_v2 = v2
        last_bases = len(targets)
        if v2 > budget:
            return "v2-overshoot", {"step": step, "v2": v2, "bases": len(targets)}
        changed = False
        for p, need in list(needs.items()):
            if targets.get(p, 0) < need:
                targets[p] = need
                changed = True
        if not changed:
            return "stable-under", {"step": step, "v2": v2, "deficit": budget - v2}
    return "rounds-limit", {"step": max_rounds, "bases": len(targets), "last_v2": last_v2}


def class_exponents(residue: int, modulus: int, max_a: int) -> List[int]:
    out = []
    a = residue
    while a <= max_a:
        out.append(a)
        a += modulus
    return out


def certify_kernel(
    label: str,
    forced: Dict[int, int],
    residue: int,
    modulus: int,
    max_a: int,
    factor_cache: Dict[int, Dict[int, int]],
    cfg: SearchConfig,
    max_rounds: int,
    max_bases: int,
    max_exp: int,
    show: int,
) -> Counter:
    counts: Counter = Counter()
    examples: Dict[str, List[str]] = defaultdict(list)
    import sys
    print(f"=== {label} (forced={format_factor(forced)}, a=={residue} mod {modulus}) ===", flush=True)
    for a in class_exponents(residue, modulus, max_a):
        if not zsigmondy_exponent_allowed(2, a, cfg):
            status, detail = "seed-zsigmondy", {}
        else:
            witness = seed_non_higgs_witness(a, factor_cache)
            if witness is not None:
                m, bad = witness
                status, detail = "seed-non-higgs", {"m": m, "bad": bad}
            else:
                seed_lower = seed_lower_bound(a, factor_cache)
                status, detail = cascade_overshoot(
                    a, forced, seed_lower, max_rounds, max_bases, max_exp
                )
        counts[status] += 1
        if len(examples[status]) < show:
            examples[status].append(f"a={a} {detail}")

    total = sum(counts.values())
    print(f"total candidates: {total}")
    for status, count in sorted(counts.items()):
        print(f"  {status}: {count}")
    for status in sorted(examples):
        print(f"  {status}_examples:")
        for line in examples[status]:
            print(f"    {line}")
    return counts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-a", type=int, default=2000)
    parser.add_argument("--factor-limit", type=int, default=240,
                        help="precompute factorization of 2^m+1 for m up to here")
    parser.add_argument("--factor-timeout", type=float, default=4.0,
                        help="per-factorization timeout in seconds (skip slow ones)")
    parser.add_argument("--use-factordb", action="store_true",
                        help="for m that time out locally, query FactorDB")
    parser.add_argument("--factordb-timeout", type=float, default=30.0)
    parser.add_argument("--max-rounds", type=int, default=20)
    parser.add_argument("--max-bases", type=int, default=600)
    parser.add_argument("--max-exp", type=int, default=200)
    parser.add_argument("--show", type=int, default=8)
    args = parser.parse_args()

    needed_ms = collect_needed_ms(args.max_a, args.factor_limit)
    print(f"Need to factor 2^m+1 for {len(needed_ms)} m values "
          f"(max m = {needed_ms[-1] if needed_ms else 0}) ...")
    factor_cache, skipped = precompute_factorizations(
        needed_ms, args.factor_timeout,
        use_factordb=args.use_factordb, factordb_timeout=args.factordb_timeout,
    )
    # Fast diagnostic: count m with small-prime non-Higgs factors (do NOT
    # invoke the expensive large-prime Higgs check here; it would multiply
    # the diagnostic cost by the cache size).
    def _fast_diag_non_higgs(factors: Dict[int, int]) -> bool:
        for p in factors:
            if p > 10**8:
                continue
            if not is_probable_prime(p):
                continue
            if not is_higgs_prime(p):
                return True
        return False
    bad_ms = [m for m in factor_cache if _fast_diag_non_higgs(factor_cache[m])]
    print(f"  cached {len(factor_cache)} factorizations, "
          f"found {len(bad_ms)} m with non-Higgs factor < 10^8 in 2^m+1")
    if skipped:
        print(f"  still skipped (timed out): {skipped}")

    cfg = SearchConfig(
        min_even_exp=1,
        max_even_exp=args.max_a,
        max_odd_exp=64,
        max_bases=48,
        max_prime=10**12,
        max_component=10**18,
        max_solutions=1,
        higgs_filter=True,
        zsigmondy_exponent_filter=True,
    )

    aggregate: Counter = Counter()
    unresolved_by_kernel: Dict[str, List[int]] = {}
    for label, forced, residue, modulus in IMPOSTORS:
        counts = certify_kernel(
            label, forced, residue, modulus, args.max_a,
            factor_cache, cfg,
            args.max_rounds, args.max_bases, args.max_exp,
            args.show,
        )
        aggregate.update(counts)
        unresolved = ["stable-under", "rounds-limit", "growth-bases", "growth-exp"]
        unresolved_by_kernel[label] = [s for s in unresolved if counts.get(s, 0) > 0]

    print()
    print("=== Aggregate summary ===")
    total = sum(aggregate.values())
    print(f"total candidates across all impostors: {total}")
    for status, count in sorted(aggregate.items()):
        print(f"  {status}: {count}")
    print()
    print("Kernels with unresolved cases:")
    for label, statuses in unresolved_by_kernel.items():
        if statuses:
            print(f"  {label}: {statuses}")


if __name__ == "__main__":
    main()
