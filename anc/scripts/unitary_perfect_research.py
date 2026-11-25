#!/usr/bin/env python3
"""
Small exact-search harness for unitary perfect numbers.

NOTE (version 3 of the paper): the --higgs-filter option assumes Hypothesis H3,
that every prime factor of a unitary perfect number is 3-Higgs, which is
unproven; the 2-adic overshoot step is withdrawn. The exact search of the paper
is in ../search/.

The equation is

    prod_i (p_i ** e_i + 1) = 2 * prod_i p_i ** e_i.

For each selected component p^e, factor p^e + 1. The odd prime
valuations accumulated from these factorizations must exactly equal the
selected exponents e. For the prime 2, the accumulated 2-adic valuation
from the odd components must be e_2 + 1.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import gcd, isqrt
import argparse
import math
import random
import time
from typing import Dict, Iterable, List, Tuple


def is_probable_prime(n: int) -> bool:
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for p in small:
        if n == p:
            return True
        if n % p == 0:
            return False
    d = n - 1
    s = 0
    while d % 2 == 0:
        s += 1
        d //= 2

    # Deterministic for n < 2^64; still a strong probable-prime test above.
    bases = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for a in bases:
        if a >= n:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                break
        else:
            return False
    return True


def pollard_rho(n: int) -> int:
    if n % 2 == 0:
        return 2
    if n % 3 == 0:
        return 3
    while True:
        c = random.randrange(1, n - 1)
        x = random.randrange(2, n - 1)
        y = x
        d = 1
        while d == 1:
            x = (x * x + c) % n
            y = (y * y + c) % n
            y = (y * y + c) % n
            d = gcd(abs(x - y), n)
        if d != n:
            return d


@lru_cache(maxsize=None)
def small_prime_powers(bound: int) -> Tuple[int, ...]:
    powers = []
    for p in primes_upto(bound):
        q = p
        while q * p <= bound:
            q *= p
        powers.append(q)
    return tuple(powers)


def pollard_pm1(n: int, bound: int = 20000) -> int | None:
    if n % 2 == 0:
        return 2
    for base in (2, 3, 5):
        if base >= n:
            continue
        a = base % n
        for q in small_prime_powers(bound):
            a = pow(a, q, n)
        d = gcd(a - 1, n)
        if 1 < d < n:
            return d
    return None


def factor_into(n: int, out: Dict[int, int]) -> None:
    if n == 1:
        return
    if is_probable_prime(n):
        out[n] += 1
        return
    # Trial division is faster for the small factors that dominate q + 1.
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            while n % p == 0:
                out[p] += 1
                n //= p
            factor_into(n, out)
            return
    d_pm1 = pollard_pm1(n)
    if d_pm1 is not None:
        factor_into(d_pm1, out)
        factor_into(n // d_pm1, out)
        return
    d = pollard_rho(n)
    factor_into(d, out)
    factor_into(n // d, out)


def factor(n: int) -> Dict[int, int]:
    out: Dict[int, int] = defaultdict(int)
    factor_into(n, out)
    return dict(sorted(out.items()))


@lru_cache(maxsize=None)
def factor_tuple(n: int) -> Tuple[Tuple[int, int], ...]:
    return tuple(factor(n).items())


@lru_cache(maxsize=None)
def q_plus_factor_tuple(p: int, e: int) -> Tuple[Tuple[int, int], ...]:
    return factor_tuple(p**e + 1)


@lru_cache(maxsize=None)
def q_plus_v2(p: int, e: int) -> int:
    return dict(q_plus_factor_tuple(p, e)).get(2, 0)


def format_factor(factors: Dict[int, int]) -> str:
    parts = []
    for p, e in sorted(factors.items()):
        parts.append(str(p) if e == 1 else f"{p}^{e}")
    return " * ".join(parts) if parts else "1"


def primes_upto(n: int) -> List[int]:
    if n < 2:
        return []
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[:2] = b"\x00\x00"
    for i in range(2, isqrt(n) + 1):
        if sieve[i]:
            start = i * i
            sieve[start : n + 1 : i] = b"\x00" * (((n - start) // i) + 1)
    return [i for i in range(n + 1) if sieve[i]]


def partial_factor_by_primes(n: int, primes: Iterable[int]) -> Tuple[Dict[int, int], int]:
    factors: Dict[int, int] = {}
    for p in primes:
        if n % p == 0:
            count = 0
            while n % p == 0:
                count += 1
                n //= p
            factors[p] = count
    return factors, n


def sigma_star_from_factor(factors: Dict[int, int]) -> int:
    result = 1
    for p, e in factors.items():
        result *= p**e + 1
    return result


def verify(n: int) -> Tuple[bool, Dict[int, int], int]:
    factors = factor(n)
    sigma = sigma_star_from_factor(factors)
    return sigma == 2 * n, factors, sigma


@dataclass(frozen=True)
class SearchConfig:
    min_even_exp: int
    max_even_exp: int
    max_odd_exp: int
    max_bases: int
    max_prime: int
    max_component: int
    max_solutions: int
    stats: bool = False
    time_limit: float = 0.0
    reset_seen_per_even: bool = False
    min_product_prune: bool = False
    forced_components: Tuple[Tuple[int, int], ...] = ()
    forbidden_primes: Tuple[int, ...] = ()
    next_strategy: str = "min-candidates"
    higgs_filter: bool = False
    zsigmondy_exponent_filter: bool = False
    seed_divisor_filter: bool = False


@dataclass
class State:
    even_exp: int
    targets: Dict[int, int]
    needs: Dict[int, int]
    v2_from_odd: int
    product: Fraction
    path: List[Tuple[int, int]]

    def key(self) -> Tuple:
        target_items = tuple(sorted(self.targets.items()))
        need_items = tuple(sorted((p, e) for p, e in self.needs.items() if e))
        return target_items, need_items, self.v2_from_odd, self.product


def add_component(state: State, p: int, e: int) -> State | None:
    q = p**e
    q_factors = dict(q_plus_factor_tuple(p, e))
    targets = dict(state.targets)
    needs = dict(state.needs)
    v2 = state.v2_from_odd

    targets[p] = e
    for r, valuation in q_factors.items():
        if r == 2:
            if p != 2:
                v2 += valuation
        else:
            needs[r] = needs.get(r, 0) + valuation

    product = state.product * Fraction(q + 1, q)
    path = state.path + [(p, e)]
    return State(state.even_exp, targets, needs, v2, product, path)


def viable(state: State, cfg: SearchConfig) -> bool:
    forbidden = set(cfg.forbidden_primes)
    if forbidden & set(state.targets):
        return False
    if forbidden & set(state.needs):
        return False
    if cfg.higgs_filter:
        if any(not is_higgs_prime(p) for p in state.targets):
            return False
        if any(not is_higgs_prime(p) for p in state.needs if p != 2):
            return False
    if len(state.targets) > cfg.max_bases:
        return False
    if state.v2_from_odd > state.even_exp + 1:
        return False
    if state.product > 2:
        return False
    unassigned = [
        (p, need)
        for p, need in state.needs.items()
        if p != 2 and p not in state.targets
    ]
    min_extra_v2 = 0
    min_forced_log = 0.0
    for p, need in unassigned:
        feasible = feasible_exponents(p, need, cfg)
        if not feasible:
            return False
        min_extra_v2 += min(q_plus_v2(p, e) for e in feasible)
        if cfg.min_product_prune:
            max_e = max(feasible)
            min_forced_log += math.log1p(1 / (p**max_e))
    if state.v2_from_odd + min_extra_v2 > state.even_exp + 1:
        return False
    if cfg.min_product_prune and (
        math.log(state.product.numerator) - math.log(state.product.denominator) + min_forced_log
        > math.log(2)
    ):
        return False
    for p, need in state.needs.items():
        if p > cfg.max_prime:
            return False
        assigned = state.targets.get(p)
        if assigned is not None and need > assigned:
            return False
        min_exp = assigned if assigned is not None else need
        if p**min_exp > cfg.max_component:
            return False
    return True


def feasible_exponents(p: int, need: int, cfg: SearchConfig) -> List[int]:
    return [
        e
        for e in range(max(need, 1), cfg.max_odd_exp + 1)
        if p**e <= cfg.max_component
        and zsigmondy_exponent_allowed(p, e, cfg)
    ]


def zsigmondy_exponent_allowed(p: int, e: int, cfg: SearchConfig) -> bool:
    if not cfg.zsigmondy_exponent_filter:
        return True
    if e <= 1:
        return True
    # Zsigmondy's only relevant exception for p^e + 1 is 2^3 + 1.
    if p == 2 and e == 3:
        return True
    for r, valuation in factor_tuple(2 * e):
        if valuation > 3:
            return False
        if cfg.higgs_filter and not is_higgs_prime(r):
            return False
    return True


@lru_cache(maxsize=None)
def odd_quotient_divisors(a: int) -> Tuple[int, ...]:
    divisors = []
    for d in range(1, isqrt(a) + 1):
        if a % d != 0:
            continue
        q = a // d
        if d < a and q % 2 == 1:
            divisors.append(d)
        if q != d and q < a and d % 2 == 1:
            divisors.append(q)
    return tuple(sorted(set(divisors)))


def seed_divisors_allowed(a: int, cfg: SearchConfig) -> bool:
    for d in odd_quotient_divisors(a):
        factors = factor(2**d + 1)
        if any(p > cfg.max_prime for p in factors):
            return False
        if any(e > cfg.max_odd_exp for e in factors.values()):
            return False
        if any(p**e > cfg.max_component for p, e in factors.items()):
            return False
        if cfg.higgs_filter and any(not is_higgs_prime(p) for p in factors):
            return False
    return True


def locally_feasible_exponents(state: State, p: int, need: int, cfg: SearchConfig) -> List[int]:
    return [
        e
        for e in feasible_exponents(p, need, cfg)
        if locally_feasible_component(state, p, e, cfg)
    ]


def locally_feasible_component(state: State, p: int, e: int, cfg: SearchConfig) -> bool:
    if p in cfg.forbidden_primes:
        return False
    if cfg.higgs_filter and not is_higgs_prime(p):
        return False
    if len(state.targets) + (0 if p in state.targets else 1) > cfg.max_bases:
        return False
    q = p**e
    if q > cfg.max_component:
        return False
    if state.product * Fraction(q + 1, q) > 2:
        return False

    next_v2 = state.v2_from_odd
    next_needs = dict(state.needs)
    for r, valuation in q_plus_factor_tuple(p, e):
        if r == 2:
            if p != 2:
                next_v2 += valuation
                if next_v2 > state.even_exp + 1:
                    return False
            continue
        if r > cfg.max_prime:
            return False
        if r in cfg.forbidden_primes:
            return False
        if cfg.higgs_filter and not is_higgs_prime(r):
            return False
        next_need = next_needs.get(r, 0) + valuation
        next_needs[r] = next_need
        assigned = state.targets.get(r)
        if assigned is not None and next_need > assigned:
            return False
        min_exp = assigned if assigned is not None else next_need
        if min_exp > cfg.max_odd_exp:
            return False
        if r**min_exp > cfg.max_component:
            return False
    return True


@lru_cache(maxsize=None)
def is_higgs_prime(p: int) -> bool:
    if p == 2:
        return True
    if p < 2 or not is_probable_prime(p):
        return False
    for r, exponent in factor_tuple(p - 1):
        if exponent > 3 or not is_higgs_prime(r):
            return False
    return True


@lru_cache(maxsize=None)
def higgs_primes_upto(limit: int) -> frozenset[int]:
    """Return 3-Higgs primes up to limit.

    A prime p is admitted if p - 1 divides the cube of the product of
    smaller admitted primes. This is OEIS A057447. Hypothesis H3 of the
    paper, which is unproven, says that every prime factor of a unitary
    perfect number is of this kind.
    """
    return frozenset(p for p in primes_upto(limit) if is_higgs_prime(p))


def is_complete(state: State) -> bool:
    if state.v2_from_odd != state.even_exp + 1:
        return False
    if state.product != 2:
        return False
    for p, e in state.targets.items():
        if p == 2:
            continue
        if state.needs.get(p, 0) != e:
            return False
    for p, need in state.needs.items():
        if p != 2 and state.targets.get(p) != need:
            return False
    return True


def next_unresolved_prime(state: State, cfg: SearchConfig) -> int | None:
    unassigned = []
    for p, need in state.needs.items():
        if p == 2:
            continue
        target = state.targets.get(p)
        if target is None:
            feasible = locally_feasible_exponents(state, p, need, cfg)
            unassigned.append((p, need, feasible))
    if not unassigned:
        return None
    if cfg.next_strategy == "smallest":
        return min(unassigned, key=lambda item: item[0])[0]
    # Choose the most constrained prime. Ties go to larger required valuation,
    # then smaller prime for reproducibility.
    return min(unassigned, key=lambda item: (len(item[2]), -item[1], item[0]))[0]


def search(cfg: SearchConfig) -> List[State]:
    solutions: List[State] = []
    seen = set()
    counts = defaultdict(int)
    start_time = time.monotonic()

    def dfs(state: State) -> None:
        counts["visited"] += 1
        if len(solutions) >= cfg.max_solutions:
            return
        if cfg.time_limit and time.monotonic() - start_time > cfg.time_limit:
            counts["timed_out"] += 1
            return
        if not viable(state, cfg):
            counts["not_viable"] += 1
            return
        key = state.key()
        if key in seen:
            counts["seen"] += 1
            return
        seen.add(key)
        if is_complete(state):
            solutions.append(state)
            counts["solutions"] += 1
            return

        p = next_unresolved_prime(state, cfg)
        if p is None:
            counts["closed_not_exact"] += 1
            return

        need = state.needs.get(p, 0)
        exponents = locally_feasible_exponents(state, p, need, cfg)
        # Low exponents maximize reciprocal product and usually close small
        # examples; try same-parity choices first because they often minimize
        # the 2-adic overshoot.
        candidates: Iterable[int] = sorted(
            exponents,
            key=lambda e: (q_plus_v2(p, e), e),
        )

        for e in candidates:
            if cfg.time_limit and time.monotonic() - start_time > cfg.time_limit:
                counts["timed_out"] += 1
                return
            nxt = add_component(state, p, e)
            if nxt is not None:
                dfs(nxt)

    for even_exp in range(cfg.min_even_exp, cfg.max_even_exp + 1):
        if not zsigmondy_exponent_allowed(2, even_exp, cfg):
            continue
        if cfg.seed_divisor_filter and not seed_divisors_allowed(even_exp, cfg):
            continue
        if cfg.reset_seen_per_even:
            seen.clear()
        base = State(
            even_exp=even_exp,
            targets={},
            needs={},
            v2_from_odd=0,
            product=Fraction(1, 1),
            path=[],
        )
        start = add_component(base, 2, even_exp)
        for p, e in cfg.forced_components:
            if start is None:
                break
            start = add_component(start, p, e)
        if start is not None:
            dfs(start)
    if cfg.stats:
        print("stats:")
        for key, value in sorted(counts.items()):
            print(f"  {key}: {value}")
        print(f"  seen_states: {len(seen)}")
        print(f"  elapsed_seconds: {time.monotonic() - start_time:.3f}")
    return solutions


def state_to_n(state: State) -> int:
    n = 1
    for p, e in state.targets.items():
        n *= p**e
    return n


def command_verify(args: argparse.Namespace) -> None:
    for n in args.n:
        ok, factors, sigma = verify(n)
        print(f"{n}: {'OK' if ok else 'NO'}")
        print(f"  factor: {format_factor(factors)}")
        print(f"  sigma*: {sigma}")


def command_search(args: argparse.Namespace) -> None:
    forced_components = tuple(parse_component_arg(item) for item in args.force_component)
    forbidden_primes = tuple(args.forbid_prime)
    cfg = SearchConfig(
        min_even_exp=args.min_even_exp,
        max_even_exp=args.max_even_exp,
        max_odd_exp=args.max_odd_exp,
        max_bases=args.max_bases,
        max_prime=args.max_prime,
        max_component=args.max_component,
        max_solutions=args.max_solutions,
        stats=args.stats,
        time_limit=args.time_limit,
        reset_seen_per_even=args.reset_seen_per_even,
        min_product_prune=args.min_product_prune,
        forced_components=forced_components,
        forbidden_primes=forbidden_primes,
        next_strategy=args.next_strategy,
        higgs_filter=args.higgs_filter,
        zsigmondy_exponent_filter=args.zsigmondy_exponent_filter,
        seed_divisor_filter=args.seed_divisor_filter,
    )
    solutions = search(cfg)
    for state in sorted(solutions, key=state_to_n):
        n = state_to_n(state)
        print(n)
        print(f"  factor: {format_factor(state.targets)}")
        print(f"  path: {state.path}")
    print(f"solutions: {len(solutions)}")


def parse_component_arg(text: str) -> Tuple[int, int]:
    if ":" in text:
        p, e = text.split(":", 1)
        return int(p), int(e)
    return int(text), 1


def command_seed_scan(args: argparse.Namespace) -> None:
    primes = primes_upto(args.max_prime) if not args.full_factor else []
    zsig_cfg = SearchConfig(
        min_even_exp=args.min_even_exp,
        max_even_exp=args.max_even_exp,
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
    for even_exp in range(args.min_even_exp, args.max_even_exp + 1):
        if not zsigmondy_exponent_allowed(2, even_exp, zsig_cfg):
            continue
        if args.seed_divisor_filter and not seed_divisors_allowed(even_exp, zsig_cfg):
            continue
        if args.full_factor:
            factors = factor(2**even_exp + 1)
            residual = 1
        else:
            factors, residual = partial_factor_by_primes(2**even_exp + 1, primes)
        if residual != 1:
            continue
        if args.higgs_filter and any(not is_higgs_prime(p) for p in factors):
            continue
        if any(p > args.max_prime for p in factors):
            continue
        if any(e > args.max_odd_exp for e in factors.values()):
            continue
        if any(p**e > args.max_component for p, e in factors.items()):
            continue
        survivors.append((even_exp, factors))

    print(f"survivors: {len(survivors)}")
    for even_exp, factors in survivors:
        print(f"{even_exp}: {format_factor(factors)}")


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(required=True)

    p_verify = sub.add_parser("verify")
    p_verify.add_argument("n", nargs="+", type=int)
    p_verify.set_defaults(func=command_verify)

    p_search = sub.add_parser("search")
    p_search.add_argument("--min-even-exp", type=int, default=1)
    p_search.add_argument("--max-even-exp", type=int, default=18)
    p_search.add_argument("--max-odd-exp", type=int, default=4)
    p_search.add_argument("--max-bases", type=int, default=12)
    p_search.add_argument("--max-prime", type=int, default=10000)
    p_search.add_argument("--max-component", type=int, default=300000)
    p_search.add_argument("--max-solutions", type=int, default=20)
    p_search.add_argument("--stats", action="store_true")
    p_search.add_argument("--time-limit", type=float, default=0.0)
    p_search.add_argument("--reset-seen-per-even", action="store_true")
    p_search.add_argument("--min-product-prune", action="store_true")
    p_search.add_argument("--force-component", action="append", default=[], help="Require odd component p:e")
    p_search.add_argument("--forbid-prime", action="append", type=int, default=[], help="Forbid a prime base")
    p_search.add_argument("--next-strategy", choices=("min-candidates", "smallest"), default="min-candidates")
    p_search.add_argument("--higgs-filter", action="store_true", help="Require all prime bases to be 3-Higgs primes")
    p_search.add_argument("--zsigmondy-exponent-filter", action="store_true", help="Use primitive-divisor Higgs constraints on exponents")
    p_search.add_argument("--seed-divisor-filter", action="store_true", help="Reject a if a smaller forced divisor 2^d+1 already violates filters")
    p_search.set_defaults(func=command_search)

    p_seed = sub.add_parser("seed-scan")
    p_seed.add_argument("--min-even-exp", type=int, default=1)
    p_seed.add_argument("--max-even-exp", type=int, default=1000)
    p_seed.add_argument("--max-odd-exp", type=int, default=4)
    p_seed.add_argument("--max-prime", type=int, default=1000000)
    p_seed.add_argument("--max-component", type=int, default=1000000)
    p_seed.add_argument("--higgs-filter", action="store_true")
    p_seed.add_argument("--zsigmondy-exponent-filter", action="store_true")
    p_seed.add_argument("--seed-divisor-filter", action="store_true")
    p_seed.add_argument("--full-factor", action="store_true", help="Fully factor 2^a+1 instead of trial division to max-prime")
    p_seed.set_defaults(func=command_seed_scan)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
