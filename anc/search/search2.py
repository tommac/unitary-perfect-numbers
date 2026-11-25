#!/usr/bin/env python3
"""Exhaustive search for unitary perfect numbers n = 2^a m, for one fixed a.

Exact facts used (and nothing else):
  (1) (2^a + 1) * prod (c + 1) = 2^(a+1) * prod c, over the odd components c = p^e
      of m (distinct odd primes p).  Equivalently prod (1 + 1/c) = T := 2^(a+1)/(2^a+1),
      and every odd prime has the same valuation on both sides.
  (2) Hence sum v_2(c + 1) = a + 1, and each odd component costs v_2(c + 1) >= 1.

Search: components in increasing order. State after choosing c_1 < ... < c_j:
  F = prod (1 + 1/c_i), S = sum v_2(c_i + 1), R = a + 1 - S (remaining budget),
  N[p] = valuation of odd p in (2^a+1) prod (c_i + 1), D[p] = chosen exponent of p,
  Q = primes with N[p] > 0 not yet chosen ("demanded").
Every future component costs >= 1, every demanded prime needs its own future
component, so |Q| <= R, and a component of a new prime needs a free slot.
The next component c satisfies F (1 + 1/c) <= T and F (1 + 1/c)^R >= T.

A branch is closed only by one of these exact arguments.  If the range of the next
component exceeds the sieve limit while a free slot remains, the branch is
reported as OPEN and the value of a is not claimed.
"""
import math
import sys
import time
from array import array
from bisect import bisect_left, bisect_right
from fractions import Fraction
from functools import lru_cache

from sympy import factorint

import os
LIMIT = int(os.environ.get("UPN_LIMIT", "20000000"))
FORCE = [tuple(int(x) for x in t.split(":")) for t in os.environ.get("UPN_FORCE", "").split(",") if t]
TIGHT = os.environ.get("UPN_TIGHT", "") not in ("", "0")


def build_spf(n):
    import numpy as np
    spf = np.zeros(n + 2, dtype=np.int32)
    for i in range(2, int((n + 1) ** 0.5) + 1):
        if spf[i] == 0:
            block = spf[i * i::i]
            block[block == 0] = i
            spf[i * i::i] = block
    idx = np.nonzero(spf == 0)[0]
    spf[idx] = idx
    return spf


def build_primes(n):
    bs = bytearray([1]) * (n + 1)
    bs[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if bs[i]:
            bs[i * i::i] = bytearray(len(range(i * i, n + 1, i)))
    return array("i", (i for i in range(3, n + 1, 2) if bs[i]))


PRIMES = build_primes(LIMIT)
SPF_LIMIT = min(LIMIT + 2, 20_000_002)
SPF = None


@lru_cache(maxsize=1_000_000)
def fac(n):
    """Factorization of n as a tuple of (prime, exponent), smallest prime first."""
    if SPF is not None and n < SPF_LIMIT:
        out = []
        while n > 1:
            p = int(SPF[n])
            k = 0
            while n % p == 0:
                n //= p
                k += 1
            out.append((p, k))
        return tuple(out)
    return tuple(sorted(factorint(n).items()))


def prime_power(c):
    """(s, g) if c = s^g with s prime, else None.  A 'composite' verdict from the
    probable-prime test is always correct, so no prime power is ever missed."""
    from sympy import integer_nthroot, isprime
    if c < 3:
        return None
    for g in range(c.bit_length(), 0, -1):
        if g == 1:
            return (c, 1) if isprime(c) else None
        r, exact = integer_nthroot(c, g)
        if exact and r > 2 and isprime(r):
            return (r, g)
    return None


def search(a):
    T = Fraction(2 ** (a + 1), 2 ** a + 1)
    B = a + 1
    N = {}
    for q, v in fac(2 ** a + 1):
        N[q] = v
    D = {}
    sols, opens = [], []
    nodes = [0]
    F0, S0, chosen0 = Fraction(1), 0, []
    for p, e in FORCE:                     # forced exact components (kernel), placed first
        c = p ** e
        D[p] = e
        for q, v in fac(c + 1):
            if q == 2:
                S0 += v
            else:
                N[q] = N.get(q, 0) + v
        F0 *= Fraction(c + 1, c)
        chosen0.append((p, e))
    if FORCE and (S0 > B or F0 > T or any(N.get(q, 0) > D[q] for q in D)):
        return [], [], 0, 0.0

    def emin(q, last, free=1):
        e = max(1, N[q])
        c = q ** e
        while c <= last:
            e += 1
            c *= q
        if TIGHT and free == 0 and q % 4 == 3 and e % 2 == 1:
            e += 1
            c *= q
        return e, c

    def rec(last, F, S, chosen):
        nodes[0] += 1
        R = B - S
        if R <= 0:
            return                                # F < T but no budget left for another component
        Q = [p for p, v in N.items() if v > 0 and p not in D]
        free = R - len(Q)
        if free < 0:
            return
        t_minus_f = T - F
        # Exact ceil(F / (T - F)) = ceil(1 / delta) with integer arithmetic only:
        num_inv = F.numerator * t_minus_f.denominator
        den_inv = F.denominator * t_minus_f.numerator
        lo = max(last + 1, (num_inv + den_inv - 1) // den_inv)
        mq = [emin(q, last, free)[1] for q in Q]
        if TIGHT and free == 0 and mq:
            lo = max(lo, min(mq))
        exp_free = free - 1 if TIGHT else free
        lhs_base = F.numerator * T.denominator
        rhs_base = T.numerator * F.denominator

        def geq_G(x):
            """Exact rational check of F * (1 + 1/x)^R >= T and F * G(x) >= T."""
            xp1 = x + 1
            if not TIGHT and lhs_base * (xp1 ** R) < rhs_base * (x ** R):
                return False
            lhs = lhs_base * xp1
            rhs = rhs_base * x
            if exp_free > 0:
                lhs *= (x + 2) ** exp_free
                rhs *= xp1 ** exp_free
            elif exp_free < 0:                    # TIGHT with free == 0
                lhs *= xp1
                rhs *= x + 2
            for m in mq:
                u = xp1 if xp1 >= m else m
                lhs *= u + 1
                rhs *= u
            return lhs >= rhs

        if not geq_G(lo):
            return                                # 100% exact prune: even x = lo cannot reach T

        # Floating-point is used ONLY for an initial bracket hint; exact geq_G certifies hi.
        guess = lo
        try:
            d = float(t_minus_f / F)
        except OverflowError:
            d = 0.0
        if d > 1e-14 and math.isfinite(d):
            u0 = int(1.0 / math.expm1(math.log1p(d) / R)) + 2
            if u0 > lo:
                target = math.log1p(d)

                def logG_hint(x):
                    return (math.log1p(1.0 / x)
                            + sum(math.log1p(1.0 / max(x + 1, m)) for m in mq)
                            + exp_free * math.log1p(1.0 / (x + 1)))

                a_h, b_h = lo, u0
                if logG_hint(a_h) >= target:
                    while b_h - a_h > 1:
                        mid = (a_h + b_h) // 2
                        if logG_hint(mid) >= target:
                            a_h = mid
                        else:
                            b_h = mid
                    guess = a_h

        if geq_G(guess):
            a_ = guess
            if not geq_G(a_ + 1):
                hi = a_
            else:
                step = 2
                while geq_G(a_ + step):
                    a_ += step
                    step *= 2
                b_ = a_ + step
                while b_ - a_ > 1:
                    mid = (a_ + b_) // 2
                    if geq_G(mid):
                        a_ = mid
                    else:
                        b_ = mid
                hi = a_
        else:
            b_ = guess
            if geq_G(b_ - 1):
                hi = b_ - 1
            else:
                step = 2
                while b_ - step > lo and not geq_G(b_ - step):
                    b_ -= step
                    step *= 2
                a_ = max(lo, b_ - step)
                while b_ - a_ > 1:
                    mid = (a_ + b_) // 2
                    if geq_G(mid):
                        a_ = mid
                    else:
                        b_ = mid
                hi = a_
        cands = []
        # powers of demanded primes
        for q in Q:
            e, c = emin(q, last)
            while c <= hi:
                if c >= lo:
                    cands.append((c, q, e))
                e += 1
                c *= q
        # components of new primes (need a free slot)
        if free == 1:
            # Exactly one slot is not reserved for a demanded prime.  A new component c
            # uses it, so c + 1 may bring no new prime and must cost exactly 1:
            # c + 1 = 2M with M odd and composed of primes already present.  Generate
            # those M directly; this is complete for every range, however large.
            caps = [(p, D[p] - N.get(p, 0)) for p in D if D[p] - N.get(p, 0) > 0]
            caps += [(q, None) for q in Q]
            mlo, mhi = (lo + 1 + 1) // 2, (hi + 1) // 2
            smooth = [1]
            for p, cap in caps:
                ext = []
                for m in smooth:
                    k, x = 0, m
                    while True:
                        k += 1
                        x *= p
                        if x > mhi or (cap is not None and k > cap):
                            break
                        ext.append(x)
                smooth += ext
            for m in smooth:
                if m < mlo:
                    continue
                c = 2 * m - 1
                pe = prime_power(c)
                if pe is None:
                    continue
                s, g = pe
                if s in D or s in N:
                    continue
                cands.append((c, s, g))
        elif free >= 2:
            if hi > LIMIT:
                opens.append((list(chosen), lo, hi, R, sorted(Q)))
            i, j = bisect_left(PRIMES, lo), bisect_right(PRIMES, min(hi, LIMIT))
            for k in range(i, j):
                s = PRIMES[k]
                if s not in D and s not in N:
                    cands.append((s, s, 1))
            for s in PRIMES:
                if s * s > hi:
                    break
                if s in D or s in N:
                    continue
                g, c = 2, s * s
                while c <= hi:
                    if c >= lo:
                        cands.append((c, s, g))
                    g += 1
                    c *= s
        cands.sort()
        for c, p, e in cands:
            f = fac(c + 1)
            s2 = f[0][1] if f[0][0] == 2 else 0
            S2 = S + s2
            if S2 > B:
                continue
            # cheap admissibility: saturated primes, and new demands must fit the budget
            newq = 0
            bad = False
            for q, v in f:
                if q == 2:
                    continue
                if q in D:
                    if N.get(q, 0) + v > D[q]:
                        bad = True
                        break
                elif N.get(q, 0) == 0:
                    newq += 1
            if bad:
                continue
            Qafter = len(Q) - (1 if p in N and p not in D and N[p] > 0 else 0) + newq
            if Qafter > B - S2:
                continue
            F2 = F * Fraction(c + 1, c)
            if F2 > T:
                continue
            D[p] = e
            for q, v in f:
                if q != 2:
                    N[q] = N.get(q, 0) + v
            if N.get(p, 0) <= e:
                if F2 == T:
                    if S2 == B and all(N.get(q, 0) == D.get(q, 0) for q in set(N) | set(D)):
                        sols.append(sorted(chosen + [(p, e)]))
                else:
                    rec(c, F2, S2, chosen + [(p, e)])
            for q, v in f:
                if q != 2:
                    N[q] -= v
                    if N[q] == 0:
                        del N[q]
            del D[p]

    t = time.time()
    if F0 == T:
        if S0 == B and all(N.get(q, 0) == D.get(q, 0) for q in set(N) | set(D)):
            sols.append(sorted(chosen0))
    else:
        rec(1, F0, S0, chosen0)
    return sols, opens, nodes[0], time.time() - t


if __name__ == "__main__":
    SPF = build_spf(SPF_LIMIT - 2)
    globals()["SPF"] = SPF
    for a in [int(x) for x in sys.argv[1].split(",")]:
        sols, opens, nodes, dt = search(a)
        print(f"a={a}: solutions={sols} OPEN={len(opens)} nodes={nodes} time={dt:.1f}s", flush=True)
        for o in opens[:3]:
            print("   open:", o, flush=True)
