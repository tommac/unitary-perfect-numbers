/*
 * upn.c -- exhaustive search for unitary perfect numbers n = 2^a m with a fixed.
 *
 * A port of search2.py (same algorithm, same pruning rules, same bounds), so that
 * both programs must visit exactly the same nodes.  Only exact facts are used:
 *   (2^a+1) prod (c+1) = 2^(a+1) prod c   over the odd components c = p^e of m.
 * Components are chosen in increasing order.  A branch is closed only by an exact
 * argument; anything the program cannot settle is printed as an OPEN line, and the
 * value of a is claimed only if no OPEN line is printed.
 *
 * Build:  cc -O2 -o upn upn.c -lgmp -lm  (with -I/-L for the GMP install)
 * Usage:  upn a "p1:e1,p2:e2,..." [worker nworkers splitdepth]
 *         (the second argument is the full factorization of 2^a+1)
 */
#include <gmp.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef unsigned long long u64;

#ifndef LIMIT
#define LIMIT 20000000ULL          /* primes up to LIMIT come from the table */
#endif
#define SPF_MAX (LIMIT + 2)        /* smallest-prime-factor table covers c+1 <= SPF_MAX */
#define BIGP (1ULL << 62)          /* primes are stored as u64 below this bound */

static uint32_t *SPF;
static uint32_t *PR;
static size_t NPR;

static int AEXP, BUD;
static mpq_t TT;

static long long NODES = 0, NOPEN = 0, NSOL = 0;
static int WORKER = 0, NWORK = 1, SPLITD = -1;
static long long SPLITCNT = 0;

#define MAXE 160
typedef struct { u64 p; int n; int d; } Ent;   /* n: numerator valuation; d: chosen exponent (0 = not chosen) */

#define MAXDEPTH 200
static u64 PATHP[MAXDEPTH];
static int PATHE[MAXDEPTH];

/* ------------------------------------------------------------------ sieve */
static void build_tables(void) {
    SPF = calloc(SPF_MAX + 1, sizeof(uint32_t));
    for (u64 i = 2; i <= SPF_MAX; i++) {
        if (SPF[i] == 0) {
            SPF[i] = (uint32_t)i;
            if (i * i <= SPF_MAX)
                for (u64 j = i * i; j <= SPF_MAX; j += i)
                    if (SPF[j] == 0) SPF[j] = (uint32_t)i;
        }
    }
    size_t cnt = 0;
    for (u64 i = 3; i <= LIMIT; i += 2) if (SPF[i] == i) cnt++;
    PR = malloc(cnt * sizeof(uint32_t));
    NPR = 0;
    for (u64 i = 3; i <= LIMIT; i += 2) if (SPF[i] == i) PR[NPR++] = (uint32_t)i;
}

/* ------------------------------------------------------------ factoring */
typedef struct { u64 p; int k; } PF;
typedef struct {
    int v2;            /* 2-adic valuation of c+1 */
    int n;             /* number of odd prime factors found */
    PF f[64];
    int unknown_lb;    /* 1 if c+1 has a prime factor above 2^62, which is then a new prime;
                          an unfactored composite part counts for nothing, since its
                          prime factors may already be in E */
} Fac;

static int is_prime_u64(u64 x) {
    if (x < 2) return 0;
    if (x <= SPF_MAX) return SPF[x] == x;
    mpz_t z; mpz_init(z); mpz_set_ui(z, 0);
    mpz_import(z, 1, 1, sizeof(u64), 0, 0, &x);
    int r = mpz_probab_prime_p(z, 30) > 0;   /* BPSW (+MR): deterministic below 2^64 */
    mpz_clear(z);
    return r;
}

static u64 mpz_to_u64(const mpz_t z) {
    u64 x = 0;
    size_t cnt = 0;
    mpz_export(&x, &cnt, -1, sizeof(u64), 0, 0, z);
    return cnt ? x : 0;
}

static void addf(Fac *F, u64 p, int k) {
    for (int i = 0; i < F->n; i++) if (F->f[i].p == p) { F->f[i].k += k; return; }
    F->f[F->n].p = p; F->f[F->n].k = k; F->n++;
}

static void spf_factor(u64 x, Fac *F) {
    while (x > 1) {
        u64 p = SPF[x]; int k = 0;
        while (x % p == 0) { x /= p; k++; }
        addf(F, p, k);
    }
}

/* Pollard-Brent rho on a composite z; returns a nontrivial factor in g (0 on failure). */
static int rho(const mpz_t z, mpz_t g) {
    mpz_t x, y, ys, q, t, c;
    mpz_inits(x, y, ys, q, t, c, NULL);
    int ok = 0;
    for (unsigned long cc = 1; cc < 20 && !ok; cc++) {
        mpz_set_ui(c, cc); mpz_set_ui(y, 2); mpz_set_ui(q, 1); mpz_set_ui(g, 1);
        unsigned long r = 1, m = 128;
        long iters = 0;
        do {
            mpz_set(x, y);
            for (unsigned long i = 0; i < r; i++) { mpz_mul(y, y, y); mpz_add(y, y, c); mpz_mod(y, y, z); }
            unsigned long k = 0;
            do {
                mpz_set(ys, y);
                unsigned long lim = (m < r - k) ? m : r - k;
                for (unsigned long i = 0; i < lim; i++) {
                    mpz_mul(y, y, y); mpz_add(y, y, c); mpz_mod(y, y, z);
                    mpz_sub(t, x, y); mpz_abs(t, t); mpz_mul(q, q, t); mpz_mod(q, q, z);
                }
                mpz_gcd(g, q, z);
                k += m;
            } while (k < r && mpz_cmp_ui(g, 1) == 0);
            r *= 2;
            iters += r;
        } while (mpz_cmp_ui(g, 1) == 0 && iters < 4000000);
        if (mpz_cmp(g, z) == 0) {
            do {
                mpz_mul(ys, ys, ys); mpz_add(ys, ys, c); mpz_mod(ys, ys, z);
                mpz_sub(t, x, ys); mpz_abs(t, t); mpz_gcd(g, t, z);
            } while (mpz_cmp_ui(g, 1) == 0);
        }
        if (mpz_cmp_ui(g, 1) > 0 && mpz_cmp(g, z) < 0) ok = 1;
    }
    mpz_clears(x, y, ys, q, t, c, NULL);
    return ok;
}

/* Factor an odd cofactor z completely into F (recursively).  Returns 0 if some part
 * could not be identified.  F->unknown_lb is set to 1 only when a prime above 2^62 is
 * found: every prime in E is below 2^62, so such a prime is certainly new. */
static int full_factor(mpz_t z, Fac *F, int depth) {
    if (mpz_cmp_ui(z, 1) == 0) return 1;
    if (mpz_sizeinbase(z, 2) <= 62 && mpz_to_u64(z) <= SPF_MAX) { spf_factor(mpz_to_u64(z), F); return 1; }
    if (mpz_probab_prime_p(z, 30) > 0) {
        if (mpz_sizeinbase(z, 2) <= 62) { addf(F, mpz_to_u64(z), 1); return 1; }
        F->unknown_lb = 1;            /* a prime too large to store: new, and reported as OPEN */
        return 0;
    }
    if (mpz_perfect_power_p(z)) {
        for (unsigned long k = mpz_sizeinbase(z, 2); k >= 2; k--) {
            mpz_t r; mpz_init(r);
            if (mpz_root(r, z, k)) {
                Fac G; memset(&G, 0, sizeof G);
                int ok = full_factor(r, &G, depth + 1);
                for (int i = 0; i < G.n; i++) addf(F, G.f[i].p, G.f[i].k * (int)k);
                if (G.unknown_lb) F->unknown_lb = 1;
                mpz_clear(r);
                return ok;
            }
            mpz_clear(r);
        }
    }
    if (depth > 40) return 0;
    mpz_t g, h; mpz_inits(g, h, NULL);
    if (!rho(z, g)) { mpz_clears(g, h, NULL); return 0; }   /* composite, factors unknown */
    mpz_divexact(h, z, g);
    int a = full_factor(g, F, depth + 1);
    int b = full_factor(h, F, depth + 1);
    mpz_clears(g, h, NULL);
    return a && b;
}

/* Factor c+1 for a component value c. */
static int factor_c1(const mpz_t c, Fac *F) {
    memset(F, 0, sizeof *F);
    mpz_t x; mpz_init(x);
    mpz_add_ui(x, c, 1);
    F->v2 = (int)mpz_scan1(x, 0);
    mpz_tdiv_q_2exp(x, x, F->v2);
    int ok;
    if (mpz_sizeinbase(x, 2) <= 62 && mpz_to_u64(x) <= SPF_MAX) { spf_factor(mpz_to_u64(x), F); ok = 1; }
    else {
        /* trial division by small primes first */
        for (size_t i = 0; i < NPR && PR[i] < 100000; i++) {
            u64 p = PR[i];
            if (mpz_divisible_ui_p(x, p)) {
                int k = 0;
                while (mpz_divisible_ui_p(x, p)) { mpz_divexact_ui(x, x, p); k++; }
                addf(F, p, k);
            }
        }
        ok = full_factor(x, F, 0);
    }
    mpz_clear(x);
    return ok;
}

/* ----------------------------------------------------------- utilities */
static int find(Ent *E, int NE, u64 p) {
    for (int i = 0; i < NE; i++) if (E[i].p == p) return i;
    return -1;
}

static void print_path(FILE *fp, int depth) {
    for (int i = 0; i < depth; i++) fprintf(fp, "%s%llu^%d", i ? "," : "", PATHP[i], PATHE[i]);
}

static int OPEN_R = -1, OPEN_NQ = -1;
static u64 OPEN_Q[64];
static void report_open(int depth, const char *why, const mpz_t lo, const mpz_t hi) {
    NOPEN++;
    printf("OPEN a=%d reason=%s path=", AEXP, why);
    print_path(stdout, depth);
    gmp_printf(" lo=%Zd hi=%Zd", lo, hi);
    if (OPEN_R >= 0) {
        printf(" R=%d nq=%d Q=", OPEN_R, OPEN_NQ);
        for (int i = 0; i < OPEN_NQ && i < 64; i++) printf("%s%llu", i ? ";" : "", OPEN_Q[i]);
    }
    printf("\n");
    fflush(stdout);
}

/* (s, g) with c = s^g, s prime (< 2^62); returns 1 on success.  'composite' verdicts are exact. */
static int prime_power(const mpz_t c, u64 *s, int *g) {
    if (mpz_cmp_ui(c, 3) < 0) return 0;
    mpz_t r; mpz_init(r);
    for (unsigned long k = mpz_sizeinbase(c, 2); k >= 1; k--) {
        int exact = (k == 1) ? (mpz_set(r, c), 1) : mpz_root(r, c, k);
        if (exact && mpz_cmp_ui(r, 2) > 0 && mpz_probab_prime_p(r, 30) > 0) {
            if (mpz_sizeinbase(r, 2) > 62) { mpz_clear(r); return -1; }
            *s = mpz_to_u64(r); *g = (int)k; mpz_clear(r); return 1;
        }
        if (k == 1) break;
    }
    mpz_clear(r);
    return 0;
}

/* ------------------------------------------- 64-bit arithmetic helpers */
static u64 mulmod(u64 a, u64 b, u64 m) { return (u64)((__uint128_t)a * b % m); }
static u64 powmod(u64 a, u64 e, u64 m) {
    u64 r = 1; a %= m;
    while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; }
    return r;
}
/* Miller-Rabin with the first 12 prime bases: deterministic for n < 3.3e24. */
static int mr64(u64 n) {
    if (n < 2) return 0;
    static const u64 B[12] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (int i = 0; i < 12; i++) { if (n == B[i]) return 1; if (n % B[i] == 0) return 0; }
    u64 d = n - 1; int r = 0;
    while ((d & 1) == 0) { d >>= 1; r++; }
    for (int i = 0; i < 12; i++) {
        u64 x = powmod(B[i], d, n);
        if (x == 1 || x == n - 1) continue;
        int comp = 1;
        for (int j = 1; j < r; j++) { x = mulmod(x, x, n); if (x == n - 1) { comp = 0; break; } }
        if (comp) return 0;
    }
    return 1;
}
static u64 iroot(u64 x, int k) {      /* floor(x^(1/k)) */
    u64 r = (u64)pow((double)x, 1.0 / k);
    for (;;) {                          /* correct the floating estimate */
        __uint128_t t = 1; int over = 0;
        for (int i = 0; i < k; i++) { t *= (r + 1); if (t > x) { over = 1; break; } }
        if (!over) { r++; continue; }
        t = 1; over = 0;
        for (int i = 0; i < k; i++) { t *= r; if (t > x) { over = 1; break; } }
        if (over) { r--; continue; }
        return r;
    }
}
/* x > 1: if x = q^k with q prime, set q, k and return 1. */
static int prime_power64(u64 x, u64 *q, int *k) {
    for (int kk = 63; kk >= 2; kk--) {
        if ((x >> kk) == 0 && kk > 1) continue;
        u64 r = iroot(x, kk);
        if (r < 3) continue;
        __uint128_t t = 1;
        for (int i = 0; i < kk; i++) t *= r;
        if (t == x && mr64(r)) { *q = r; *k = kk; return 1; }
    }
    if (mr64(x)) { *q = x; *k = 1; return 1; }
    return 0;
}

static u64 HMAX = 20000000000ULL;      /* new primes are scanned up to HMAX */
static int TIGHT = 0;                  /* UPN_TIGHT=1 uses exact R-factor bound + mod-4 parity */

/* ---------------------------------------------------------------- search */
typedef struct { mpz_t c; u64 p; int e; } Cand;

static void rec(int depth, const mpz_t last, const mpq_t F, int S, Ent *E, int NE);

static void process_fac(int depth, const mpz_t c, u64 p, int e, const mpq_t F, int S,
                        Ent *E, int NE, int nq, const Fac *given);

static void process(int depth, const mpz_t c, u64 p, int e, const mpq_t F, int S,
                    Ent *E, int NE, int nq) {
    process_fac(depth, c, p, e, F, S, E, NE, nq, NULL);
}

static void process_fac(int depth, const mpz_t c, u64 p, int e, const mpq_t F, int S,
                        Ent *E, int NE, int nq, const Fac *given) {
    Fac fc;
    int ok;
    if (given) { fc = *given; ok = 1; } else ok = factor_c1(c, &fc);
    int S2 = S + fc.v2;
    if (S2 > BUD) return;
    int newq = 0, bad = 0;
    for (int i = 0; i < fc.n; i++) {
        int j = find(E, NE, fc.f[i].p);
        if (j >= 0 && E[j].d > 0) { if (E[j].n + fc.f[i].k > E[j].d) { bad = 1; break; } }
        else if (j < 0 || E[j].n == 0) newq++;
    }
    if (bad) return;
    int jp = find(E, NE, p);
    int pinq = (jp >= 0 && E[jp].d == 0 && E[jp].n > 0);
    int Qafter = nq - pinq + newq + fc.unknown_lb;
    if (Qafter > BUD - S2) return;
    mpq_t F2; mpq_init(F2);
    mpq_set_num(F2, c); mpz_add_ui(mpq_denref(F2), c, 0);
    mpz_add_ui(mpq_numref(F2), c, 1); mpz_set(mpq_denref(F2), c);
    mpq_canonicalize(F2);
    mpq_mul(F2, F2, F);
    int cmp = mpq_cmp(F2, TT);
    if (cmp > 0) { mpq_clear(F2); return; }
    if (!ok) {   /* incomplete factorization of c+1 that the cheap tests did not exclude */
        mpz_t z; mpz_init_set_ui(z, 0);
        PATHP[depth] = p; PATHE[depth] = e;
        report_open(depth + 1, "unfactored", c, z);
        mpz_clear(z); mpq_clear(F2);
        return;
    }
    Ent E2[MAXE]; int NE2 = NE;
    memcpy(E2, E, sizeof(Ent) * NE);
    int j = find(E2, NE2, p);
    if (j < 0) { j = NE2++; E2[j].p = p; E2[j].n = 0; }
    E2[j].d = e;
    for (int i = 0; i < fc.n; i++) {
        int k = find(E2, NE2, fc.f[i].p);
        if (k < 0) { k = NE2++; E2[k].p = fc.f[i].p; E2[k].n = 0; E2[k].d = 0; }
        E2[k].n += fc.f[i].k;
    }
    if (NE2 >= MAXE - 70) { fprintf(stderr, "MAXE overflow\n"); exit(2); }
    if (E2[j].n <= e) {
        PATHP[depth] = p; PATHE[depth] = e;
        if (cmp == 0) {
            if (S2 == BUD) {
                int good = 1;
                for (int i = 0; i < NE2; i++) if (E2[i].n != E2[i].d) { good = 0; break; }
                if (good) {
                    NSOL++;
                    printf("SOLUTION a=%d components=", AEXP);
                    print_path(stdout, depth + 1);
                    printf("\n"); fflush(stdout);
                }
            }
        } else {
            rec(depth + 1, c, F2, S2, E2, NE2);
        }
    }
    mpq_clear(F2);
}

static int take(int depth);

/* New primes c in [L, H] (L > LIMIT), found by a segmented sieve.  With f free
 * slots, a new component c may cost v_2(c+1) <= f and bring at most f - v_2(c+1)
 * new primes (Lemma prefix (c)); its other odd prime factors must be present
 * primes within their remaining capacity (Lemma prefix (b)).  These are exactly
 * the tests process() applies, done here with 64-bit arithmetic; a candidate that
 * passes is handed to process_fac() together with its complete factorization. */
static long long SCANNED = 0;
static void scan_big(int depth, u64 L, u64 H, const mpq_t F, int S, Ent *E, int NE, int nq, int f) {
    u64 ap[MAXE]; int acap[MAXE]; int na = 0;
    for (int i = 0; i < NE; i++) {
        if (E[i].d > 0) { ap[na] = E[i].p; acap[na] = E[i].d - E[i].n; na++; }
        else if (E[i].n > 0) { ap[na] = E[i].p; acap[na] = -1; na++; }
    }
    const u64 SEG = 1ULL << 24;
    uint8_t *comp = malloc(SEG);
    mpz_t c; mpz_init(c);
    if ((L & 1) == 0) L++;
    for (u64 base = L; base <= H; base += SEG) {
        u64 end = base + SEG - 1; if (end > H) end = H;
        memset(comp, 0, end - base + 1);
        for (size_t k = 0; k < NPR; k++) {
            u64 p = PR[k];
            if (p * p > end) break;
            u64 st = ((base + p - 1) / p) * p;
            if (st < p * p) st = p * p;
            for (u64 x = st; x <= end; x += p) comp[x - base] = 1;
        }
        for (u64 n = base; n <= end; n += 2) {
            if (comp[n - base]) continue;
            SCANNED++;
            if (find(E, NE, n) >= 0) continue;
            if (!take(depth)) continue;
            u64 x = n + 1;
            int cost = __builtin_ctzll(x);
            if (cost > f) continue;
            x >>= cost;
            int newmax = f - cost;
            Fac fc; memset(&fc, 0, sizeof fc); fc.v2 = cost;
            int bad = 0;
            for (int i = 0; i < na && x > 1; i++) {
                if (x % ap[i] == 0) {
                    int k = 0;
                    while (x % ap[i] == 0) { x /= ap[i]; k++; }
                    if (acap[i] >= 0 && k > acap[i]) { bad = 1; break; }
                    addf(&fc, ap[i], k);
                }
            }
            if (bad) continue;
            if (x > 1) {
                if (newmax < 1) continue;
                u64 q; int k;
                if (prime_power64(x, &q, &k)) addf(&fc, q, k);
                else if (newmax >= 2) {          /* two or more new primes: factor fully */
                    mpz_set_ui(c, n);
                    process(depth, c, n, 1, F, S, E, NE, nq);
                    continue;
                } else continue;                 /* x has >= 2 distinct primes, only 1 allowed */
            }
            mpz_set_ui(c, n);
            process_fac(depth, c, n, 1, F, S, E, NE, nq, &fc);
        }
    }
    mpz_clear(c);
    free(comp);
}

static int take(int depth) {
    if (depth != SPLITD) return 1;
    return (SPLITCNT++ % NWORK) == WORKER;
}

static double logG_hint(double x, const double *mq, int nq, int exp_free) {
    double s = log1p(1.0 / x);
    for (int i = 0; i < nq; i++) { double m = mq[i] > x + 1 ? mq[i] : x + 1; s += log1p(1.0 / m); }
    return s + exp_free * log1p(1.0 / (x + 1));
}

/* Exact rational check of F * (1 + 1/x)^R >= TT and F * G(x) >= TT. */
static int geq_G_mpz(const mpz_t x, const mpq_t F, int R, const mpz_t *mq, int nq, int exp_free) {
    mpz_t lhs, rhs, xp1, xp2, u, up1;
    mpz_inits(lhs, rhs, xp1, xp2, u, up1, NULL);
    mpz_add_ui(xp1, x, 1);
    if (!TIGHT) {
        mpz_pow_ui(u, xp1, (unsigned long)R);
        mpz_mul(lhs, mpq_numref(F), mpq_denref(TT));
        mpz_mul(lhs, lhs, u);
        mpz_pow_ui(u, x, (unsigned long)R);
        mpz_mul(rhs, mpq_numref(TT), mpq_denref(F));
        mpz_mul(rhs, rhs, u);
        if (mpz_cmp(lhs, rhs) < 0) {
            mpz_clears(lhs, rhs, xp1, xp2, u, up1, NULL);
            return 0;
        }
    }
    mpz_add_ui(xp2, x, 2);
    mpz_mul(lhs, mpq_numref(F), mpq_denref(TT));
    mpz_mul(lhs, lhs, xp1);
    mpz_mul(rhs, mpq_numref(TT), mpq_denref(F));
    mpz_mul(rhs, rhs, x);
    if (exp_free > 0) {
        mpz_pow_ui(u, xp2, (unsigned long)exp_free);
        mpz_mul(lhs, lhs, u);
        mpz_pow_ui(u, xp1, (unsigned long)exp_free);
        mpz_mul(rhs, rhs, u);
    } else if (exp_free < 0) {   /* TIGHT with free_ == 0 */
        mpz_mul(lhs, lhs, xp1);
        mpz_mul(rhs, rhs, xp2);
    }
    for (int i = 0; i < nq; i++) {
        if (mpz_cmp(xp1, mq[i]) >= 0) mpz_set(u, xp1);
        else mpz_set(u, mq[i]);
        mpz_add_ui(up1, u, 1);
        mpz_mul(lhs, lhs, up1);
        mpz_mul(rhs, rhs, u);
    }
    int res = mpz_cmp(lhs, rhs) >= 0;
    mpz_clears(lhs, rhs, xp1, xp2, u, up1, NULL);
    return res;
}

static void rec(int depth, const mpz_t last, const mpq_t F, int S, Ent *E, int NE) {
    NODES++;
    int R = BUD - S;
    if (R <= 0) return;
    u64 Qp[MAXE]; int nq = 0;
    for (int i = 0; i < NE; i++) if (E[i].n > 0 && E[i].d == 0) Qp[nq++] = E[i].p;
    int free_ = R - nq;
    if (free_ < 0) return;

    mpq_t delta, inv; mpq_inits(delta, inv, NULL);
    mpq_sub(delta, TT, F); mpq_div(delta, delta, F);          /* T/F - 1 */
    mpq_inv(inv, delta);
    mpz_t lo, hi, t; mpz_inits(lo, hi, t, NULL);
    mpz_cdiv_q(lo, mpq_numref(inv), mpq_denref(inv));         /* exact ceil(1/delta) */
    mpz_add_ui(t, last, 1);
    if (mpz_cmp(t, lo) > 0) mpz_set(lo, t);

    {
        /* m_q: least admissible power of q above last (even exponent if TIGHT, free_==0, q=3 mod 4) */
        double mqd[MAXE];
        mpz_t *mq = malloc(sizeof(mpz_t) * (nq + 1));
        int *mqe = malloc(sizeof(int) * (nq + 1));
        for (int i = 0; i < nq; i++) {
            int j = find(E, NE, Qp[i]);
            int e = E[j].n > 1 ? E[j].n : 1;
            mpz_init(mq[i]); mpz_ui_pow_ui(mq[i], Qp[i], e);
            while (mpz_cmp(mq[i], last) <= 0) { mpz_mul_ui(mq[i], mq[i], Qp[i]); e++; }
            if (TIGHT && free_ == 0 && (Qp[i] & 3ULL) == 3ULL && (e & 1)) {
                mpz_mul_ui(mq[i], mq[i], Qp[i]);
                e++;
            }
            mqe[i] = e;
            mqd[i] = mpz_get_d(mq[i]);
        }
        if (TIGHT && free_ == 0 && nq > 0) {
            mpz_set(t, mq[0]);
            for (int i = 1; i < nq; i++) if (mpz_cmp(mq[i], t) < 0) mpz_set(t, mq[i]);
            if (mpz_cmp(t, lo) > 0) mpz_set(lo, t);
        }
        int exp_free = TIGHT ? (free_ - 1) : free_;

        /* 100% exact prune: if F * G(lo) < TT, no x >= lo can reach TT. */
        if (!geq_G_mpz(lo, F, R, mq, nq, exp_free)) {
            for (int i = 0; i < nq; i++) mpz_clear(mq[i]);
            free(mq); free(mqe);
            goto done;
        }

        /* Floating-point hint ONLY for initial bracket; exact geq_G_mpz certifies hi. */
        mpz_t guess, a2, b2, mid, step;
        mpz_inits(guess, a2, b2, mid, step, NULL);
        mpz_set(guess, lo);
        double d = mpq_get_d(delta);
        if (d > 1e-14 && isfinite(d)) {
            double U = 1.0 / expm1(log1p(d) / R);
            if (isfinite(U) && U > 0.0 && U < 9e18) {
                mpz_set_d(b2, U); mpz_add_ui(b2, b2, 2);
                if (mpz_cmp(b2, lo) > 0 && mpz_sizeinbase(b2, 2) < 63) {
                    u64 ah = mpz_get_ui(lo), bh = mpz_get_ui(b2);
                    double target = log1p(d);
                    if (logG_hint((double)ah, mqd, nq, exp_free) >= target) {
                        while (bh - ah > 1) {
                            u64 mh = ah + (bh - ah) / 2;
                            if (logG_hint((double)mh, mqd, nq, exp_free) >= target) ah = mh; else bh = mh;
                        }
                        mpz_set_ui(guess, ah);
                    }
                }
            }
        }

        if (geq_G_mpz(guess, F, R, mq, nq, exp_free)) {
            mpz_set(a2, guess);
            mpz_add_ui(t, a2, 1);
            if (!geq_G_mpz(t, F, R, mq, nq, exp_free)) {
                mpz_set(hi, a2);
            } else {
                mpz_set_ui(step, 2);
                for (;;) {
                    mpz_add(t, a2, step);
                    if (!geq_G_mpz(t, F, R, mq, nq, exp_free)) break;
                    mpz_set(a2, t);
                    mpz_mul_2exp(step, step, 1);
                }
                mpz_add(b2, a2, step);
                for (;;) {
                    mpz_sub(mid, b2, a2);
                    if (mpz_cmp_ui(mid, 1) <= 0) break;
                    mpz_fdiv_q_2exp(mid, mid, 1); mpz_add(mid, mid, a2);
                    if (geq_G_mpz(mid, F, R, mq, nq, exp_free)) mpz_set(a2, mid); else mpz_set(b2, mid);
                }
                mpz_set(hi, a2);
            }
        } else {
            mpz_set(b2, guess);
            mpz_sub_ui(t, b2, 1);
            if (geq_G_mpz(t, F, R, mq, nq, exp_free)) {
                mpz_set(hi, t);
            } else {
                mpz_set_ui(step, 2);
                for (;;) {
                    mpz_sub(t, b2, step);
                    if (mpz_cmp(t, lo) <= 0 || geq_G_mpz(t, F, R, mq, nq, exp_free)) break;
                    mpz_set(b2, t);
                    mpz_mul_2exp(step, step, 1);
                }
                mpz_sub(a2, b2, step);
                if (mpz_cmp(a2, lo) < 0) mpz_set(a2, lo);
                for (;;) {
                    mpz_sub(mid, b2, a2);
                    if (mpz_cmp_ui(mid, 1) <= 0) break;
                    mpz_fdiv_q_2exp(mid, mid, 1); mpz_add(mid, mid, a2);
                    if (geq_G_mpz(mid, F, R, mq, nq, exp_free)) mpz_set(a2, mid); else mpz_set(b2, mid);
                }
                mpz_set(hi, a2);
            }
        }
        mpz_clears(guess, a2, b2, mid, step, NULL);

        /* (a) powers of demanded primes */
        mpz_t c; mpz_init(c);
        for (int i = 0; i < nq; i++) {
            int e = mqe[i];
            mpz_set(c, mq[i]);
            while (mpz_cmp(c, hi) <= 0) {
                if (mpz_cmp(c, lo) >= 0 && take(depth)) process(depth, c, Qp[i], e, F, S, E, NE, nq);
                mpz_mul_ui(c, c, Qp[i]); e++;
            }
        }
        for (int i = 0; i < nq; i++) mpz_clear(mq[i]);
        free(mq); free(mqe);

        if (free_ == 1) {
            /* (b) smooth shortcut: c + 1 = 2M, M odd and composed of present primes */
            u64 ap[MAXE]; int acap[MAXE]; int na = 0;
            for (int i = 0; i < NE; i++) {
                if (E[i].d > 0 && E[i].d - E[i].n > 0) { ap[na] = E[i].p; acap[na] = E[i].d - E[i].n; na++; }
            }
            for (int i = 0; i < NE; i++) if (E[i].n > 0 && E[i].d == 0) { ap[na] = E[i].p; acap[na] = -1; na++; }
            mpz_t mlo, mhi; mpz_inits(mlo, mhi, NULL);
            mpz_add_ui(mlo, lo, 2); mpz_fdiv_q_2exp(mlo, mlo, 1);
            mpz_add_ui(mhi, hi, 1); mpz_fdiv_q_2exp(mhi, mhi, 1);
            size_t cap = 1024, cnt = 1;
            mpz_t *sm = malloc(sizeof(mpz_t) * cap);
            mpz_init_set_ui(sm[0], 1);
            for (int i = 0; i < na; i++) {
                size_t base = cnt;
                for (size_t k = 0; k < base; k++) {
                    mpz_t x; mpz_init_set(x, sm[k]);
                    for (int kk = 1;; kk++) {
                        mpz_mul_ui(x, x, ap[i]);
                        if (mpz_cmp(x, mhi) > 0 || (acap[i] >= 0 && kk > acap[i])) break;
                        if (cnt == cap) { cap *= 2; sm = realloc(sm, sizeof(mpz_t) * cap); }
                        mpz_init_set(sm[cnt++], x);
                    }
                    mpz_clear(x);
                }
            }
            for (size_t k = 0; k < cnt; k++) {
                if (mpz_cmp(sm[k], mlo) >= 0) {
                    mpz_mul_2exp(c, sm[k], 1); mpz_sub_ui(c, c, 1);
                    u64 s; int g;
                    int pp = prime_power(c, &s, &g);
                    if (pp < 0) { report_open(depth, "bigprime", c, c); }
                    else if (pp == 1 && find(E, NE, s) < 0) {
                        if (take(depth)) process(depth, c, s, g, F, S, E, NE, nq);
                    }
                }
                mpz_clear(sm[k]);
            }
            free(sm);
            mpz_clears(mlo, mhi, NULL);
        } else if (free_ >= 2) {
            /* (c) new primes by scanning, and their powers */
            if (mpz_cmp_ui(hi, LIMIT) > 0) {
                if (mpz_sizeinbase(hi, 2) > 62 || mpz_get_ui(hi) > HMAX || mpz_get_ui(hi) / LIMIT > LIMIT) {
                    OPEN_R = R; OPEN_NQ = nq;
                    for (int i = 0; i < nq && i < 64; i++) OPEN_Q[i] = Qp[i];
                    report_open(depth, "range", lo, hi);
                    OPEN_R = -1;
                } else {
                    u64 L = mpz_sizeinbase(lo, 2) < 63 ? mpz_get_ui(lo) : (u64)-1;
                    if (L < LIMIT + 1) L = LIMIT + 1;
                    scan_big(depth, L, mpz_get_ui(hi), F, S, E, NE, nq, free_);
                }
            }
            u64 lo64 = mpz_sizeinbase(lo, 2) < 63 ? mpz_get_ui(lo) : (u64)-1;
            u64 hiS = mpz_cmp_ui(hi, LIMIT) > 0 ? LIMIT : mpz_get_ui(hi);
            if (lo64 <= hiS) {
                size_t lo_i = 0, hi_i = NPR;
                while (lo_i < hi_i) { size_t m = (lo_i + hi_i) / 2; if (PR[m] < lo64) lo_i = m + 1; else hi_i = m; }
                for (size_t k = lo_i; k < NPR && PR[k] <= hiS; k++) {
                    u64 s = PR[k];
                    if (find(E, NE, s) >= 0) continue;
                    if (!take(depth)) continue;
                    mpz_set_ui(c, s);
                    process(depth, c, s, 1, F, S, E, NE, nq);
                }
            }
            for (size_t k = 0; k < NPR; k++) {
                u64 s = PR[k];
                mpz_set_ui(c, s); mpz_mul_ui(c, c, s);
                if (mpz_cmp(c, hi) > 0) break;
                if (find(E, NE, s) >= 0) continue;
                for (int g = 2; mpz_cmp(c, hi) <= 0; g++) {
                    if (mpz_cmp(c, lo) >= 0 && take(depth)) process(depth, c, s, g, F, S, E, NE, nq);
                    mpz_mul_ui(c, c, s);
                }
            }
        }
        mpz_clear(c);
    }
done:
    mpz_clears(lo, hi, t, NULL);
    mpq_clears(delta, inv, NULL);
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: upn a factorization [worker nworkers splitdepth]\n"); return 1; }
    AEXP = atoi(argv[1]); BUD = AEXP + 1;
    if (argc >= 6) { WORKER = atoi(argv[3]); NWORK = atoi(argv[4]); SPLITD = atoi(argv[5]); }
    if (getenv("UPN_HMAX")) HMAX = strtoull(getenv("UPN_HMAX"), NULL, 10);
    if (getenv("UPN_TIGHT") && strcmp(getenv("UPN_TIGHT"), "0") != 0) TIGHT = 1;
    build_tables();
    mpq_init(TT);
    mpz_t num, den; mpz_inits(num, den, NULL);
    mpz_ui_pow_ui(num, 2, AEXP + 1);
    mpz_ui_pow_ui(den, 2, AEXP); mpz_add_ui(den, den, 1);
    mpq_set_num(TT, num); mpq_set_den(TT, den); mpq_canonicalize(TT);

    /* seed factorization, checked against 2^a+1 */
    Ent E[MAXE]; int NE = 0;
    mpz_t prod; mpz_init_set_ui(prod, 1);
    char *s = strdup(argv[2]);
    for (char *tok = strtok(s, ","); tok; tok = strtok(NULL, ",")) {
        u64 p = strtoull(tok, NULL, 10);
        char *colon = strchr(tok, ':');
        int k = colon ? atoi(colon + 1) : 1;
        if (p >= BIGP || !is_prime_u64(p)) { fprintf(stderr, "seed factor %llu not accepted\n", p); return 3; }
        E[NE].p = p; E[NE].n = k; E[NE].d = 0; NE++;
        for (int i = 0; i < k; i++) mpz_mul_ui(prod, prod, p);
    }
    if (mpz_cmp(prod, den) != 0) { fprintf(stderr, "seed factorization does not multiply to 2^a+1\n"); return 3; }

    mpq_t F; mpq_init(F); mpq_set_ui(F, 1, 1);
    mpz_t last; mpz_init_set_ui(last, 1);
    int S0 = 0;
    /* Optional forced components (UPN_FORCE="p:e,..."): they are exact components of
     * n, placed before the search; the remaining components are searched as usual. */
    int depth0 = 0;
    if (getenv("UPN_FORCE")) {
        char *fs = strdup(getenv("UPN_FORCE"));
        for (char *tok = strtok(fs, ","); tok; tok = strtok(NULL, ",")) {
            u64 p = strtoull(tok, NULL, 10);
            char *colon = strchr(tok, ':');
            int e = colon ? atoi(colon + 1) : 1;
            if (p >= BIGP || !is_prime_u64(p) || p == 2) { fprintf(stderr, "bad forced prime\n"); return 3; }
            mpz_t c; mpz_init(c); mpz_ui_pow_ui(c, p, e);
            Fac fc; if (!factor_c1(c, &fc)) { fprintf(stderr, "cannot factor forced component\n"); return 3; }
            int j = find(E, NE, p);
            if (j < 0) { j = NE++; E[j].p = p; E[j].n = 0; }
            E[j].d = e;
            for (int i = 0; i < fc.n; i++) {
                int k = find(E, NE, fc.f[i].p);
                if (k < 0) { k = NE++; E[k].p = fc.f[i].p; E[k].n = 0; E[k].d = 0; }
                E[k].n += fc.f[i].k;
            }
            S0 += fc.v2;
            mpq_t r; mpq_init(r);
            mpz_add_ui(mpq_numref(r), c, 1); mpz_set(mpq_denref(r), c); mpq_canonicalize(r);
            mpq_mul(F, F, r); mpq_clear(r); mpz_clear(c);
            PATHP[depth0] = p; PATHE[depth0] = e; depth0++;
        }
        int consistent = (S0 <= BUD) && (mpq_cmp(F, TT) <= 0);
        for (int i = 0; i < NE; i++) if (E[i].d > 0 && E[i].n > E[i].d) consistent = 0;
        if (!consistent) {
            printf("DONE a=%d forced=%s nodes=0 solutions=0 open=0 note=forced-components-inconsistent\n", AEXP, getenv("UPN_FORCE"));
            return 0;
        }
        if (mpq_cmp(F, TT) == 0) {
            int good = (S0 == BUD);
            for (int i = 0; i < NE; i++) if (E[i].n != E[i].d) good = 0;
            if (good) { printf("SOLUTION a=%d components=", AEXP); print_path(stdout, depth0); printf("\n"); }
            printf("DONE a=%d forced=%s nodes=0 solutions=%d open=0\n", AEXP, getenv("UPN_FORCE"), good);
            return 0;
        }
    }
    clock_t t0 = clock();
    /* Resuming an open node: UPN_FORCE holds the node's prefix and UPN_LAST its largest
     * component, so that exactly the subtree of that node is searched. */
    if (getenv("UPN_LAST")) mpz_set_str(last, getenv("UPN_LAST"), 10);
    if (SPLITD >= 0) SPLITD += depth0;
    rec(depth0, last, F, S0, E, NE);
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    printf("DONE a=%d worker=%d/%d nodes=%lld solutions=%lld open=%lld scanned_big=%lld hmax=%llu seconds=%.1f\n",
           AEXP, WORKER, NWORK, NODES, NSOL, NOPEN, SCANNED, HMAX, secs);
    return 0;
}
