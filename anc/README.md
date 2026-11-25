# Ancillary files (version 3)

This directory contains the non-TeX material used by the manuscript.

## New in version 3: `search/`

The exact search of Section 3, which proves Theorem 1 and the proposition on the
ten kernel candidates in Section 4 without any hypothesis. Every pruning decision is
made in exact integer arithmetic.

- `search2.py`: the search in Python with exact rational arithmetic (needs `sympy`, `numpy`).
- `upn.c`: the same search in C with GMP. The two programs visit the same number of
  nodes on every value of `a` where both were run.
- `run_all.py`: factors `2^a+1` (with `sympy`), proves any prime factor above `2^64` with PARI/GP,
  runs `./upn` for each `a`, and reports `SETTLED` only when no branch is open.
- `logs/`: the output for every `a` in Theorem 1 (`aN.log`), and for the kernel runs
  of Section 4 (`logs/kernel/`). The last line of each C log is
  `DONE ... nodes=... solutions=... open=0`. The case `a = 117` was run with the
  Python program, because 2^117+1 has a prime factor above 2^64; its log ends
  `OPEN=0`.

To rerun:

```bash
cd anc/search
cc -O2 -o upn upn.c -lgmp -lm          # add -I/-L flags for your GMP install
python3 run_all.py 1-15,17,19,21,23,25 --jobs 4 --log mylogs
UPN_FORCE="3:4,41:1" ./upn 45 "3:3,11:1,19:1,331:1,18837001:1"   # one kernel run
```

The largest runs, `a = 16` and `a = 25`, take hours on one core. For `a = 16` the
scan limit for new primes is raised from 2e10 to 3e10, so that no branch is left open:

```bash
UPN_HMAX=30000000000 python3 run_all.py 16 --log mylogs
UPN_FORCE="3:4,41:1" python3 search2.py 117                        # the Python-only kernel run
```

## Material from version 2

The rest of this directory is unchanged from version 2. It supports Sections 4
and 5; the impostor-kernel certificate (Filters Z and N) is valid only under
Hypothesis H3, and its Filter O is withdrawn.

## Contents

- `factor_cache.json`: JSON map `m -> {prime_factor: exponent}` for known prime factors used by the verification scripts. Some entries are partial factorizations; the scripts treat partial entries conservatively.
- `data/`: current frontier list through `m <= 50000`, the through-40000 community list, and notes about the `40000 < m <= 50000` retry pass.
- `scripts/`: Python verification scripts and their shared research-code dependency.
- `aprcl_transcripts/`: PARI/GP inputs, decimal prime records, success markers, and logs for the APR-CL primality-verification runs used in the deep-Pratt lemma.

## Smoke test

From this directory's parent source root:

```bash
cd anc/scripts
python3 verify_Heven_rigorous.py --k-min 1 --k-max 600 --deep-secs 1
```

Expected final line:

```text
Verified H_even ∩ [2, 1200]: [2, 6, 10, 18, 26, 30, 46, 62, 82, 122]
```

## Candidate frontier

The current undecided frontier is `data/candidate_frontier_50000.tsv`. It has 262 rows, matching the bound in the paper: 10 verified elements of `H_even` through 1200 plus 262 undecided candidates gives `|H_even ∩ [2,50000]| <= 272`.

The `candidates_for_forum.*` files contain the current through-40000 unresolved community list after removing the closed case `m=2446`.

## Main impostor-kernel certificate

```bash
cd anc/scripts
python3 impostor_certificate.py \
    --max-a 10000 --factor-limit 2000 --factor-timeout 3.0 \
    --max-rounds 200 --max-bases 8000 --max-exp 5000
```

The script uses `../factor_cache.json` by default. Set `FACTOR_CACHE=/path/to/cache.json` to override this. The full command is intentionally much longer than the smoke test.
A completed run log is bundled at `data/full_certificate_run.log`.
Use `--use-factordb` only when intentionally extending the cache with live FactorDB lookups.

## APR-CL transcripts

The APR-CL directory contains reproducibility transcripts, not standalone certificate objects. The PARI/GP scripts call `isprime(n, 2)` on the recorded decimal primes; rerunning them with the same PARI/GP version should reproduce the success markers.
