# APR-CL primality-verification transcripts for Lemma 23

For each of the six m values in Lemma 23 (Deep-Pratt closures), this directory contains:
- prp_m{m}.txt: the decimal representation of the prime divisor p* | 2^m + 1
- check_m{m}.gp: the PARI/GP script that ran the APR-CL proof (`isprime(n, 2)`)
- aprcl_done_m{m}.txt: contains `1` on success
- aprcl_proofs.log: timing log

| m     | digits | wall-clock | result |
|-------|--------|-----------:|-------:|
|  2446 |   368  |    3.4 s   |   1    |
| 10294 |  1549  |     79 s   |   1    |
| 10958 |  1649  |     88 s   |   1    |
| 17398 |  2612  |    615 s   |   1    |
| 19066 |  2870  |    679 s   |   1    |
| 20282 |  3053  |    958 s   |   1    |

PARI/GP version: 2.17.3 (released Nov 2025, GMP 6.3.0 backend).
Date of proof: 2026-05-18.
Hardware: Apple Silicon, multi-threaded (11-core parallel).

To reproduce a single proof:
```
echo 'print(isprime(eval(readstr("prp_m10294.txt")[1]), 2));' | gp -q \
  --default parisizemax=8000000000 \
  --default threadsizemax=4000000000
```
For the larger candidates (m >= 17398), use:
```
--default parisizemax=32000000000 --default threadsizemax=16000000000
```

## Additional ECPP certificates (downloaded from FactorDB)

For two additional medium-large PRPs that did not yield deep-Pratt closures
(no non-3-Higgs witness was found in p* - 1 below 10^9), Primo-format ECPP
certificates already exist in FactorDB and are bundled here for completeness:

| m     | p* digits | source                                                | file                  |
|-------|-----------|-------------------------------------------------------|-----------------------|
|  4166 |    620    | FactorDB id 1100000000000442093 (Primo 3.0.9)         | primo_cert_m4166.txt  |
| 15758 |   2358    | FactorDB id 1100000000619254778 (Primo 4.0.4-LX64, uploaded 2013-11-17) | primo_cert_m15758.txt |

Both certificates are status "Verified" on FactorDB. Readers can run
Primo's verifier or any compatible ECPP-certificate checker on these
files to independently confirm primality of each p*.

These two primes are NOT currently used in Lemma 23 (the deep-Pratt sweep
up to 10^9 found no non-3-Higgs witness in either p* - 1); they are
included so that any future closure that does find such a witness becomes
rigorous immediately, without the need for new primality work.

