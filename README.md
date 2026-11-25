# Unitary perfect numbers with a small power of two, and a conditional bounded-box reduction

**Tom Maciejewski** · Preprint, 2026 · [Read the paper (PDF)](paper.pdf)

## Abstract

A unitary perfect number is a positive integer $n$ with $\sigma^*(n)=2n$, where $\sigma^*(n)$ is the sum of the divisors $d$ of $n$ with $\gcd(d,n/d)=1$. Only five are known, and it has been conjectured that there are only finitely many. Every unitary perfect number has the form $n=2^am$ with $m$ odd and $a\ge1$. We prove that the only unitary perfect numbers with $a\le 17$ or $a\in\{19, 21, 23, 25\}$ are $6$, $60$, $90$ and $87360$; so a sixth unitary perfect number has $a\ge 18$ and $a\notin\{19, 21, 23, 25\}$. The previously published bound was $a\ge 11$. The proof is an exhaustive search over the odd prime-power components of $m$, taken in increasing order, which uses nothing but the balance $(2^a+1)\prod(p^e+1)=2^{a+1}\prod p^e$: every odd prime must occur to the same power on both sides, and the power of $2$ allows at most $a+1$ components. The search recovers the four known unitary perfect numbers with $a\le 6$, and with the component $5^4$ fixed it finds Wall's $24$-digit example at $a=18$. This version corrects the earlier ones. Version 1 asserted, and Proposition 1 of version 2 claimed to prove, that every prime divisor of a unitary perfect number is a 3-Higgs prime. That proof is invalid and the statement is open, so the bounded-box elimination of “impostor” kernels for $a\le 10000$ claimed there holds only under it. We state the hypothesis explicitly, repair a second gap (a 2-adic budget step assumed that $v_2(p^e+1)$ cannot decrease as $e$ grows), settle unconditionally every case in which that step had been used, and withdraw several conjectures and claims of version 2 about the numbers $2^m+1$ that were false or unsupported.

## Contents

| File | What it is |
|---|---|
| [`paper.pdf`](paper.pdf) | The paper (16 pages) |
| [`paper.tex`](paper.tex) | Its LaTeX source |
| [`anc/README.md`](anc/README.md) | Guide to the ancillary files and how to rerun them. |
| [`anc/search/`](anc/search/) | The exact search in C and Python, with the log of every run behind the theorem. |
| [`anc/scripts/`](anc/scripts/) | Python verification scripts for the conditional results and for numbers of the form 2^m+1. |
| [`anc/data/`](anc/data/) | Candidate lists and the log of the full certificate run. |
| [`anc/factor_cache.json`](anc/factor_cache.json) | Known prime factors used by the scripts. |
| [`anc/aprcl_transcripts/`](anc/aprcl_transcripts/) | PARI/GP primality-proof transcripts. |

## Citing this paper

```bibtex
@misc{maciejewski2026upn,
  author = {Maciejewski, Tom},
  title  = {Unitary perfect numbers with a small power of two, and a conditional bounded-box reduction},
  year   = {2026},
  note   = {Preprint},
  url    = {https://github.com/tommac/unitary-perfect-numbers},
}
```

## Status

This is a preprint and has not been peer reviewed. If you find an error, please open an issue.

## License

The paper is licensed under [CC BY 4.0](LICENSE). The code in `anc/` is under the [MIT License](LICENSE-CODE).
