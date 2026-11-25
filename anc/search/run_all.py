#!/usr/bin/env python3
"""Run the exact search (upn.c) for a list of exponents a and summarise the results.

    cc -O2 -o upn upn.c -lgmp -lm        (add -I/-L flags for your GMP install)
    python3 run_all.py 1-15,17,19,21,23   [--jobs 8] [--log logs]

For each a it factors 2^a+1 (sympy), proves every prime factor above 2^64 with
PARI/GP if available, runs ./upn a <factorization>, and records the result.  A
value of a counts as settled only if the run ends with "open=0".
"""
import argparse
import concurrent.futures as cf
import os
import re
import shutil
import subprocess
import sys

from sympy import factorint


def parse(spec):
    out = []
    for part in spec.split(","):
        if "-" in part:
            lo, hi = map(int, part.split("-"))
            out += list(range(lo, hi + 1))
        else:
            out.append(int(part))
    return out


def prove_prime(p):
    if p < 2 ** 64:
        return "BPSW (deterministic below 2^64)"
    if not shutil.which("gp"):
        return "UNPROVEN (PARI/GP not found)"
    r = subprocess.run(["gp", "-q"], input=f"isprime({p})\n", capture_output=True, text=True).stdout.strip()
    return "proved by PARI/GP isprime" if r == "1" else f"NOT PRIME? gp says {r}"


def run(a, logdir, force):
    f = factorint(2 ** a + 1)
    seed = ",".join(f"{p}:{k}" for p, k in sorted(f.items()))
    notes = [f"{p}: {prove_prime(p)}" for p in f if p >= 2 ** 64]
    env = dict(os.environ)
    if force:
        env["UPN_FORCE"] = force
    out = subprocess.run(["./upn", str(a), seed], capture_output=True, text=True, env=env).stdout
    name = f"a{a}" + (f"_force_{force.replace(':', '^').replace(',', '_')}" if force else "")
    with open(os.path.join(logdir, name + ".log"), "w") as fh:
        fh.write(f"# a={a} seed 2^a+1 = {seed}\n")
        for n in notes:
            fh.write(f"# seed prime {n}\n")
        fh.write(out)
    done = [l for l in out.splitlines() if l.startswith("DONE")]
    sols = [l.split("components=")[1] for l in out.splitlines() if l.startswith("SOLUTION")]
    m = re.search(r"nodes=(\d+).*open=(\d+)", done[-1]) if done else None
    return a, (int(m.group(1)) if m else None), (int(m.group(2)) if m else None), sols, notes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("avalues")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--log", default="logs")
    ap.add_argument("--force", default="", help="forced components, e.g. 3:4,41:1")
    args = ap.parse_args()
    os.makedirs(args.log, exist_ok=True)
    with cf.ThreadPoolExecutor(args.jobs) as ex:
        results = sorted(ex.map(lambda a: run(a, args.log, args.force), parse(args.avalues)))
    ok = True
    for a, nodes, opens, sols, notes in results:
        status = "SETTLED" if opens == 0 else "NOT SETTLED"
        ok &= opens == 0
        print(f"a={a:4d}  nodes={nodes}  open={opens}  {status}  solutions={sols or 'none'}")
        for n in notes:
            print(f"        seed prime {n}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
