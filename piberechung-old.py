#!/usr/bin/env python3
"""
piberechung-old.py – Pi via Chudnovsky-Algorithmus, einfache iterative Version.

Kein Binary Splitting, kein gmpy2, kein Multiprocessing.
Nur Python-Bordmittel (decimal).

Vergleichsscript zu piberechung.py.
"""

import argparse
import sys
import time
from decimal import Decimal, getcontext


def berechne_pi(stellen: int) -> str:
    getcontext().prec = stellen + 20

    C = 426880 * Decimal(10005).sqrt()
    K = Decimal(6)
    M = Decimal(1)
    X = Decimal(1)
    S = Decimal(13591409)

    for k in range(1, stellen // 14 + 2):
        M = M * (K ** 3 - 16 * K) / (k ** 3)
        X *= -262537412640768000
        S += M * (13591409 + 545140134 * k) / X
        K += 12

    pi = C / S
    getcontext().prec = stellen + 1
    return str(+pi)


def main():
    parser = argparse.ArgumentParser(
        description="Berechnet Pi (einfache iterative Version, nur decimal)"
    )
    parser.add_argument(
        "stellen", nargs="?", type=int, default=100,
        help="Anzahl der Dezimalstellen (Standard: 100)"
    )
    parser.add_argument(
        "-b", "--benchmark", action="store_true",
        help="Zeigt Laufzeit an"
    )
    args = parser.parse_args()

    if args.stellen < 1:
        print("Fehler: Mindestens 1 Stelle erforderlich.", file=sys.stderr)
        sys.exit(1)

    start = time.perf_counter()
    pi_str = berechne_pi(args.stellen)
    dauer = time.perf_counter() - start

    print(pi_str[:args.stellen + 2])

    if args.benchmark:
        print(f"\n--- Benchmark ---", file=sys.stderr)
        print(f"Stellen:  {args.stellen}", file=sys.stderr)
        print(f"Backend:  decimal (iterativ, kein Binary Splitting)", file=sys.stderr)
        print(f"Laufzeit: {dauer:.4f}s", file=sys.stderr)


if __name__ == "__main__":
    main()
