#!/usr/bin/env python3
"""
pi_berechnung.py – Pi via Chudnovsky-Algorithmus mit Binary Splitting.

Backends (Priorität):
  1. mpmath + gmpy2  – schnellstes (GMP-Arithmetik, C-Implementierung)
  2. mpmath          – schnell dank internem Binary Splitting in C
  3. decimal         – reines Python, Binary Splitting statt Iteration

Binary Splitting transformiert die Chudnovsky-Reihe in eine einzige
Ganzzahl-Berechnung und vermeidet teure Hochpräzisions-Gleitkommaarithmetik
pro Schritt. Alle Zwischenwerte sind Integer; nur ein einziges sqrt am Ende.

  pi = 426880 * sqrt(10005) * Q / T
  (P, Q, T) via rekursivem Teile-und-Herrsche über Termindizes

Installieren für maximale Geschwindigkeit:
  pip install mpmath gmpy2
"""

import argparse
import math
import sys
import time
from decimal import Decimal, getcontext

# ── Konstante ──────────────────────────────────────────────────────────────────
_A = 13591409
_B = 545140134
_C3_OVER_24 = 640320 ** 3 // 24  # 10939058860032000


# ── Binary Splitting (reines Python) ──────────────────────────────────────────

def _bs(a: int, b: int) -> tuple:
    """Liefert (P, Q, T) für Chudnovsky-Terme [a, b)."""
    if b - a == 1:
        if a == 0:
            P, Q = 1, 1
        else:
            P = (6*a - 5) * (2*a - 1) * (6*a - 1)
            Q = _C3_OVER_24 * a * a * a
        T = P * (_A + _B * a)
        if a & 1:
            T = -T
        return P, Q, T
    m = (a + b) >> 1
    Pl, Ql, Tl = _bs(a, m)
    Pr, Qr, Tr = _bs(m, b)
    return Pl * Pr, Ql * Qr, Tl * Qr + Pl * Tr


def _pi_decimal(stellen: int) -> str:
    n = math.ceil(stellen / 14.18) + 10
    _, Q, T = _bs(0, n)
    getcontext().prec = stellen + 20
    pi = Decimal(426880) * Decimal(10005).sqrt() * Decimal(Q) / Decimal(T)
    getcontext().prec = stellen + 1
    return str(+pi)


# ── mpmath-Backend ─────────────────────────────────────────────────────────────

def _pi_mpmath(stellen: int) -> str:
    from mpmath import mp
    mp.dps = stellen + 5
    return mp.nstr(mp.pi, stellen + 1)


# ── Backend-Auswahl ────────────────────────────────────────────────────────────

def _backend() -> tuple:
    """Gibt (name, funktion, gmpy2_verfuegbar) zurück."""
    try:
        import mpmath  # noqa: F401
        has_gmpy2 = False
        try:
            import gmpy2  # noqa: F401
            has_gmpy2 = True
        except ImportError:
            pass
        label = "mpmath + gmpy2" if has_gmpy2 else "mpmath"
        return label, _pi_mpmath, has_gmpy2
    except ImportError:
        return "decimal (binary splitting)", _pi_decimal, False


# ── CLI ────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Berechnet Pi auf beliebig viele Dezimalstellen"
    )
    parser.add_argument(
        "stellen", nargs="?", type=int, default=100,
        help="Anzahl der Dezimalstellen (Standard: 100)"
    )
    parser.add_argument(
        "-b", "--benchmark", action="store_true",
        help="Zeigt Laufzeit und Backend an"
    )
    args = parser.parse_args()

    if args.stellen < 1:
        print("Fehler: Mindestens 1 Stelle erforderlich.", file=sys.stderr)
        sys.exit(1)

    backend_name, fn, has_gmpy2 = _backend()

    start = time.perf_counter()
    pi_str = fn(args.stellen)
    dauer = time.perf_counter() - start

    print(pi_str[:args.stellen + 2])

    if args.benchmark:
        print(f"\n--- Benchmark ---", file=sys.stderr)
        print(f"Stellen:  {args.stellen}", file=sys.stderr)
        print(f"Backend:  {backend_name}", file=sys.stderr)
        print(f"Laufzeit: {dauer:.4f}s", file=sys.stderr)
        if not has_gmpy2:
            print("Tipp:     pip install gmpy2  →  GMP-Arithmetik, deutlich schneller", file=sys.stderr)


if __name__ == "__main__":
    main()
