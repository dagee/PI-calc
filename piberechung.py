#!/usr/bin/env python3
"""
pi_berechnung.py – Pi via Chudnovsky-Algorithmus mit Binary Splitting.

Backends (Priorität):
  1. mpmath + gmpy2  – schnellstes (GMP-Arithmetik, C-Implementierung)
  2. mpmath          – schnell dank internem Binary Splitting in C
  3. decimal         – reines Python, Binary Splitting statt Iteration

Parallelisierung (--workers N):
  Verteilt Binary Splitting auf N Prozesse; --workers 0 = alle CPU-Kerne.
  Empfohlen ab ~10 Mio. Stellen. Der Speedup ist nicht linear: die finalen
  Multiplikationen von Milliarden-Digit-Zahlen bleiben sequenziell.

  pi = 426880 * sqrt(10005) * Q / T
  (P, Q, T) via rekursivem Teile-und-Herrsche über Termindizes

Installieren für maximale Geschwindigkeit:
  pip install mpmath gmpy2
"""

import argparse
import math
import multiprocessing
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from decimal import Decimal, getcontext

# ── Integer-Backend ────────────────────────────────────────────────────────────
try:
    from gmpy2 import mpz as _int
    _HAS_GMPY2 = True
except ImportError:
    _int = int
    _HAS_GMPY2 = False

# ── Konstanten ─────────────────────────────────────────────────────────────────
_A = _int(13591409)
_B = _int(545140134)
_C3_OVER_24 = _int(640320 ** 3 // 24)  # 10939058860032000


# ── Binary Splitting ───────────────────────────────────────────────────────────

def _bs(a: int, b: int) -> tuple:
    """Liefert (P, Q, T) für Chudnovsky-Terme [a, b)."""
    if b - a == 1:
        ia = _int(a)
        if a == 0:
            P, Q = _int(1), _int(1)
        else:
            P = (6*ia - 5) * (2*ia - 1) * (6*ia - 1)
            Q = _C3_OVER_24 * ia * ia * ia
        T = P * (_A + _B * ia)
        if a & 1:
            T = -T
        return P, Q, T
    m = (a + b) >> 1
    Pl, Ql, Tl = _bs(a, m)
    Pr, Qr, Tr = _bs(m, b)
    return Pl * Pr, Ql * Qr, Tl * Qr + Pl * Tr


def _bs_range(args: tuple) -> tuple:
    """Worker: Binary Splitting für Teilbereich [a, b)."""
    a, b = args
    return _bs(a, b)


def _merge(left: tuple, right: tuple) -> tuple:
    Pl, Ql, Tl = left
    Pr, Qr, Tr = right
    return Pl * Pr, Ql * Qr, Tl * Qr + Pl * Tr


def _bs_multi(n: int, workers: int) -> tuple:
    """Verteilt n Terme auf `workers` Prozesse und merged die Ergebnisse."""
    chunk = max(1, (n + workers - 1) // workers)
    ranges = [(i, min(i + chunk, n)) for i in range(0, n, chunk)]

    with ProcessPoolExecutor(max_workers=len(ranges)) as executor:
        parts = list(executor.map(_bs_range, ranges))

    # Sequenzieller Tree-Merge; bei großen n dominieren diese Multiplikationen
    while len(parts) > 1:
        parts = [
            _merge(parts[i], parts[i + 1]) if i + 1 < len(parts) else parts[i]
            for i in range(0, len(parts), 2)
        ]

    return parts[0]


# ── Pi-Berechnung ──────────────────────────────────────────────────────────────

def _pi_mpmath(stellen: int) -> str:
    from mpmath import mp
    mp.dps = stellen + 5
    return mp.nstr(mp.pi, stellen + 1)


def _pi_decimal(stellen: int) -> str:
    n = math.ceil(stellen / 14.18) + 10
    _, Q, T = _bs(0, n)
    getcontext().prec = stellen + 20
    pi = Decimal(426880) * Decimal(10005).sqrt() * Decimal(int(Q)) / Decimal(int(T))
    getcontext().prec = stellen + 1
    return str(+pi)


def _pi_parallel(stellen: int, workers: int) -> str:
    n = math.ceil(stellen / 14.18) + 10
    _, Q, T = _bs_multi(n, workers)
    # Finales sqrt via mpmath/GMP wenn verfügbar, sonst Decimal
    try:
        from mpmath import mp
        mp.dps = stellen + 5
        pi = 426880 * mp.sqrt(10005) * mp.mpf(int(Q)) / mp.mpf(int(T))
        return mp.nstr(pi, stellen + 1)
    except ImportError:
        getcontext().prec = stellen + 20
        pi = Decimal(426880) * Decimal(10005).sqrt() * Decimal(int(Q)) / Decimal(int(T))
        getcontext().prec = stellen + 1
        return str(+pi)


# ── Backend-Auswahl ────────────────────────────────────────────────────────────

def _backend(workers: int) -> tuple:
    """Gibt (name, funktion, gmpy2_verfuegbar) zurück."""
    if workers > 1:
        gmpy2_tag = "gmpy2=ja" if _HAS_GMPY2 else "gmpy2=nein"
        label = f"parallel binary splitting ({workers} Kerne, {gmpy2_tag})"
        return label, lambda stellen: _pi_parallel(stellen, workers), _HAS_GMPY2

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
    parser.add_argument(
        "-w", "--workers", type=int, default=1, metavar="N",
        help="Parallele Worker-Prozesse (0 = alle CPUs, Standard: 1 = sequenziell)"
    )
    args = parser.parse_args()

    if args.stellen < 1:
        print("Fehler: Mindestens 1 Stelle erforderlich.", file=sys.stderr)
        sys.exit(1)

    workers = args.workers if args.workers > 0 else multiprocessing.cpu_count()
    backend_name, fn, has_gmpy2 = _backend(workers)

    start = time.perf_counter()
    pi_str = fn(args.stellen)
    dauer = time.perf_counter() - start

    print(pi_str[:args.stellen + 2])

    if args.benchmark:
        print(f"\n--- Benchmark ---", file=sys.stderr)
        print(f"Stellen:  {args.stellen}", file=sys.stderr)
        print(f"Backend:  {backend_name}", file=sys.stderr)
        print(f"Laufzeit: {dauer:.4f}s", file=sys.stderr)
        if workers == 1 and not has_gmpy2:
            print("Tipp:     pip install gmpy2  →  GMP-Arithmetik, deutlich schneller", file=sys.stderr)
        if workers == 1 and args.stellen >= 10_000_000:
            print(f"Tipp:     --workers 0  →  alle {multiprocessing.cpu_count()} Kerne nutzen", file=sys.stderr)


if __name__ == "__main__":
    main()
