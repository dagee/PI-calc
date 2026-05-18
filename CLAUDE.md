# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Single-file Python script that computes Pi to an arbitrary number of decimal places using the Chudnovsky algorithm with Binary Splitting.

## Running

```bash
# Default: 100 decimal places
./piberechung.py

# Custom precision
./piberechung.py 1000

# With benchmark output (elapsed time + backend printed to stderr)
./piberechung.py 1000 --benchmark
```

## Dependencies & Performance

Backends are selected automatically in priority order:

| Backend | Install | Speedup vs. baseline |
|---|---|---|
| mpmath + gmpy2 | `pip install mpmath gmpy2` | fastest (GMP) |
| mpmath only | `pip install mpmath` | ~2000–100000× |
| decimal (fallback) | built-in | ~10× vs. old iterative |

`--benchmark` prints the active backend and a `pip install gmpy2` hint when gmpy2 is missing.

## Algorithm

**Binary Splitting** converts the Chudnovsky series into a single large-integer computation via divide-and-conquer over term indices, producing integers `(P, Q, T)` such that `pi = 426880 * sqrt(10005) * Q / T`. Only one high-precision `sqrt` is needed at the very end. This avoids per-term high-precision floating-point arithmetic and is asymptotically much faster than the iterative approach.

mpmath uses the same algorithm internally in C; the pure-Python `decimal` fallback in `_bs()` / `_pi_decimal()` applies the same Binary Splitting directly.
