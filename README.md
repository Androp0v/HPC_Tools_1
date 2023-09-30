# Gauss-Jordan Elimination
This repo implements Gauss-Jordan elimination in C. Loosely based in _Numerical Recipes in C_, Chapter 2.1.

## Observations
On a M1 Pro (8P cores + 2E cores), for a 2048x2048 matrix, the LAPACK version executes in ~180ms while this version takes ~4700ms. 

The reason for this performance difference is obvious after looking at a performance trace of both the LAPACK version and the custom one:
- LAPACK uses a single CPU core for a brief period, and then launches a 10 thread computation. This is obviously because LAPACK's `dgsev` is not based on Gauss-Jordan, but rather on LU Decomposition. So the decomposition is computed first (on a single thread), and then a solver is run for every row of matrix B (concurrently, on all threads).
- Our version, based on the Gauss-Jordan elimination, can't be multithreaded (mostly, there's a small loop section that could be made multicore, shaving off ~25% of the execution time).

Not only is our version slower, it's also less numerically accurate (since it's performing a lot of sequential operations, slowly adding precision errors over time). Because of this, `check_result` may return `Result is wrong!` unless a bigger `epsilon` is chosen (`1e-3` seems to work fine for 4096x406 matrices).

## Tools
Compiled using LLVM + Clang 14.0.3 (Apple's version, bundled with Xcode). Because of that, code is annotated with `#pragma clang loop...` instead of `#pragma omp simd...` (Apple's Clang does not have OpenMP support). Not that it matters: decompiling the binary shows that the code uses ARM's vector extensions even if you don't include the compiler hints.

If an Apple platform is detected (`#if defined(__APPLE__)` for C source files, `ifeq ($(OS), Darwin)` for the makefile), the code will use Apple's Grand Central Dispatch (`<dispatch/dispatch.h>`) for a slight multicore boost in one of the loops.

The repo contains a simple Xcode project configured with an external build system preset in order to ensure the makefile is used for compilation.

# Bibliography
1. _Numerical Recipes in C_, Chapter 2.1.
2. _CS 267 Dense Linear Algebra: Parallel Gaussian Elimination_, University of California, Berkeley [link](https://people.eecs.berkeley.edu/~demmel/cs267_Spr14/Lectures/lecture13_densela_2_jwd14_4pp.pdf).
