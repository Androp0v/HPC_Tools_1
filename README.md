# Gauss-Jordan Elimination
This repo implements Gauss-Jordan elimination in C. Loosely based in _Numerical Recipes in C_, Chapter 2.1, although this version uses the (much faster) partial pivoting discussed in 2.1.2 instead of the full-pivoting of the algorithm in the book.

## Observations
On a M1 Pro (8P cores + 2E cores), for a 2048x2048 matrix, the LAPACK version executes in ~170ms while this version takes ~1600ms. 

Due to numerical inacurracies inherent to the algorithms, for big matrices `check_result` may return `Result is wrong!` unless a bigger `epsilon` is chosen (`1e-3` seems to work fine for 4096x406 matrices).

## Tools
Compiled using LLVM + Clang 14.0.3 (Apple's version, bundled with Xcode). Because of that, code is annotated with `#pragma clang loop...` instead of `#pragma omp simd...` (Apple's Clang does not have OpenMP support). Not that it matters: decompiling the binary shows that the code uses ARM's vector extensions even if you don't include the compiler hints.

If an Apple platform is detected (`#if defined(__APPLE__)` for C source files, `ifeq ($(OS), Darwin)` for the makefile), the code will use Apple's Grand Central Dispatch (`<dispatch/dispatch.h>`) for a slight multicore boost in one of the loops.

The repo contains a simple Xcode project configured with an external build system preset in order to ensure the makefile is used for compilation.

# Bibliography
1. _Numerical Recipes in C_, Chapter 2.1.
2. _CS 267 Dense Linear Algebra: Parallel Gaussian Elimination_, University of California, Berkeley [link](https://people.eecs.berkeley.edu/~demmel/cs267_Spr14/Lectures/lecture13_densela_2_jwd14_4pp.pdf).
