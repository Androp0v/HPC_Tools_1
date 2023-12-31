# Gaussian Elimination
This repo implements Gaussian elimination in C. Originally based in _Numerical Recipes in C_, Chapter 2.1, with several changes:
- This version uses the (much faster) partial pivoting, discussed in Chapter 2.1.2 [1], instead of the full-pivoting of the algorithm in the book.
- This version implements a _delayed updates_ mechanism, inspired by the suggestions in [2].

# Results
On a MacBook Pro with an M1 Pro CPU (8P cores + 2E cores):
| Matrix size | LAPACK (ms) | Ours (ms) | Ours slower by |
| ----------- | ----------- | --------- | -------------- |
| 512x512     | 11          | 30        | 2.7x           |
| 1024x1024   | 35          | 129       | 3.6x           |
| 2048x2048   | 175         | 800       | 4.6x           |
| 4096x4096   | 1600        | 6500      | 4.1x           |

## Profiling
See below a trace captured using Xcode Instruments, which has been extremely useful in both diagnosing the bottlenecks in the algorithm's execution and estimating the CPU utilization.

![Instruments-Trace](Figures/Xcode-Instruments-Trace.png)

## Tools
Compiled using LLVM + Clang 14.0.3 (Apple's version, bundled with Xcode). Because of that, code is annotated with `#pragma clang loop...` instead of `#pragma omp simd...` (Apple's Clang does not have OpenMP support). Not that it matters: decompiling the binary shows that the code uses ARM's vector extensions even if you don't include the compiler hints.

If an Apple platform is detected (`#if defined(__APPLE__)` for C source files, `ifeq ($(OS), Darwin)` for the makefile), the code will use Apple's Grand Central Dispatch (`<dispatch/dispatch.h>`) to parallelize the loops. This is because Apple Silicon uses a heterogeneous CPU architecture, which complicates load balancing, and Grand Central Dispatch does the load balancing for us (also, this bypasses the need for `#pragma omp parallel for` which is not supported using Apple's Clang anyway).

For profiling, singposting was introduced using Apple's [`OSSignposter`](https://developer.apple.com/documentation/os/ossignposter). This can be disabled by removing `#define USE_SIGNPOSTING` (this is only defined by default for Apple platforms).

The repo contains a simple Xcode project configured with an external build system preset in order to ensure the makefile is used for compilation.

# Bibliography
1. _Numerical Recipes in C_, Chapter 2.1.
2. _CS 267 Dense Linear Algebra: Parallel Gaussian Elimination_, University of California, Berkeley [link](https://people.eecs.berkeley.edu/~demmel/cs267_Spr14/Lectures/lecture13_densela_2_jwd14_4pp.pdf).
