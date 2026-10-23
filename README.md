# Matrix Multiplication Optimization

C experiments with matrix layout, loop transformations, cache blocking, GNU vector extensions, AVX/FMA, and OpenMP.

Coursework by **Ali Amiri**, Sharif University of Technology. Semester 6 — Microprocessor / Going Faster.

## Experiments

- `01-baseline/`: indexed and pointer-based multiplication.
- `02-vector-layout/`: vector layout experiment.
- `03-loop-and-layout/`: loop and memory-layout variants.
- `04-avx/`: AVX implementations.
- `05-openmp/`: parallel/vector implementation.

Build an individual source, for example:

```bash
gcc -std=gnu11 -O2 03-loop-and-layout/main.c -lm -o matrix
gcc -std=gnu11 -O2 -mavx -mfma -fopenmp 05-openmp/main.c -lm -o matrix-omp
```

Compiler extensions and platform-specific allocation appear in these coursework programs. AVX/FMA examples require compatible hardware. Benchmark defaults may allocate large matrices. No new speedup or correctness claims are made from compilation checks.

## Scope and attribution

This repository is a curated coursework archive. Original exercise structure and existing credits are retained; course-provided interfaces and testbenches may be present. Bundled textbooks, assignment PDFs, vendor IP implementations, binaries, and IDE caches are excluded. No blanket open-source license has been added to material with mixed course provenance.

See `VALIDATION.md` for checks and known limits.
