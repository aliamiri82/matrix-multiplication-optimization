# Validation

Checked on 2026-10-23.

- All five C sources passed GNU C11 syntax checks with GCC/MinGW 8.1 and AVX/FMA/OpenMP enabled.
- Performance benchmarks and arithmetic comparisons were not executed. Compiler syntax checks do not establish numerical equivalence, valid matrix-size constraints, or speedup.
- Sources use GNU extensions and some platform-specific allocation; review platform assumptions before building on another system.
