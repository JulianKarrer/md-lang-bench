# Minimal Rust vs. C++ Benchmark for Molecular Dynamics

Benchmarks the averge runtime across 100 subsequent time steps of `N` atoms on a regular lattice with Lennard-Jones potentials implemented through direct summation.

Mark `run_benchmarks.sh` as executable and run it to generate `runtimes.png`.
Requires `rustc, python3, meson`

![Runtime Plot](https://raw.githubusercontent.com/JulianKarrer/md-lang-bench/master/runtimes.png)
