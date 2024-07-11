#!/bin/bash

# compile and execute C++ benchmark
cd ./cpp
meson setup build --buildtype=release
cd ./build
meson compile
cd ..
./build/bench


# compile and execute Rust benchmark
cd ../rust
RUSTFLAGS='-C target-cpu=native' cargo run -r
cd ..
python3 plot_results.py