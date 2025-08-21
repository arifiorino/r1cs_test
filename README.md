# R1CS Matrix Multiplication

This repository proves that AB=C for NxN matrices A,B,C in zero knowledge.

This is done via R1CS constraints, and is proven via the `ark_marlin` crate.

The proof is single-threaded.

Do `cargo run 5` to run a proof with N=5.

Do `cargo run 64` to run a proof with N=64.
