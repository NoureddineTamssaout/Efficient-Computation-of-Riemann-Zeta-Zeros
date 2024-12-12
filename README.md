# An Optimization Approach to the Riemann Hypothesis

## Overview

This project implements an optimized numerical approach to compute the zeros of the Riemann Zeta function using high-performance computing (HPC) techniques, this implementation is designed for accuracy, scalability, and efficiency.

## Features

- **Efficient Computation**: Optimized with Horner's method for polynomial evaluation and precomputed constants.
- **Parallelization**: Utilizes OpenMP for multi-threaded execution.
- **Precision**: Employs the Riemann-Siegel formula for accurate evaluation.
- **Scalability**: Handles large computational ranges with high resolution.

## Requirements

### Compilers
- GCC with OpenMP support.
- Intel `icpx` Compiler (recommended for best performance).


### Libraries
- OpenMP
- Standard C++17 or higher

## Usage

### Compilation

To compile with GCC:  
`g++ -O3 -fopenmp CodeRiemannOptimise.cpp -o CodeRiemannOptimise`

To compile with Intel `icpx`:  
`icpx -qopenmp -O3 CodeRiemannOptimise.cpp -o CodeRiemannOptimise_icpx`


### Execution

Run the compiled program with:  
`./CodeRiemannOptimise <START> <END> <SAMPLING>`

- **START**: Lower bound of the range.  
- **END**: Upper bound of the range.  
- **SAMPLING**: Sampling rate (resolution).
