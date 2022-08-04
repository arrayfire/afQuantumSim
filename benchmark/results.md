Tests
========

The QFT and Grover algorithms were run with 10 qubits and Shor with 12 qubits
All tests average measurements were done with 100 consecutive measurements of executing the algorithm.
The first number lists the average time in milliseconds, and the second number lists the sample standard deviation in milliseconds for the 100 measurements.


The experiments were run on a Macbook Pro with MacOSX12.5 Monterey, 16GB RAM, Intel Core i9-9980H, AMD Radeon Pro 5500M

## [ArrayFire Quantum Simulator](https://github.com/edwinsolisf/afQuantumSim) - C++14, [ArrayFire 3.8.2](https://github.com/arrayfire/arrayfire)
Compiled with Clang 13.1.6 x86_64-apple-darwin21.5.0

- QFT: (11.29 ± 4.18) ms
- Grover: (247 ± 43) ms
- Shor: (3927 ± 13) ms

## [Qiskit 0.20.2](https://github.com/Qiskit/qiskit) - Python 3.8.6

- QFT: (10.68 ± 4.59) ms 
- Grover: (26.40 ± 0.04) ms
- Shor: (1776 ± 1132) ms

## [Quantum++ 3.1](https://github.com/softwareQinc/qpp) - C++17
Compiled with Clang 13.1.6 x86_64-apple-darwin21.5.0

- QFT: (12.78 ± 2.10) ms
- Grover: (1431 ± 56) ms
- Shor: (2149 ± 1184) ms

|Tests\\Library|afQuantum|Qiskit|Quantum++|
|---|---|---|---|
|QFT| (11.29 ± 4.18) ms | (10.68 ± 4.59) ms | (12.78 ± 2.10) ms|
|Grover| (247 ± 43) ms| (26.40 ± 0.04) ms | (1431 ± 56) ms |
|Shor| (3927 ± 13) ms | (1776 ± 1132) ms | (2149 ± 1184) ms |

## Results
### QFT
For the Quantum Fourier Transform, afQuantumSim is 5.5% slower than Qiskit and 13.2% faster than Quantum++.

### Grover Search
For the Grover Search Algorithm, afQuantumSim is 9.4x slower than Qiskit and around 5.78x faster than Quantum++.

### Shor Prime Factorization Algorithm
For the Shor algorithm, afQuantumSim is 55% slower and 45% slower than Qiskit and Quantum++ respectively.