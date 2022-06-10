ArrayFire Quantum Simulator
=============

## Description
Implementation of a Quantum Computer Simulation using ArrayFire as Backend

Currently supports the simulation of up to 30 qubits (theoretically), but has only been tested up to 13 qubits.

It has been tested to run correctly with ArrayFire CPU and ArrayFire OpenCL.

It requires at least C++14 and ArrayFire 3.8 is sufficient to run.

## To Do

### High Priority
* Allow for caching of the generated fundamental gates

* Optimize the generation of the fundamental `X`, `Y`, `Z`, `Hadamard`, `Phase` gates.

* Optimize the addition of custom gates by use of `CircuitGate` and `ControlCircuitGate`

* Implement the addition of same gate multiple times with ranges to decrease computation

* Allow gates to add control qubits at any location

* Implement elaborated examples (see [link](https://qiskit.org/textbook/ch-applications/algs_for_apps_index.html))

* Optimize the profiling of multiple measurements using `profile_measure` and `profile_measure_all`

* Implement the adjoint operator

* Add a QMeasurement class to deal with profiling and measurements

### Medium Priority
* Restructure the code to support more qubits and test

* Implement deferral of simulation (currently the addition of matrix requires a matrix multiplication on the spot)

* Implement an optimizer that using the deferral of simulation finds optimizations of the circuit created (removal of gates that cancel each other, combine the addition of the same gate multiple times in a pre-generated matrix)

* Implement tests for new features

### Low Priority
* Implement rotation operations: `ROTx`, `ROTy`, `ROTz`

* Implement a way to regenerate the state of individual qubits from the global state

* Restructure the implementation of the visualization of gates

## Bugs
* To report...