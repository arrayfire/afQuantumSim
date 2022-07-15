ArrayFire Quantum Simulator
=============

Implementation of a Quantum Computer Simulation using ArrayFire as Backend

Currently supports the simulation of up to 30 qubits (theoretically), but has only been tested up to 13 qubits.

## Requirements

### [ArrayFire](https://github.com/arrayfire/arrayfire)
It requires at least C++14 and ArrayFire 3.8 is sufficient to run.
It has been tested to run correctly with ArrayFire CPU and ArrayFire OpenCL.

### [NLOPT](https://github.com/stevengj/nlopt)
It requires the lastest version of NLOPT optimization library for the Variational Quantum Eigensolvers optimizers

## To Do

### Completed
* Optimize the generation of the fundamental `X`, `Y`, `Z`, `Hadamard`, `Phase` gates.

* Optimize the addition of custom gates by use of `Gate` and `ControlGate`

* Optimize the profiling of multiple measurements using `profile_measure` and `profile_measure_all`

* Implement deferral of simulation (evaluation of the circuit matrix is controlled by user)

* Implement rotation operations: `RotX`, `RotY`, `RotZ` and its controlled equivalents

* Allow gates to add control qubits at any location

* Implement tests for new features

* Project restructuring

### WIP

* Implement elaborated examples (see [link](https://qiskit.org/textbook/ch-applications/algs_for_apps_index.html))

* Implement a rewire circuit feature

### High Priority

* Implement the adjoint operator and unitary gate

* Add a QMeasurement class to deal with profiling and measurements

### Medium Priority
* Allow for caching of the generated fundamental gates

* Restructure the code to support more qubits and test

* Implement an optimizer (removal of gates that cancel each other, combine the addition of the same gate multiple times in a pre-generated matrix)

### Low Priority
* Implement a way to regenerate the state of individual qubits from the global state

* Restructure the implementation of the visualization of gates

### Scrapped

* Implement the addition of same gate multiple times with ranges to decrease computation

## Bugs
* To report...