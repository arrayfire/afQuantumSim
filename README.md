ArrayFire Quantum Simulator
=============

Implementation of a Quantum Computer Simulation using ArrayFire as Backend

Currently supports the simulation of up to 30 qubits (theoretically), but has only been tested up to 13 qubits.

## Features
...
## Example

A simple example to simulate entaglement of two qubits:
```c++
    // Initialize library
    aqs::initialize(argc, argv);

    // Create a 2-qubit Quantum Circuit
    aqs::QCircuit qc{ 2 };

    // Add gates to the circuit
    qc << aqs::H{0} << aqs::CX{ 0 , 1 };

    // Compile the circuit
    qc.compile();

    // Create a 2-qubit Simulator with qubits initialize to the |1> state
    aqs::QSimulator qs{ 2 , aqs::QState::one() };

    // Simulate the circuit with the simulator
    qs.simulate(qc);

    // Profile the simulation for 100 simulations
    aqs::print_profile(qs.profile_measure_all(100));

```

## Build Steps
...
## Requirements

### [ArrayFire](https://github.com/arrayfire/arrayfire)
It requires at least C++14 and ArrayFire 3.8 is sufficient to run.
It has been tested to run correctly with ArrayFire CPU and ArrayFire OpenCL.

### [NLOPT](https://github.com/stevengj/nlopt)
It requires the lastest version of NLOPT optimization library for the Variational Quantum Eigensolvers optimizers

