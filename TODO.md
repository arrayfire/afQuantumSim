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

* Restructure the implementation of the visualization of gates

* Implement a rewire circuit feature

* Implement the adjoint operator

### WIP

* Implement elaborated examples (see [link](https://qiskit.org/textbook/ch-applications/algs_for_apps_index.html))

### High Priority

* Implement a Quantum Noise Simulator

* Add a QMeasurement class to deal with profiling and measurements

### Medium Priority

* Implement unitary gate

* Allow for caching of the generated fundamental gates

* Restructure the code to support more qubits and test

* Implement an optimizer (removal of gates that cancel each other, combine the addition of the same gate multiple times in a pre-generated matrix)

### Low Priority

* Implement a way to regenerate the state of individual qubits from the global state

### Scrapped

* Implement the addition of same gate multiple times with ranges to decrease computation

## Fix Bugs
* To report...