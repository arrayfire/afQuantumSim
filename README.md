ArrayFire Quantum Simulator
=============

ArrayFire Quantum Simulator AQS is a C++14 Quantum Computer Simulator library using ArrayFire
as the backend for high performance simulators across various devices (CPUs and GPUs).

The library supports simulations of up to 30 qubits, but it's limited to the hardware's memory resources.

## Features
- Fast and optimized simulations of circuits
- Lightweight, portable, and extendable library 
- Intuitive to use with its high level of abstraction
- Allows granular control through its api allowing for easy access of the underlying arrays through ArrayFire
- Multiple special gates that implement complex features without the need of manual implementations
- Includes a circuit displayer from constructed circuits

## Example

A simple example to simulate entaglement of two qubits:
```c++
    // Initialize library
    aqs::initialize(argc, argv);

    // Create a 2-qubit Quantum Circuit
    aqs::QCircuit qc{ 2 };

    // Add gates to the circuit
    qc << aqs::H{0} << aqs::CX{ 0 , 1 };

    // Compile the circuit (Optional)
    qc.compile();

    // Create a 2-qubit Simulator with qubits initialized to the |1> state
    aqs::QSimulator qs{ 2 , aqs::QState::one() };

    // Simulate the circuit with the simulator
    qs.simulate(qc);

    // Profile the simulation for 100 simulations
    aqs::print_profile(qs.profile_measure_all(100));

```

You can display your circuits from a schematic:
```c++
    std::string schematic =
    // Declare number of qubits
    "2;"
    // Initialize their states
    "0,1;"
    "1,1;"
    // Add Gates
    "H,0,1: 0;"
    "X,1,1: 0 , 1;";

    // Print
    std::cout << aqs::gen_circuit_text_image(schematic);
```
and get an output in your terminal like this:
```sh
     ┌───┐           
|1⟩──┤ H ├──────█────
     └───┘      │    
              ┌─┴─┐  
|1⟩───────────┤ X ├──
              └───┘  
```


## Build Steps
1. Clone the library
```sh
git clone https://github.com/edwinsolisf/afQuantumSim.git my_path
```
2. Compile the library
    - CMake:
    ```sh
    cd my_path
    mkdir build && cd build
    cmake ..
    make
    ```
    - Manually: compile the source files `quantum.h`, `quantum_algo.cpp`, `quantum_gates.cpp`, `quantum_visual.cpp`, and `utils.cpp`
    with its respective headers into a library

3. Add the correct include path and link the library into your program

    - CMake: link with the target `afquantum`

You can compile the examples, tests, and benchmarks in CMake by enabling the options `AQS_BUILD_EXAMPLES`, `AQS_BUILD_TESTS`, and `AQS_BUILD_BENCHMARKS`, respectively.

## Requirements

### [ArrayFire](https://github.com/arrayfire/arrayfire)
It requires at least C++14 and ArrayFire 3.8 is sufficient to run.
It has been tested to run correctly with ArrayFire CPU and ArrayFire OpenCL.

### [NLOPT](https://github.com/stevengj/nlopt)
It requires the lastest version of NLOPT optimization library for the Variational Quantum Eigensolvers optimizers

