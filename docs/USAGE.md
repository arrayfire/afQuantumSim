Usage
======

ArrayFire Quantum Simulator AQS is a C++14 Quantum Computer Simulator library using ArrayFire
as the backend for high-performance simulators across various devices.

The main goal of this library is to provide an easy workflow for designing and simulating quantum circuits through a high level API,
but also allowing low level access and fast, high performance computations through the use of ArrayFire.

## Table of Contents
1. [Initialization](#initialization)
2. [Quantum States](#quantum-states)
    1. [Construction](#construction)
    2. [Gates](#gates)
    2. [Operations](#operations)
3. [Quantum Circuit](#quantum-circuit)
    1. [Construction](#construction-1)
    3. [Addition of Gates](#addition-of-gates)
    4. [Compilation and Execution](#compilation-and-execution)
4. [Quantum Simulator](#quantum-simulator)
    1. [Construction](#construction-2)
    2. [Modifying Initial States](#modifying-initial-states)
    3. [Simulation](#simulation)
    4. [Results](#results)
5. [Workflow](#workflow)
6. [Quantum Algorithms](#quantum-algorithms)
    1. [Quantum Fourier Transform](#quantum-fourier-transform)
    2. [Inverse Quantum Fourier Transform](#inverse-quantum-fourier-transform)
    3. [Grover Search](#grover-search)
    4. [Hamiltonian Decomposition](#hamiltonian-decomposition)
    5. [Hamiltonian Evolution](#hamiltonian-evolution)
    6. [Variational Quantum Eigensolver](#variational-quantum-eigensolver)
7. [Special Gates](#special-gates)
8. [Visuals](#visuals)
    1. [State](#state)
    2. [Profile Results](#profile-results)
    3. [Circuit Matrix](#circuit-matrix)
    4. [Circuit Image](#circuit-image)
9. [Comprehensive Example](#comprehensive-example)


## Initialization

All the functions, classes, and operators are declared inside the namespace `aqs`.

To use the library, it must be first initialized with the function `initialize`.
It must receive the command line arguments, and optionally it can receive the specific ArrayFire backend to use:

```c++
    int main(int argc, char** argv) {
        ...
        aqs::initialize(argc, argv /*, Optional: af::Backend::AF_BACKEND_CUDA*/);
        ...
    }
```

It uses the default backend set by ArrayFire and the first device from that backend if not specified in the command line arguments.

## Quantum States

In this library, `QState` represents the quantum state of a qubit expressed with $|0\rangle$ and $|1\rangle$ in the computational basis.
It stores the normalized components of the basis vectors. It is located in the [`quantum.h`](../include/quantum.h) file.

### Construction
To create the `QState` object, you must pass the complex coefficients of each component:
```c++
    aqs::QState zero_state{ 1.f , 0.f };
    aqs::QState one_state{ 0.f , 1.f };
    aqs::QState some_state{ { 0.5f, 0.5f } , { 0.5f , 0.5f } };
```

For the common $|0\rangle$ and $|1\rangle$ states, you can use `QState::zero()` and `QState::one()` for those states respectively.

Additionally, for creating states with certain angles on the bloch sphere, there is the function `QState::create_state` which creates
an array with the complex coefficients of the described state:

```c++
    float azimut = aqs::pi / 2.f;
    float polar  = aqs::pi;

    aqs::QState state{ aqs::QState::create_state(azimut, polar) };
    // Equivalent to aqs::QState{ { 1/sqrt(2) , 0.f } , { -1/sqrt(2) , 0.f } };
```

### Operations
For single-qubit operations there are the functions which return the applied the given operations:
- `X_op`: Returns the result of applying pauli X matrix
- `Y_op`: Returns the result of applying pauli Y matrix
- `Z_op`: Returns the result of applying pauli Z matrix
- `Hadamard_op`: Returns the result of applying the hadamard matrix
- `RotX_op`: Returns the result of applying a rotation around the X axis
- `RotY_op`: Returns the result of applying a rotation around the Y axis
- `RotZ_op`: Returns the result of applying a rotation around the Z axis
- `Phase_op`: Returns the result of applying a phase rotation on the $|1\rangle$ of the state

Example:
```c++
    aqs::QState result = aqs::Hadamard_op(aqs::QState::one());
    // Results in QState{ { 1/sqrt(2) , 0.f } , { -1/sqrt(2) , 0.f } };
```

## Quantum Circuit
While it is possible to simulate a quantum computer with qubits using individual `QState`s, it limits the usefulness
of taking advantage of the entangling multiple qubits needed for quantum algorithms.

So for simulating a quantum circuit, the library provides the `QCircuit` class which manages the creation of a circuit with
multiple quantum gates.

### Construction
To create a `QCircuit` object, you must pass the number of qubits (i.e. quantum registers) that the circuit contains:
```c++
    uint32_t qubit_count = 4;
    aqs::QCircuit qc{ qubit_count };
```

### Gates
The usefulness of a Quantum Computer comes from the gates contained in it.
Therefore, the library provides a structure for implementing gates, adding them to the circuit, and for executting
the simulation's calculations.

In this library, a quantum gate is represented by any class that implements the `QGate` abstract class.
This allows for each gate to implement its logic in a predictable way despite each gate's actions.

The most important part of the gate is the modification it does to the circuit. This is done with the function `QGate::operator()(QCircuit&)`.
In here, each gate must include its logic to modify the internal circuit matrix that represents the simulation of the quantum circuit.

For the majority (and likely all) of quantum programs to be simulated with the circuit, they will require the basic quantum gates.
This is the list of all the gates provided by the library:

- `X`: Pauli X gate
- `Y`: Pauli Y gate
- `Z`: Pauli Z gate
- `H`: Pauli Hadamard gate
- `RotX`: Pauli X rotation
- `RotY`: Pauli Y rotation
- `RotZ`: Pauli Z rotation
- `Phase`: Phase rotation of $|1\rangle$
- `Swap`: Swaps two qubit states
- `CX`: Controlled Pauli X gate
- `CY`: Controlled Pauli Y gate
- `CZ`: Controlled Pauli Z gate
- `CH`: Controlled Hadamard gate
- `CRotX`: Controlled Pauli X rotation gate
- `CRotY`: Controlled Pauli Y rotation gate
- `CRotZ`: Controlled Pauli Z rotation gate
- `CPhase`: Controlled Phase rotation
- `CSwap`: Controlled swap
- `CCX`: Double Controlled Pauli X gate, also known as the Toffoli gate
- `Or`: Or of two qubits 
- `Gate`: Allows for embedding a previously created circuit into another one
- `ControlGate`: Allows for embedding a previously created circuit into another one with a control qubit for the whole gate

### Addition of Gates
Adding gates to the circuit is done by using the `<<` operator with a `QGate` derived object. This results on inserting
the gate into a list of gates stored by the circuit.

Example:
```c++
    aqs::QCircuit qc{ 2 };

    qc << aqs::H{ 0 } // Add a Hadamard gate in qubit 0
       << aqs::CX{ 0 , 1 }; // Adds a X gate in qubit 1 controlled by qubit 0
```

### Compilation and Execution
Adding gates to the circuit just stores them in a list, but does not compute anything. To ensure the gates are added into the circuit through matrix computation
the `QCircuit` object must call the `compile` function to start the computations. This is to allow the program to defer this calculation when needed:

```c++
    aqs::QCircuit qc{ 2 };

    // ... Addition of gates

    qc.compile();
```

After this function is called, all the gates that have been added modify the `QCircuit` in the way the gate its expected to do it.

The result of this computation is stored in an internal `af::array` that stores the matrix representation of the circuit. You can retrieve this
result by calling the `circuit()` function.

Note that this operation is only required if a the matrix representation of the circuit is needed. You are able to simulate the circuit directly
without compilation.

## Quantum Simulator
While it is nice to get the representation of a circuit in matrix form, the main usefulness of a Quantum Computer is to measure the way the circuit
affects the state of the qubits, and then meaasure them. In addition, it is desirable to be able to test the same circuit with different inputs and
compare the outputs for each simulation.

For these needs, the library provides the `QSimulator` class which manages the initial statevector for the quantum computer simulation.

### Construction
`QSimulator` has many constructors depending on the needs of the program, but the first argument is always the number of qubits. This number must match
the circuit that it will simulate:

```c++
    uint32_t qubit_count = 5;

    aqs::QSimulator qs{ qubit_count }; // Sets all qubits to |0>
```

You can set the initial state of all qubits by passing a single `QState` object to the constructor. By default, the initial state for all qubits is the $|0\rangle$ state (same as `QState::zero`).

```c++
    aqs::QState initial_state{ {0.1f , 0.7f}, {0.5f , -0.5f} };

    aqs::QSimulator qs{ qubit_count , initial_state }; // Sets all qubits to (0.1 + 0.7i) |0> + (0.5 - 0.5i) |1>
```

It is also possible to set the state of each qubit individually by passing an `std::vector<QState>` with each of the states.

```c++
    std::vector<aqs::QState> states{ aqs::QState::zero() , aqs::QState::one() , aqs::QState{ 0.6f , 0.8f} };

    aqs::QSimulator qs{ 3 , states };
```
Note that the size of the vector must match the number of qubits in the simulator.

Ultimately, there is also the possibility of initializing a simulator with already created statevector.

```c++
    af::array statevector = create_statevector();

    aqs::QSimulator{ qubit_count , statevector };
```

The statevector must be an `af::array` and be a non-zero vector.

### Modifying initial states

The `QSimulator` object allows for modification of the states after construction by modifying the initial states of each qubit.
Each qubit's state can be obtained by calling `QSimulator::qubit(uint32_t qubit)` function, and with this reference, you modify the state.

After any modification in the state of the qubits, for the changes to take effect it is required that you update the statevector by
calling the member function `generate_statevector`. Once this function is called, it will overwrite the statevector with the result of the
tensor product of the initial states stored in the `QState` object of each qubit.

```c++
    aqs::QSimulator qs{ 2 };
    af::array statevector = qs.statevector(); // [1 0 0 0]
    
    qs.qubit(0) = aqs::QState::one();
    statevector = qs.statevector(); // Still [1 0 0 0]

    qs.generate_statevector();
    statevector = qs.statevector(); // Now it is [0 0 1 0]
```

### Simulation
Finally to simulate a circuit with the statevector prepared, you call the member function `simulate` and pass the circuit to simulate:

```c++
    aqs::QCircuit qc{ qubit_count };
    ...

    aqs::QSimulate qs{ qubit_count };
    ...

    qs.simulate(qc);
```

After executing this function, the computation will start and once finished, the result will be stored in the statevector of the `QSimulator` object.
This way, it is possible to chain computations of different circuits from the result of a simulation to the input of another, by just calling
the `simulate` function with the next circuit in line.

Note that when simulating, the computation will use the circuit matrix if its has been compiled; otherwise, it will compute the gates individually on the statevector.

### Results
After the simulation is done, we would like to obtain results in the form of measurements, probabilities, or even the actual resulting statevector.
For that, the `QSimulator` class contains useful methods for accomplishing these tasks:

- `statevector`: This function returns the internal statevector resulting from the circuit simulation. Useful when it is needed to see the actual states of the registers.
- `probabilities`: Returns the probabilities of all the possible measurements that can be done in the computational basis from the statevector.
- `qubit_probability_true` and `qubit_probability_false`: Returns the probability of measuring a single qubit `true` $|1\rangle$ and `false` $|0\rangle$
in the compuatational basis, respectively.
- `state_probability`: Returns the probability of measuring that state in the compuatational basis from the statevector.
- `measure`: Makes a measurement on a single qubit, collapses the qubit's state in the statevector and returns the result of the measurement.
- `measure_all`: Makes a measurement of the whole statevector, updates the statevector array to this measurement, and returns the result of the statevector measurement.
- `peek_measure`: Makes a measurement on a single qubit and returns the result of the measurement. Does not collapse (modify) the statevector.
- `peek_measure_all`: Makes a measurement of the whole statevector and returns the result of the statevector measurement. Does not collapse (modify) the statevector.
- `profile_measure`: Makes multiple measurements on a single qubit and returns the distribution of the measurement results. Does not collapse (modify) the statevector.
- `profile_measure_all`: Makes multiple measurements of the whole statevector and returns the distribution of the statevector measurements. Does not collapse (modify) the statevector.

## Workflow

The general workflow for simulating quantum circuits with this library can be described the following pseudocode:
```c++
    // Create the circuit
    aqs::QCircuit qc{ qubit_count };

    // Add the gates
    qc << Gates;

    // Compile the circuit (Optional)
    qc.compile();

    // Create the simulator
    aqs::QSimulator qs{ qubit_count, initial_states };

    // Simulate the circuit with the initial state
    qs.simulate(qc);

    // Make some measurements
    auto results = qs.profile_measure_all(counts);
```

This workflow is advantageous as it allows great flexibility for the many needs for simulations:
- Add all the gates at once and defer the computation when needed, or add gates and compute by parts
- Separate a circuit into multiple step circuits and compute them separately and then add them to a single circuit
or simulate each of them by step and get results at each step.
- Simulate the same circuit with different initial statevectors and simulate the same statevector with different circuits

All of these possibilities with little code and low computational overhead.

## Quantum Algorithms
For commonly used quantum algorithms, the library provides functions that create these algorithms in circuits.

All implemented algorithms are contained in the file [`quantum_algo.h`](../include/quantum.h).

### Quantum Fourier Transform
The Quantum Fourier Transform is a variation of the Fourier Transform which transforms the statevector from the computational basis to
the fourier basis. It is usually represented as _QFT_.
This algorithm consists in encoding the overall statevector into phase rotations in each qubit.

To utilize this algorithm, the library provides the function `fourier_transform` which returns a `QCircuit` with the Quantum Fourier Transform
algorithm encoded inside it:

```c++
    aqs::QCircuit qft_circuit = aqs::fourier_transform(qubit_count);
```

### Inverse Quantum Fourier Transform

In conjuction to the Quantum Fourier Transform, there is the _Inverse Quantum Fourier Transform_ which does the reverse operation, that is, transforming
the statevector from the fourier basis to the computational basis. It is usually represented as _QFT†_.

The library provides the algorithm with the similarly named function `inverse_fourier_transform` which returns the adjoint of the _QFT_ circuit.

```c++
    aqs::QCircuit inverse_qft_circuit = aqs::inverse_fourier_transform(qubit_count);
```

This algorithm is really useful for obtaining the computational basis representation for qubit rotations such as in Quantum Phase Estimation.

### Grover Search
A well-known quantum algorithm is the grover search algorithm which gets its fame for executing in $ O( \sqrt{n} ) $ complexity.

The library provides an implementation of the search algorithm in the function `grover_search`. Given an oracle, that is, the gate that marks the search state, the grover gate executes the
grover algorithm (oracle + amplification). The library also provides the `grover_oracle` function that returns a custom oracle for making one unique state.

For example, searching the state `0010` using the algorithm would be written as
```c++
    uint32_t qubit_count = 4;
    uint32_t marked_state = 0b0010;
    uint32_t iterations = 1;
    aqs::QCircuit oracle = aqs::grover_oracle(qubit_count, marked_state);
    aqs::QCircuit grover = aqs::grover_search(qubit_count, oracle, iterations);
```

There are other algorithms that require grover iterations instead of the whole grover search algorithm, so for that there is the `grover_iteration`
function.

Note that the grover oracle circuit marks the state by changing the sign of the $|1\rangle$ component for only that state. Similarly, the grover iteration
and grover search functions amplify only those states which were affected by the pauli-Z gate.

### Hamiltonian Decomposition
It is common in the area of variational methods of wanting to encode a matrix into a circuit. The most straight forward way is to decompose the hamiltonian matrix
as a sum of pauli gates operations on the circuit. For this case, the library provides the `decompose_hamiltonian` function which does exactly this. It returns a
vector with the pauli gates layout and the coefficient for that term.

```c++
    af::array matrix;
    ...

    std::vector<std::pair<std::string, af::cfloat>> decomposition = aqs::decompose_hamiltonian(qubit_count, matrix);
```

The string part is a description of the pauli gates. It will contain the same amount of letters as number of qubits, and the letters can be `i`,
`x`, `y`, and `z` which stand for no gate, pauli X gate, pauli Y gate, and pauli Z gate, respectively.

There is also the reverse operation, in which given a description of the pauli decomposition for the matrix, the function `compose_hamiltonian` reconstructs the matrix
and stores it into an `af::array`.

```c++
    af::array matrix = aqs::compose_hamiltonian(description);
```

### Hamiltonian Evolution
Another famous use for a quantum computer is its ability to simulate a quantum system in polynomial time.
This is another algorithm that the library provides through a circuit with the function `hamiltonian_evolution_circuit` which returns
the evolved version of the hamiltonian matrix passed for a given time step.

```c++
    uint32_t step_count;
    af::array matrix;

    ...

    aqs::QCircuit evolution_circuit = aqs::hamiltonian_evolution_circuit(matrix, step_count);
```

With the circuit it is possible to find the evolved state at any give time point as follows:
```c++
    float time;

    ...

    for (uint32_t i = 0; i < time * step_count; ++i)
        qs.simulate(qc);
```

### Variational Quantum Eigensolver
Using quantum computers has been proving to be very useful when combined with variational methods in order to solve linear algebra problems.
Among those problems, there is the common problem of finding the smallest eigenvalue of a matrix. This is mainly useful for determining the ground state
of quantum systems.

Therefore, the library provides an implemented version of using the Variational Quantum Eigensolver algorithm for these purposes. It utilizes
the local non-gradient `Cobyla` optimization algorithm implemented by the `NLopt` library in conjuction to the simulations of the matrix through the quantum circuits.

The function `variational_quantum_eigensolver` executes the previously mentioned algorithm and returns the found minimum eigenvalue and the parameters found
for the state generator circuit. It takes in the hamiltonian matrix to determine the smallest eigenvalue of, the range to search the eigenvalues in,
the type of state generator used in the variation of parameters, the precision or tolerance of the search, and the maximum number of iterations for the search.

```c++
    af::array matrix;

    ...

    std::pair<float, std::vector<float>> result = aqs::variational_quantum_eigensolver(matrix, range);

    float minimum_eigenvalue = result.first;
    const std::vector<float>& state_parameters = result.second;
```

The library provides two variational quantum state generators as circuits: `VQE::LINEAR` AND `VQE::FULL`. The difference between both is the capabilities of
entaglement in the state that is generated and the number of gates for the generator. `VQE::LINEAR` has a linear number of gates while `VQE::FULL` has
a squared number of qubits gates.

With the type of state generator and the parameters obtained from the eigensolver, it is possible to recover the eigenstate as follows:
```c++
    std::pair<float, std::vector<float>> result = aqs::variational_quantum_eigensolver(matrix, range, VQE::FULL);
    const std::vector<float>& state_parameters = result.second;

    aqs::QCircuit state_generator = aqs::full_entanglemente_varstate(qubits, qubits, state_parameters);
```

## Special Gates

To facilitate creation of complex gates like Multiple Control Phase Gate, the Adjoint of a Gate,
or adding multiple gates at the same time, the library facilitates multiple functions that create
those circuits.

All implemented special gates are contained in the file [`quantum_gates.h`](../include/quantum.h)

- `Group_Gate`: Adds a given gate at multiple specified locations
- `Control_Group_Gate`: Adds a given gate at multiple specified locations all which are controlled by a given qubit
- `NControlGate`: Adds a gate controlled by multiple qubits
- `Rewire_Gate`: Reconnects the inputs of a gate to the specified targets
- `Adjoint_Gate`: Creates the adjoint (inverse) of a given gate

## Visuals

For displaying measurements or circuits, the library provides many functions to print the information to the standard outstream (console/terminal).
This functions are located in [`quantum_visuals.h`](../include/quantum_visuals.h).

### State
In a real quantum computer it is not possible to obtain the complete quantum state of a qubit nor for the whole statevector for the quantum computer.
However, for the simulators it is possible to obtain them which may be useful for checking mathematical models or inspecting an algorithm.

That is why to visualize these quantities the library provides the functions `print_state` and `print_statevector`.

### Profile Results
In general, to obtain results from quantum computers, qubits are measured multiple times to get an idea of the probability distributions for each state.
As this process is random, the functions `profile_measure` and `profile_measure_all` give you a vector with the results of each measurement.
In order to visualize the probability distributions for the states from these measurements, the library provides the `print_profile` function which displays
the counts and overall probability for each state from the random measurements.

```c++
    aqs::QSimulator qs{ qubit_count };

    ...

    auto profile = qs.profile_measure_all(counts);

    aqs::print_profile(profile);
```


### Circuit Matrix
As quantum circuits are represented internally by the matrices that operate on the different states. For this case, one can display the matrix representation
of the circuit by using the function `print_circuit_matrix`.

```c++
    aqs::QCircuit qc{ qubit_count };

    ...

    qc.compile();
    aqs::print_circuit_matrix(qc);
```


### Circuit Image
Another useful tool for developing, debugging, and visualizing quantum circuits, its displaying the drawing of all gates and connections in a circuit.
This is why the library contains a text circuit displayer. To use it, you use the function `gen_circuit_text_image` with the created circuit and simulator:

```c++
    aqs::QCircuit qc{ qubit_count };
    aqs::QSimulator qs{ qubit_count };

    // Add gates...

    std::string image = aqs::gen_circuit_text_image(qc, qs);
```

This function returns the text image as a string of UTF-8 characters which can then be dumped into a file or to the console.

As a shortcut, to display the circuit directly to the console (standard output stream), you can use the `print_circuit_text_image` function
with the same parameters as before.

For example, the following code:
```c++
    aqs::QCircuit qc{ 2 };

    qc << aqs::H{0} << aqs::CX{ 0 , 1 };

    aqs::QSimulator qs{ 2 , aqs::QState::one() };

    aqs::print_circuit_text_image(qc, qs);
```

would produce the output:
```sh
     ┌───┐           
|1⟩──┤ H ├──────█────
     └───┘      │    
              ┌─┴─┐  
|1⟩───────────┤ X ├──
              └───┘  
```

In addition to displaying circuits, it also possible to write an schematic for the circuit image to display.
To produce the text image given a circuit schematic, you use the function `gen_circuit_text_image` and pass a `std::string` with
the circuit schematic.

For example, the process for the circuit above with an schematic would be:
```c++
    std::string schematic =
    "2;" // Declare number of qubits

    "0,1;" // Initialize their states
    "1,1;"

    "H,0,1: 0;" // Add Gates
    "X,1,1: 0 , 1;";

    std::cout << aqs::gen_circuit_text_image(schematic);
```

The general structure for the schematic to create a circuit text image is the following:

1. Declare the number of qubits in the circuit.

2. Optionally, declare the state of each of the qubits by starting with the index of the qubit, and separated by a comma, the state of the qubit.
    Note that only currently $|0\rangle$ and $|1\rangle$ are supported for initial states.

3. Add the gates in the circuit. For each gate added you must specify:
    1. Name of the gates (can be separated by spaces and lines by not by commas, colons, or semicolons)
    2. Number of control qubits controlling the gate
    3. Number of target qubit the gate is affecting
    4. A colon to separate the list of qubits
    5. Comma separated list of control qubits
    6. Comma separated list of target qubits

    For example: `QFT,2,4:6,1,2,3,4,5` represents a gate called `QFT` that is controlled by the qubits 1 and 6, and targets/affects
    qubits 2,3,4, and 5.


Note that each _statement_ is terminated by a semicolon.

## Comprehensive Example
One of the most challenging scientific tasks done in quantum physics is simulating a quantum system in a timely manner. With classical computers
due to the many body interactions and the computationally intesive job of solving the Schrödinger equation, the time complexity is usually exponential.
However, with the emergence of quantum computers, that is expected to change.

Among the goals of simulating quantum systems, there is determining the **Ground Energy State** of a system, that is, the minimum energy state of the system.
This quantity is useful because it gives us insight in the behavior of the system at certain conditions and also tells us how it evolves through time.
Knowing the ground state is very useful in Material Science as well as in Chemistry.

For example, let's take as system a hydrogen molecule $H_2$. This is a symple quantum system containing two atoms, both made up of one proton and one electron.
Effectively, there are 4 quantum particles in this system.

What we would like to find is the ground state for this system. To do this, we will use the [**Hartree-Fock Method**](https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method) to accomplish this.
In simple terms, what the algorithm does is to approximate the solution to the Schrödinger equation for the given hamiltonian of the system.
In this case, our hamiltonian is a molecular hamiltonian, which in our case is just a matrix which represents the energy of the particles in the system.
This matrix changes depending on the molecular structure, that is where the particles are located.

For our simulations, let's suppose the hydrogen atoms in the molecule are separated by $1.322  \overset{\circ}{A}$. Then after computing the molecular hamiltonian for this setup
we obtain the matrix in atomic units:

$$
    \small
    \begin{bmatrix}
         0.756 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
         0 & 0.3077 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
         0 & 0 & 0.3077 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
         0 & 0 & 0 & 0.5645 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0.1790 & 0 & 0 & 0\\
         0 & 0 & 0 & 0 & -0.5219& 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
         0 & 0 & 0 & 0 & 0 & -0.4784& 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
         0 & 0 & 0 & 0 & 0 & 0 & -0.2994& 0 & 0 & -0.1790& 0 & 0 & 0 & 0 & 0 & 0\\
         0 & 0 & 0 & 0 & 0 & 0 & 0 & 0.4491 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
         0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & -0.5219& 0 & 0 & 0 & 0 & 0 & 0 & 0\\
         0 & 0 & 0 & 0 & 0 & 0 & -0.1790& 0 & 0 & -0.2994& 0 & 0 & 0 & 0 & 0 & 0\\
         0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 &-0.4784 & 0 & 0 & 0 & 0 & 0\\
         0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0.4491 & 0 & 0 & 0 & 0\\
         0 & 0 & 0 & 0.179 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & -1.1173& 0 & 0 & 0\\
         0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & -0.4032& 0 & 0\\
         0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 &-0.4032 & 0\\
         0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1.0161
    \end{bmatrix}
$$

Quite the handful of a matrix. Using the `decompose_hamiltonian` function we can obtain the pauli decomposition for the hamiltonian, which will allow us to see how we can express the hamiltonian with
basic pauli matrices. This allows us to represent this matrix into a more compact way:

```c++
    auto hamiltonian_description = aqs::compose_hamiltonian(hamiltonian_matrix);
    /* Returns

        std::vector<std::pair<std::string, af::cfloat>> hamiltonian_description = {
        {"iizi", {-0.24274501250395486f}},
        {"iiiz", {-0.24274501250395486f}},
        {"iiii", {-0.04207255204090424f}},
        {"ziii", {0.17771358235540047f}},
        {"izii", {0.1777135823554005f}},
        {"zizi", {0.12293330446049033f}},
        {"iziz", {0.12293330446049033f}},
        {"ziiz", {0.16768338851167847f}},
        {"izzi", {0.16768338851167847f}},
        {"zzii", {0.17059759275420894f}},
        {"iizz", {0.1762766138632343}},
        {"yyxx", {-0.044750084051188126f}},
        {"xxyy", {-0.044750084051188126f}},
        {"yxxy", {0.044750084051188126f}},
        {"xyyx", {0.044750084051188126f}}
        };
    */
```

However, the actual usefulness comes from the simulation of evolution of pauli gates is straight forward to implement in
quantum circuits.
You can observe the circuit representation of the evolved matrix by using `hamiltonian_evolution_circuit`:
```c++
    aqs::QCircuit evolved_circuit = aqs::hamiltonian_evolution_circuit(hamiltonian_matrix, 1);

    aqs::print_circuit_text_image(evolved_circuit, aqs::QSimulator{4});
```

This would output:

```
     ┌──────┐    ┌───────┐    ┌──────┐                                                                                                                      ┌───┐    ┌──────┐                                                                      ┌───┐    ┌──────┐    ┌───┐                                                                                  ┌───┐                ┌───┐    ┌──────┐                                                                      ┌───┐    ┌──────┐    ┌───┐                                                                         ┌───┐                                                                                                           
|0⟩──┤ RotZ ├────┤ Phase ├────┤ RotZ ├──────────────────█────────────────────█────────────────────█────────────────────█────────────────────────────────────┤ Y ├────┤ RotX ├──────█────────────────────────────────────────────────────────█──────┤ Y ├────┤ RotX ├────┤ H ├──────────────────█────────────────────────────────────────────────────────█──────┤ H ├────────────────┤ Y ├────┤ RotX ├──────█────────────────────────────────────────────────────────█──────┤ Y ├────┤ RotX ├────┤ H ├──────█───────────────────────────────────────────────────────────█──────┤ H ├─────────────────────█────────────────────█────────────────────────────────────────────────────────────────
     └──────┘    └───────┘    └──────┘                  │                    │                    │                    │                                    └───┘    └──────┘      │                                                        │      └───┘    └──────┘    └───┘                  │                                                        │      └───┘                └───┘    └──────┘      │                                                        │      └───┘    └──────┘    └───┘      │                                                           │      └───┘                     │                    │                                                                
     ┌──────┐    ┌───────┐                ┌──────┐    ┌─┴─┐    ┌──────┐    ┌─┴─┐                  │                    │                                    ┌───┐    ┌──────┐    ┌─┴─┐                                                    ┌─┴─┐    ┌───┐    ┌──────┐    ┌───┐    ┌──────┐    ┌─┴─┐                                                    ┌─┴─┐    ┌───┐    ┌──────┐    ┌───┐                ┌─┴─┐                                                    ┌─┴─┐    ┌───┐                ┌───┐    ┌─┴─┐                                                       ┌─┴─┐    ┌───┐                     │                    │                                                                
|0⟩──┤ RotZ ├────┤ Phase ├────────────────┤ RotZ ├────┤ X ├────┤ RotZ ├────┤ X ├──────────────────┼────────────────────┼────────█────────────────────█──────┤ Y ├────┤ RotX ├────┤ X ├──────█──────────────────────────────────────█──────┤ X ├────┤ Y ├────┤ RotX ├────┤ Y ├────┤ RotX ├────┤ X ├──────█──────────────────────────────────────█──────┤ X ├────┤ Y ├────┤ RotX ├────┤ H ├────────────────┤ X ├──────█──────────────────────────────────────█──────┤ X ├────┤ H ├────────────────┤ H ├────┤ X ├─────────█──────────────────────────────────────█──────┤ X ├────┤ H ├─────────────────────┼────────────────────┼────────█────────────────────█──────────────────────────────────
     └──────┘    └───────┘                └──────┘    └───┘    └──────┘    └───┘                  │                    │        │                    │      └───┘    └──────┘    └───┘      │                                      │      └───┘    └───┘    └──────┘    └───┘    └──────┘    └───┘      │                                      │      └───┘    └───┘    └──────┘    └───┘                └───┘      │                                      │      └───┘    └───┘                └───┘    └───┘         │                                      │      └───┘    └───┘                     │                    │        │                    │                                  
     ┌──────┐    ┌───────┐                                                          ┌──────┐    ┌─┴─┐    ┌──────┐    ┌─┴─┐    ┌─┴─┐    ┌──────┐    ┌─┴─┐    ┌───┐                         ┌─┴─┐                                  ┌─┴─┐    ┌───┐                         ┌───┐    ┌──────┐             ┌─┴─┐                                  ┌─┴─┐    ┌───┐    ┌──────┐             ┌───┐                         ┌─┴─┐                                  ┌─┴─┐    ┌───┐                         ┌───┐    ┌──────┐    ┌─┴─┐                                  ┌─┴─┐    ┌───┐    ┌──────┐                  │                    │        │                    │                                  
|0⟩──┤ RotZ ├────┤ Phase ├──────────────────────────────────────────────────────────┤ RotZ ├────┤ X ├────┤ RotZ ├────┤ X ├────┤ X ├────┤ RotZ ├────┤ X ├────┤ H ├─────────────────────────┤ X ├──────█────────────────────█──────┤ X ├────┤ H ├─────────────────────────┤ Y ├────┤ RotX ├─────────────┤ X ├──────█────────────────────█──────┤ X ├────┤ Y ├────┤ RotX ├─────────────┤ H ├─────────────────────────┤ X ├──────█────────────────────█──────┤ X ├────┤ H ├─────────────────────────┤ Y ├────┤ RotX ├────┤ X ├──────█────────────────────█──────┤ X ├────┤ Y ├────┤ RotX ├──────────────────┼────────────────────┼────────┼────────────────────┼────────█────────────────────█────
     └──────┘    └───────┘                                                          └──────┘    └───┘    └──────┘    └───┘    └───┘    └──────┘    └───┘    └───┘                         └───┘      │                    │      └───┘    └───┘                         └───┘    └──────┘             └───┘      │                    │      └───┘    └───┘    └──────┘             └───┘                         └───┘      │                    │      └───┘    └───┘                         └───┘    └──────┘    └───┘      │                    │      └───┘    └───┘    └──────┘                  │                    │        │                    │        │                    │    
     ┌──────┐    ┌───────┐                                                                                                                                  ┌───┐                                  ┌─┴─┐    ┌──────┐    ┌─┴─┐    ┌───┐                                  ┌───┐                                  ┌─┴─┐    ┌──────┐    ┌─┴─┐    ┌───┐                                  ┌───┐    ┌──────┐                      ┌─┴─┐    ┌──────┐    ┌─┴─┐    ┌───┐    ┌──────┐                      ┌───┐    ┌──────┐             ┌─┴─┐    ┌──────┐    ┌─┴─┐    ┌───┐    ┌──────┐             ┌──────┐    ┌─┴─┐    ┌──────┐    ┌─┴─┐    ┌─┴─┐    ┌──────┐    ┌─┴─┐    ┌─┴─┐    ┌──────┐    ┌─┴─┐  
|0⟩──┤ RotZ ├────┤ Phase ├──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤ H ├──────────────────────────────────┤ X ├────┤ RotZ ├────┤ X ├────┤ H ├──────────────────────────────────┤ H ├──────────────────────────────────┤ X ├────┤ RotZ ├────┤ X ├────┤ H ├──────────────────────────────────┤ Y ├────┤ RotX ├──────────────────────┤ X ├────┤ RotZ ├────┤ X ├────┤ Y ├────┤ RotX ├──────────────────────┤ Y ├────┤ RotX ├─────────────┤ X ├────┤ RotZ ├────┤ X ├────┤ Y ├────┤ RotX ├─────────────┤ RotZ ├────┤ X ├────┤ RotZ ├────┤ X ├────┤ X ├────┤ RotZ ├────┤ X ├────┤ X ├────┤ RotZ ├────┤ X ├──
     └──────┘    └───────┘                                                                                                                                  └───┘                                  └───┘    └──────┘    └───┘    └───┘                                  └───┘                                  └───┘    └──────┘    └───┘    └───┘                                  └───┘    └──────┘                      └───┘    └──────┘    └───┘    └───┘    └──────┘                      └───┘    └──────┘             └───┘    └──────┘    └───┘    └───┘    └──────┘             └──────┘    └───┘    └──────┘    └───┘    └───┘    └──────┘    └───┘    └───┘    └──────┘    └───┘  
```

We care about evolving the hamiltonian matrix for finding the minimum eigenvalue because it allows us to represent any hermitian matrix with unitary matrix operations
and use the property of eigenstates to find the smallest eigenvalue.

In essence, if we are searching the smallest eigenvalue $\lambda$ with eigenstate $|\psi\rangle$ of the matrix $\hat{H}$, then if we call $\hat{T}$ the evolution operator, for
the evolution circuit we would have
$$\hat{T}|\psi\rangle = e^{i\hat{H}t}|\psi\rangle = e^{i\lambda t}|\psi\rangle$$

In a quantum computer, we can't get the actual $|\psi\rangle$ state, but we can obtain the expectation value for this operation:
$$\langle\psi|\hat{T}|\psi\rangle = \langle\psi|e^{i\lambda t}|\psi\rangle = e^{i\lambda t}$$
which is just a complex number from which we can obtain the value $\lambda$ and use and optimizer to minimize this quantity.

All of this procedure is done inside the function `variational_quantum_eigensolver` so we just need to give the matrix we are trying to find the minimum eigenvalue of
with some other parameters and it will do the work for us.

For our case, we can set the search space of the eigenvalue to be $[-10, 10]$ by setting `range = 10`. For the variational state generator we can set either one as it is a simple
matrix, let's leave it as the default of `AQS::LINEAR`, for the tolerance let's leave it at 0 for maximum precision, and let's use 1000 iterations for the optimization algorithm.

```c++
    float range = 10.f;
    float tolerance = 0.f;
    uint32_t iterations = 1000;

    auto result = aqs::variational_quantum_eigensolver(hamiltonian_matrix, range, aqs::VQE::LINEAR, tolerance, iterations);

    auto ground_state = result.first;
```

After the algorithm runs, it will output a `std::pair` contain the minimum eigenvalue in the first result of the pair.
Thus the ground state energy for the $H_2$ is that minimum eigenvalue. You can expect the output to be around the value of $-1.125\text{ Ha}$
which is close to the expected $-1.1362 \text{ Ha}$ for the given arrangement of the molecule.

While running this algorithm in a classical computer may not be as efficient as other method, this goes to show the usefulness of quantum methods and algorithms that
can be useful to solve these kinds of problems, and having a library which gives you the tools to work with this kind of algorithms and develop their own algorithm goes the extra mile in
researching Quantum Computing.
