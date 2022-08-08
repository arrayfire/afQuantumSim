Usage
======

# Table of Contents
1. [Initialization](#initialization)
2. [Quantum States](#quantum-states)
    1. [Construction](#construction)
    2. [Operations](#operations)
3. [Quantum Circuit](#quantum-circuit)
    1. [Construction](#construction-1)
    2. [Gates](#gates)
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

In this library, `QState` represents the quantum state of a qubit expressed with $|0\rangle$ and $|1\rangle$ in the computational basis
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

So for simulating a quantum circuit, the library provides `QCircuit` class which manages the creation of a circuit with
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

### Results
After the simulation is done, we would like to obtain results in the form of measurements, probabilities, or even the actual resulting statevector.
For that, the `QSimulator` class contains useful methods for accomplishing these tasks:

- `statevector`: This function returns the internal statevector resulting from the circuit simulation. Useful when it is needed to see the actual states of the registers.
- `probabilities`: Returns the probabilities of all the possible measurements that can be done in the computational basis from the statevector.
- `qubit_probability_true`/`qubit_probability_false`: Returns the probability of measuring a single qubit `true` ($|1\rangle$) or `false` ($|0\rangle$)
in the compuatational basis
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

    // Compile the circuit
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


### Hamiltonian Decomposition


### Hamiltonian Evolution
Another famous use for a quantum computer is its ability to simulate a quantum system in polynomial time by 

### Variational Quantum Eigensolver



## Special Gates

To facilitate creation of complex gates like Multiple Control Phase Gate, the Adjoint of a Gate,
or adding multiple gates at the same time, the library facilitates multiple functions that create
those circuits.

All implemented special gates are contained in the file [`quantum_gates.h`](../include/quantum.h)

- `Group_Gate`:
- `Control_Group_Gate`:
- `NControlGate`:
- `Rewire_Gate`:
- `Adjoint_Gate`:

## Visuals

For displaying measurements or circuits, the library provides many functions to print the information to the standard outstream (console/terminal).
This functions are located in [`quantum_visuals.h`](../include/quantum_visuals.h).

### State


### Profile Results


### Circuit Matrix


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