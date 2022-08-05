Usage
======

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

### Construction
In this library, `QState` represents the quantum state of a qubit expressed with $|0\rangle$ and $|1\rangle$ in the computational basis
It stores the normalized components of the basis vectors.

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
of taking advantage of the entangling multiple qubits to take advantage.

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
affects the state of the qubits, and then meaasure them.

## Quantum Algorithms
For commonly used quantum algorithms, the library provides functions that create these algorithms in circuits.

All implemented algorithms are contained in the file [`quantum_algo.h`](include/quantum.h).

### Quantum Fourier Transform


## Special Gates

To facilitate creation of complex gates like Multiple Control Phase Gate, the Adjoint of a Gate,
or adding multiple gates at the same time, the library facilitates multiple functions that create
those circuits.

All implemented special gates are contained in the file [`quantum_gates.h`](include/quantum.h)

### Group_Gate

### Control_Group_Gate

### NControlGate

### Rewire_Gate

### Adjoint_Gate

## Visuals

For displaying measurements or circuits, the library provides many functions to print the information to the standard outstream (console/terminal).
This functions are located in [`quantum_visuals.h`](include/quantum_visuals.h).

### 