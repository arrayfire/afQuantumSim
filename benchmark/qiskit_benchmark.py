#######################################
# Code for the algorithms taken from the Qiskit Textbook
# See license in docs/licenses/qiskit
#######################################

import time
import numpy as np
from numpy.random import randint
from fractions import Fraction
from math import gcd

from qiskit import Aer, assemble, transpile
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.circuit.controlledgate import ControlledGate
from qiskit.visualization import plot_histogram

def grover(qubits : int):
    solution_count = 1
    search_space = 1 << qubits
    iterations = int(np.pi * np.sqrt(search_space / solution_count) / 4)

    marked_state = 5
    grover_circuit = QuantumCircuit(qubits)

    qubit_list = []
    for i in range(qubits - 1):
        qubit_list.append(i)

    for i in range(qubits):
        grover_circuit.h(i)

    for j in range(iterations):
        for i in range(qubits):
            if marked_state & (1 << i):
                grover_circuit.x(i)
        
        grover_circuit.mcp(np.pi, qubit_list , qubits - 1)

        for i in range(qubits):
            if marked_state & (1 << i):
                grover_circuit.x(i)

        for i in range(qubits):
            grover_circuit.h(i)
            grover_circuit.x(i)

        grover_circuit.mcp(np.pi, qubit_list, qubits - 1)

        for i in range(qubits):
            grover_circuit.x(i)
            grover_circuit.h(i)

    sim = Aer.get_backend('aer_simulator')

    grover_circuit.save_statevector()
    qobj = assemble(grover_circuit)
    result = sim.run(qobj).result()
    #print(result.get_counts())


def qft(qubits : int):
    qft_circuit = QuantumCircuit(qubits)

    for i in range(qubits):
        qft_circuit.h(qubits - 1 - i)
        for j in range(qubits - 1 - i):
            qft_circuit.cp(np.pi / (1 << (qubits - 1 - i - j)), j, qubits - 1 - i)

    for i in range(qubits // 2):
        qft_circuit.swap(i, qubits - 1 - i)

    sim = Aer.get_backend('aer_simulator')

    qft_circuit.save_statevector()
    qobj = assemble(qft_circuit)
    result = sim.run(qobj).result()


def shor(qubits : int):
    N = 15
    np.random.seed(1) # This is to make sure we get reproduceable results

    def c_amod15(a, power):
        """Controlled multiplication by a mod 15"""
        if a not in [2,4,7,8,11,13]:
            raise ValueError("'a' must be 2,4,7,8,11 or 13")
        U = QuantumCircuit(4)        
        for iteration in range(power):
            if a in [2,13]:
                U.swap(0,1)
                U.swap(1,2)
                U.swap(2,3)
            if a in [7,8]:
                U.swap(2,3)
                U.swap(1,2)
                U.swap(0,1)
            if a in [4, 11]:
                U.swap(1,3)
                U.swap(0,2)
            if a in [7,11,13]:
                for q in range(4):
                    U.x(q)
        U = U.to_gate()
        U.name = "%i^%i mod 15" % (a, power)
        c_U = U.control()
        return c_U

    def qft_dagger(n):
        """n-qubit QFTdagger the first n qubits in circ"""
        qc = QuantumCircuit(n)
        # Don't forget the Swaps!
        for qubit in range(n//2):
            qc.swap(qubit, n-qubit-1)
        for j in range(n):
            for m in range(j):
                qc.cp(-np.pi/float(2**(j-m)), m, j)
            qc.h(j)
        qc.name = "QFT†"
        return qc

    def qpe_amod15(a):
        n_count = 8
        qc = QuantumCircuit(4+n_count, n_count)
        for q in range(n_count):
            qc.h(q)     # Initialize counting qubits in state |+>
        qc.x(3+n_count) # And auxiliary register in state |1>
        for q in range(n_count): # Do controlled-U operations
            qc.append(c_amod15(a, 2**q), 
                    [q] + [i+n_count for i in range(4)])
        qc.append(qft_dagger(n_count), range(n_count)) # Do inverse-QFT
        qc.measure(range(n_count), range(n_count))
        # Simulate Results
        aer_sim = Aer.get_backend('aer_simulator')
        # Setting memory=True below allows us to see a list of each sequential reading
        t_qc = transpile(qc, aer_sim)
        qobj = assemble(t_qc, shots=1)
        result = aer_sim.run(qobj, memory=True).result()
        readings = result.get_memory()
        # print("Register Reading: " + readings[0])
        phase = int(readings[0],2)/(2**n_count)
        # print("Corresponding Phase: %f" % phase)
        return phase

    a = 7
    factor_found = False
    attempt = 0
    while not factor_found:
        attempt += 1
        # print("\nAttempt %i:" % attempt)
        phase = qpe_amod15(a) # Phase = s/r
        frac = Fraction(phase).limit_denominator(N) # Denominator should (hopefully!) tell us r
        r = frac.denominator
        # print("Result: r = %i" % r)
        if phase != 0:
            # Guesses for factors are gcd(x^{r/2} ±1 , 15)
            guesses = [gcd(a**(r//2)-1, N), gcd(a**(r//2)+1, N)]
            # print("Guessed Factors: %i and %i" % (guesses[0], guesses[1]))
            for guess in guesses:
                if guess not in [1,N] and (N % guess) == 0: # Check to see if guess is a factor
                    # print("*** Non-trivial factor found: %i ***" % guess)
                    factor_found = True

def benchmark(qubits : int, count : int, func):
    total = 0
    total_squared = 0
    for _ in range(count):
        start_time = time.perf_counter()
        func(qubits)
        end_time = time.perf_counter()

        duration = end_time - start_time
        total += duration
        total_squared += duration * duration

    average_time = total / count
    average_squared_time = total_squared / count
    sstddev = (average_squared_time - average_time * average_time) * count / (count - 1)

    print(func.__name__)
    print(f"Total Time: {total} s; Avg. Time: {average_time * 1000} ms, S StdDev: {sstddev * 1000} ms; Test Count: {count}")

def main():
    qubits = 10
    count = 100

    benchmark(qubits, count, qft)
    benchmark(qubits, count, grover)
    benchmark(qubits, count, shor)

if __name__ == '__main__':
    main()