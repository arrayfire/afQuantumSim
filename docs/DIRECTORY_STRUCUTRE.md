Directory Structure
===================

```
.
|-- .gitignore
|-- LICENSE.txt
|-- README.md
|
|-- CMakeLists.txt
|
|-- benchmark/
|   |-- CMakeLists.txt
|   |-- benchmark.cpp
|   |-- qpp_benchmark.cpp
|   |-- qiskit_benchmark.py
|   |-- results.md
|
|-- docs/
|   |-- DIRECTORY_STRUCTURE.md
|   |-- REQUIREMENTS.md
|   |-- TODO.md
|   |-- USAGE.md
|   |-- licenses/
|       |-- nlopt/
|           |-- <licenses and author files>
|       |-- qiskit/
|           |-- <licenses and author files>
|       |-- Quantum++/
|           |-- <licenses and author files>
|
|-- examples/
|   |-- CMakeLists.txt
|   |-- <C++ source files of examples>
|
|-- include/
|   |-- quantum_algo.h
|   |-- quantum_gates.h
|   |-- quantum_visuals.h
|   |-- quantum.h
|   |-- utils.h
|
|-- src/
|   |-- quantum_algo.cpp
|   |-- quantum_gates.cpp
|   |-- quantum_visuals.cpp
|   |-- quantum.cpp
|   |-- utils.cpp
|
|-- test/
|   |-- CMakeLists.txt
|   |-- stress_test.cpp
|   |-- stress_test.h
|   |-- tests.cpp
|
```