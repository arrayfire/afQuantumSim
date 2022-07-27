/*******************************************************
 * Copyright (c) 2022, ArrayFire
 * All rights reserved.
 *
 * This file is distributed under 3-clause BSD license.
 * The complete license agreement can be obtained at:
 * http://arrayfire.com/licenses/BSD-3-Clause
 ********************************************************/

#include "quantum.h"
#include "quantum_visuals.h"

#include "utils.h"

#include <iostream>

int main(int argc, char** argv)
{
    aqs::initialize(argc, argv);
    std::cout << af::infoString() << std::endl;

    aqs::QCircuit qc{2};
    qc << aqs::X{0} << aqs::H{0} << aqs::H{1};
    qc.compile();

    aqs::QSimulator qs(2);
    qs.simulate(qc);

    std::cout << "Statevector in the z-basis (|0> and |1>): \n";
    aqs::print_statevector(qs);

    qs.set_basis(aqs::QSimulator::Basis::X);

    std::cout << "Statevector in the x-basis (|+> and |->x): \n";
    aqs::print_statevector(qs);

    qs.set_basis(aqs::QSimulator::Basis::Y);

    std::cout << "Statevector in the y-basis (|+>_y and |->_y): \n";
    aqs::print_statevector(qs);

    return 0;
}