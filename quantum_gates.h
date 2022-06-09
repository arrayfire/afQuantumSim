/*******************************************************
 * Copyright (c) 2022, ArrayFire
 * All rights reserved.
 *
 * This file is distributed under 3-clause BSD license.
 * The complete license agreement can be obtained at:
 * http://arrayfire.com/licenses/BSD-3-Clause
 ********************************************************/

#pragma once

#include "quantum.h"

namespace aqs
{
    QCircuit Control_GroupPhase(int qubits, int control_qubit, int target_qubit_begin, float angle);
    QCircuit NControl_Gate(int qubits, int control_qubit_begin, int control_qubit_count, int target_qubit_begin, const QCircuit& gate);
}