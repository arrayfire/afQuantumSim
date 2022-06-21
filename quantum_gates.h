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
    QCircuit Control_GroupPhase(uint32_t qubits, uint32_t control_qubit, uint32_t target_qubit_begin, float angle);
    QCircuit NControl_Gate(uint32_t qubits, uint32_t control_qubit_begin, uint32_t control_qubit_count, uint32_t target_qubit_begin, const QCircuit& gate);
}