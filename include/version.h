/*******************************************************
 * Copyright (c) 2022, ArrayFire
 * All rights reserved.
 *
 * This file is distributed under 3-clause BSD license.
 * The complete license agreement can be obtained at:
 * http://arrayfire.com/licenses/BSD-3-Clause
 ********************************************************/

/**
 * @file version.h
 * @brief Defines the current version of the AQS library code
 */

#pragma once

/**
 * @def AQS_VERSION
 * @showinitializer
 */
#define AQS_VERSION "1.0.0"

/**
 * @def AQS_VERSION_MAJOR
 * @showinitializer
 */
#define AQS_VERSION_MAJOR 1

/**
 * @def AQS_VERSION_MINOR
 * @showinitializer
 */
#define AQS_VERSION_MINOR 0

/**
 * @def AQS_VERSION_PATCH
 * @showinitializer
 */
#define AQS_VERSION_PATCH 0

namespace aqs
{
    /**
     * @brief Returns a string of the AQS version in major, minor, patch order separated by dots
     * 
     */
    inline const char* get_version()
    {
        return AQS_VERSION;
    }
}