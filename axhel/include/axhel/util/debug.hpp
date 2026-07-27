// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "axhel/util/defines.hpp"

#ifdef AXHEL_DEBUG

#include <iostream>

#define AXHEL_LOG(message)            \
    do {                                \
        std::cerr                       \
            << "[AXHEL] "               \
            << message                  \
            << '\n';                    \
    } while (false)

#else

#define AXHEL_LOG(message) \
    do {                     \
    } while (false)

#endif