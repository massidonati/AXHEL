// Copyright (C) 2026 University of Pisa - Dept. of Information Engineering
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "axhel/util/defines.hpp"

#define AXHEL_PRAGMA_IMPL(x) _Pragma(#x)
#define AXHEL_PRAGMA(x) AXHEL_PRAGMA_IMPL(x)

#if defined(AXHEL_USE_GNU)

#define AXHEL_UNROLL(N) AXHEL_PRAGMA(GCC unroll N)

#elif defined(AXHEL_USE_CLANG)

#define AXHEL_UNROLL(N) AXHEL_PRAGMA(clang loop unroll_count(N))

#else

#define AXHEL_UNROLL(N)

#endif