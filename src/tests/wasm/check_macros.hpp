/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2026, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 * Minimal test helpers for the WASM correctness suite. Mirrors the bare-
 * main+CHECK style used in scikit-bio-binaries' WASM tests rather than
 * pulling in a framework. Intentionally tiny: the WASM build deliberately
 * has no test framework dependency.
 */

#ifndef UNIFRAC_WASM_CHECK_MACROS_HPP
#define UNIFRAC_WASM_CHECK_MACROS_HPP

#include <cmath>
#include <cstdio>
#include <cstdlib>

#define CHECK(cond)                                                     \
    do {                                                                \
        if (!(cond)) {                                                  \
            std::fprintf(stderr, "FAIL %s:%d: %s\n",                    \
                         __FILE__, __LINE__, #cond);                    \
            std::exit(1);                                               \
        }                                                               \
    } while (0)

#define CHECK_EQ(a, b)                                                  \
    do {                                                                \
        auto _a = (a);                                                  \
        auto _b = (b);                                                  \
        if (!(_a == _b)) {                                              \
            std::fprintf(stderr, "FAIL %s:%d: " #a " == " #b            \
                         " (got %ld vs %ld)\n",                         \
                         __FILE__, __LINE__,                            \
                         (long)_a, (long)_b);                           \
            std::exit(1);                                               \
        }                                                               \
    } while (0)

template <class T>
static inline bool almost_equal(T a, T b, T abs_tol) {
    return std::fabs(a - b) <= abs_tol;
}

#endif /* UNIFRAC_WASM_CHECK_MACROS_HPP */
