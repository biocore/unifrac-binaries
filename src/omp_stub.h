/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2026, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#ifndef UNIFRAC_OMP_STUB_H
#define UNIFRAC_OMP_STUB_H

#if defined(_OPENMP)
#include <omp.h>
#else

#ifdef __cplusplus
extern "C" {
#endif

static inline int omp_get_max_threads(void) { return 1; }
static inline int omp_get_thread_num(void)  { return 0; }
static inline int omp_get_num_devices(void) { return 0; }

#ifdef __cplusplus
}
#endif

#endif /* _OPENMP */
#endif /* UNIFRAC_OMP_STUB_H */
