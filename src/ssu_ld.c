/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2023-2025, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <dlfcn.h>
#include <pthread.h>

/* Handle pointing to the approriate libssu implementing the functionality
 * Initialized on first use. */
static void *dl_handle = NULL;

static void ssu_load(const char *fncname,
                     void **dl_ptr) {
   char *error;

   if (dl_handle==NULL) {
       const char* lib_name = ssu_get_lib_name();
       const char* env_cpu_info = getenv("UNIFRAC_CPU_INFO");
       if ((env_cpu_info!=NULL) && (env_cpu_info[0]=='Y')) {
           printf("INFO (unifrac): Using shared library %s\n",lib_name);
       }
       dl_handle = dlopen(lib_name, RTLD_LAZY);
       if (!dl_handle) {
          fputs(dlerror(), stderr);
          exit(1);
       }
   }

   *dl_ptr = dlsym(dl_handle, fncname);
   if ((error = dlerror()) != NULL)  {
       fputs(error, stderr);
       exit(1);
   }
}

static pthread_mutex_t dl_mutex = PTHREAD_MUTEX_INITIALIZER;

static void cond_ssu_load(const char *fncname,
                     void **dl_ptr) {

   pthread_mutex_lock(&dl_mutex);
   if ((*dl_ptr)==NULL) ssu_load(fncname,dl_ptr);
   pthread_mutex_unlock(&dl_mutex);
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
/* This function is sometimes used when included, so prevent the compiler warning */
static bool ssu_load_check() {

   pthread_mutex_lock(&dl_mutex);
   if (dl_handle==NULL) {
       const char* lib_name = ssu_get_lib_name();
       dl_handle = dlopen(lib_name, RTLD_LAZY);
       if (!dl_handle) {
          /* Callers turn this into "GPU not detected", which conflates "no such
           * library" with "the library is there but would not load". Surface
           * dlerror() under either info flag so the two are distinguishable;
           * without it, diagnosing a runner that stops seeing its GPU takes
           * several CI runs of correlation.
           */
          const char* dl_msg = dlerror();
          const char* env_cpu_info = getenv("UNIFRAC_CPU_INFO");
          const char* env_gpu_info = getenv("UNIFRAC_GPU_INFO");
          if (((env_cpu_info!=NULL) && (env_cpu_info[0]=='Y')) ||
              ((env_gpu_info!=NULL) && (env_gpu_info[0]=='Y'))) {
              printf("INFO (unifrac): Could not load shared library %s: %s\n",
                     lib_name, (dl_msg!=NULL) ? dl_msg : "unknown error");
          }
          pthread_mutex_unlock(&dl_mutex);
	  return false;
       }
       // only print out if the library exists
       const char* env_cpu_info = getenv("UNIFRAC_CPU_INFO");
       if ((env_cpu_info!=NULL) && (env_cpu_info[0]=='Y')) {
           printf("INFO (unifrac): Using shared library %s\n",lib_name);
       }
   }
   pthread_mutex_unlock(&dl_mutex);
   return true;
}
#pragma GCC diagnostic push


