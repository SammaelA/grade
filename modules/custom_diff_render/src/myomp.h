#pragma once

#ifdef USE_OMP
#include <omp.h>
#else
static int omp_get_thread_num() {return 1;}
#endif