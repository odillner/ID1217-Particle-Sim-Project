#ifndef PARTICLE_SIM_PTHREADS_PTHREAD_H
#define PARTICLE_SIM_PTHREADS_PTHREAD_H

#include <vector>
#include "common.h"

#define NOF_THREADS 4

typedef std::vector <particle_t *>  particle_vector_t;

void perform_steps(int, bool);      /* performs n steps of simulation */

#endif


