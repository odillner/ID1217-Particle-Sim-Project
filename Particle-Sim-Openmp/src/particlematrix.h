#ifndef PARTICLEMATRIX_H
#define PARTICLEMATRIX_H

#include "common.h"
#include <vector>


class ParticleMatrix {
public:
    ParticleMatrix (int);               /* constructor, accepts number of particles and simulation size */
    void perform_steps(int, bool);      /* performs n steps of simulation */

    void print();
private:
    typedef std::vector <particle_t *>  particle_vector_t;

    particle_t * particles;             /* main particle array */
    particle_vector_t * particle_matrix;
    int nof_slices;                     /* number of slices */
    int nof_particles;                  /* number of particles */
    double size;

};

#endif //PARTICLEMATRIX_H

