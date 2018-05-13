#include <cstdio>
#include <cmath>
#include <omp.h>
#include "particlematrix.h"
#include "common.h"


double i_start, i_total = 0;
double c_start, c_total = 0;
double m_start, m_total = 0;
double s_start, s_total = 0;

ParticleMatrix::ParticleMatrix(int n) {
    nof_particles = n;
    nof_slices = (n > 10) ? n/10 : 1;

    size = sqrt(0.0005 * n);

    particle_matrix = new particle_vector_t[nof_slices];
    particles = new particle_t[nof_particles];

    init_particles(nof_particles, particles);
}

void ParticleMatrix::perform_steps(int n, bool perform_save) {

#pragma omp parallel num_threads(4)
    for (int step = 0; step < n; step ++) {

        #pragma omp master
        {
            i_start = omp_get_wtime();

            for (int i = 0; i < nof_slices; i++) {
                particle_matrix[i].clear();
            }

            for (int i = 0; i < nof_particles; i++) {
                particle_matrix[(int) (particles[i].x * floor(nof_slices / size))].emplace_back(&particles[i]);
            }

            i_total += omp_get_wtime() - i_start;

            c_start = omp_get_wtime();
        }
        #pragma omp barrier

        #pragma omp for schedule(guided)
        for (int slice = 0; slice < nof_slices; slice++) {
            int start = ((slice > 0) ? (slice - 1) : slice);
            int end = ((slice < nof_slices - 1) ? (slice + 2) : (slice + 1));

            for (particle_t * curr_particle: particle_matrix[slice]) {
                (*curr_particle).ax = (*curr_particle).ay = 0;

                for (int j = start; j < end; j++) {
                    for (particle_t * neigh_particle: particle_matrix[j]) {
                        apply_force((*curr_particle), (*neigh_particle));
                    }
                }
            }
        }

#pragma omp master
        {
            c_total += omp_get_wtime() - c_start;

            m_start = omp_get_wtime();
        }

#pragma omp for
        for (int i = 0; i < nof_particles; i++) {
            move(particles[i]);
        }


#pragma omp master
        {

            m_total += omp_get_wtime() - m_start;

            if(perform_save) {
                s_start = omp_get_wtime();
                save(particles);

                s_total += omp_get_wtime() - s_start;
            }

        }

    }
}
void ParticleMatrix::print(){
    printf("i time: %f, c time: %f, m time: %f, save_time: %f\n", i_total, c_total, m_total, s_total);

}
