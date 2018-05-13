#include <cstdio>
#include <cmath>
#include "particlematrix.h"

double i_start, i_total = 0;
double c_start, c_total = 0;
double m_start, m_total = 0;
double s_start, s_total = 0;

ParticleMatrix::ParticleMatrix(int n) {
    nof_particles = n;
    nof_slices = n/10;

    size = sqrt(0.0005 * n);

    particle_matrix = new particle_vector_t[nof_slices];
    particles = new particle_t[nof_particles];
    init_particles(nof_particles, particles);
}

void ParticleMatrix::index_particles() {

    for (int i = 0; i < nof_slices; i++) {
        particle_matrix[i].clear();
    }

    for (int i = 0; i < nof_particles; i++) {
        int curr_particle_slice = (particles[i].x*floor(nof_slices/size));

        particle_matrix[curr_particle_slice].emplace_back(&particles[i]);
    }

}
void ParticleMatrix::collision_check() {
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


}

void ParticleMatrix::perform_steps(int n, bool perform_save) {

    for (int steps = 0; steps < n; steps++) {
        i_start = read_timer( );
        index_particles();
        i_total += read_timer( ) - i_start;

        c_start = read_timer( );
        collision_check();
        c_total += read_timer( ) - c_start;

        m_start = read_timer( );
        for (int i = 0; i < nof_particles; i++) {
            move(particles[i]);
        }
        m_total += read_timer( ) - m_start;

        if(perform_save) {
            s_start = read_timer( );
            save(particles);

            s_total += read_timer( ) - s_start;
        }
    }

}
void ParticleMatrix::print(){
    printf("i time: %f, c time: %f, m time: %f, save_time: %f\n", i_total, c_total, m_total, s_total);
}