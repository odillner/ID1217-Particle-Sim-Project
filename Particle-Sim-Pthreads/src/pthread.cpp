#include <cmath>
#include <unistd.h>
#include "pthread.h"

particle_t *particles;               /* main particle array */
particle_vector_t *particle_matrix;  /* particle data structure */

int nof_slices;                      /* number of slices */
int nof_particles;                   /* number of particles */
double size;

pthread_t threads[NOF_THREADS];

pthread_mutex_t barrier_lock;
pthread_cond_t barrier_cond;

int num_arrived = 0;

bool p_save = false;

void Barrier() {
    pthread_mutex_lock(&barrier_lock);
    num_arrived++;
    if (num_arrived == NOF_THREADS) {
        num_arrived = 0;
        pthread_cond_broadcast(&barrier_cond);
    } else {
        pthread_cond_wait(&barrier_cond, &barrier_lock);
    }
    pthread_mutex_unlock(&barrier_lock);
}

void *Worker(void *arg) {
    int pid = (int) arg;

    int slice_length = (nof_slices + NOF_THREADS - 1) / NOF_THREADS;
    int start = pid * slice_length;
    int end = (pid == NOF_THREADS - 1) ? (nof_slices) : start + slice_length;

    for (int step = 0; step < NSTEPS; step++) {

        if (pid == 0) {

            for (int i = 0; i < nof_slices; i++) {
                particle_matrix[i].clear();
            }

            for (int i = 0; i < nof_particles; i++) {
                particle_matrix[(int) (particles[i].x * floor(nof_slices / size))].emplace_back(&particles[i]);
            }

        }

        Barrier();


        for (int i = start; i < end; i++) {
            int comp_start = ((i > 0) ? (i - 1) : i);
            int comp_end = ((i < nof_slices - 1) ? (i + 2) : (i + 1));

            for (particle_t *curr_particle : particle_matrix[i]) {
                (*curr_particle).ax = (*curr_particle).ay = 0;

                for (int j = comp_start; j < comp_end; j++) {
                    for (particle_t *neigh_particle: particle_matrix[j]) {
                        apply_force((*curr_particle), (*neigh_particle));
                    }
                }
            }
        }

        for (int i = start * 10; i < end * 10; i++) {
            move(particles[i]);
        }

        Barrier();

        if (pid == 0) {
            if (p_save) {
                save(particles);
            }
        }

        Barrier();
    }
}

void perform_steps(int n, bool perform_save) {
    nof_particles = n;
    nof_slices = (n < 10) ? 1 : n / 10;
    size = sqrt(0.0005 * n);

    p_save = perform_save;

    particle_matrix = new particle_vector_t[nof_slices];
    particles = new particle_t[nof_particles];

    init_particles(nof_particles, particles);

    pthread_cond_init(&barrier_cond, NULL);
    pthread_mutex_init(&barrier_lock, NULL);

    for (int i = 0; i < NOF_THREADS; i++) {
        pthread_create(&threads[i], NULL, Worker, (void *) i);
    }
    for (int i = 0; i < NOF_THREADS; i++) {
        pthread_join(threads[i], NULL);
    }
}