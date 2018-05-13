#include <mpi.h>
#include <cmath>
#include <vector>
#include "common.h"



//
//  benchmarking program
//
void master(int nof_proc, int nof_slices, int nof_particles, int *partition_sizes, int *actual_partition_sizes, FILE *fsave) {
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous(6, MPI_DOUBLE, &PARTICLE);
    MPI_Type_commit(&PARTICLE);

    double size = get_size();
    particle_t *particles = new particle_t[nof_particles];
    particle_vector_t *particle_matrix = new particle_vector_t[nof_slices];
    int *slice_sizes = new int[nof_slices];
    int temp;

    init_particles(nof_particles, particles);

    for (int step = 0; step < NSTEPS; step++) {

        for (int i = 0; i < nof_slices; i++) {
            particle_matrix[i].clear();
            slice_sizes[i] = 0;
        }
        for (int i = 0; i < nof_particles; i++) {
            temp = (int) (particles[i].x * floor(nof_slices / size));
            particle_matrix[temp].emplace_back(&particles[i]);
            slice_sizes[temp]++;
        }
        temp = 0;

        for (int i = 0; i < nof_proc; i++) {

            for (int j = 0; j < partition_sizes[i]; j++) {

                int curr_slice = temp + j;
                MPI_Send(&slice_sizes[curr_slice], 1, MPI_INT, i + 1, 0, MPI_COMM_WORLD);
                for (particle_t *curr_particle : particle_matrix[curr_slice]) {
                    MPI_Send(curr_particle, 1, PARTICLE, i + 1, 0, MPI_COMM_WORLD);
                }
            }

            temp += partition_sizes[i] - 2;

            if (i != 0 && i != nof_proc-1) {
                temp--;
            }

        }

        if( fsave && (step%SAVEFREQ) == 0 )
            save( fsave, nof_particles, particles );

        temp = 0;
        for (int i = 0; i < nof_proc; i++) {
            for (int k = 0; k < actual_partition_sizes[i]; k++) {
                MPI_Recv(&particles[temp], 1, PARTICLE, i + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                temp++;
            }
        }

    }

}

void worker(int rank, int nof_proc, int partition){
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );

    std::vector<particle_t> local;
    particle_vector_t * local_matrix = new particle_vector_t[partition];
    int * slice_sizes = new int[partition];
    particle_t curr_particle;
    int start_slice = (rank == 1) ? 0 : 1;
    int end_slice = (rank == nof_proc-1) ? partition : partition-1;
    int temp = 0;


    for (int step = 0; step < NSTEPS; step++) {
        local.clear();

        for (int i = 0; i < partition; i++) {
            MPI_Recv(&slice_sizes[i], 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int j = 0; j < slice_sizes[i]; j++) {
                MPI_Recv(&curr_particle, 1, PARTICLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                local.emplace_back(curr_particle);
            }
        }

        temp = 0;

        for (int i = 0; i < partition; i++) {
            local_matrix[i].clear();
            for (int j = 0; j < slice_sizes[i]; j++) {
                local_matrix[i].emplace_back(&local[j + temp]);
            }
            temp += slice_sizes[i];
        }

        for (int slice = start_slice; slice < end_slice; slice++) {
            for (particle_t *curr_particle: local_matrix[slice]) {
                (*curr_particle).ax = (*curr_particle).ay = 0;

                for (int j = 0; j < partition; j++) {
                    for (particle_t *neigh_particle: local_matrix[j]) {
                        apply_force((*curr_particle), (*neigh_particle));
                    }
                }
            }
        }


        temp = 0;
        for (int slice = start_slice; slice < end_slice; slice++) {
            temp += slice_sizes[slice];
        }

        int start_index = start_slice * slice_sizes[0];
        int end_index = temp + start_index;


        for (int i = start_index; i < end_index; i++) {
            move(local[i]);
        }


        for (int i = start_slice; i < end_slice; i++) {
            for (int j = 0; j < slice_sizes[i]; j++) {
                MPI_Send(local_matrix[i][j], 1, PARTICLE, 0, 0, MPI_COMM_WORLD);
            }
        }

    }

}

int main( int argc, char **argv )
{
    //
    //  process command line parameters
    //
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );
    char * savename = read_string( argc, argv, "-o", NULL );

    //
    //  set up MPI
    //
    int n_proc, rank;

    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );


    //
    //  allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;

    //
    //  set up the data partitioning across processors
    //

    int nof_slices = n / 10;

    n_proc--;

    int particle_per_proc = (nof_slices + n_proc - 1) / n_proc;
    int *partition_offsets = (int*) malloc( (n_proc+1) * sizeof(int) );
    for( int i = 0; i < n_proc+1; i++ )
        partition_offsets[i] = min( i * particle_per_proc, nof_slices );

    int *partition_sizes = (int*) malloc( n_proc * sizeof(int) );
    int *actual_partition_sizes = (int*) malloc( n_proc * sizeof(int) );
    for( int i = 0; i < n_proc; i++ ) {
        actual_partition_sizes[i] = partition_offsets[i + 1] - partition_offsets[i];
        partition_sizes[i] = actual_partition_sizes[i];
        if (n_proc > 1) {
            partition_sizes[i]++;
            if (i != 0 && i != n_proc - 1) {
                partition_sizes[i]++;
            }
        }
    }

    //
    //  allocate storage for local partition
    //
    int nlocal = partition_sizes[rank-1];

    set_size(n);

    double simulation_time = read_timer( );

    if (rank == 0) {
        producer(n_proc, nof_slices, n, partition_sizes, actual_partition_sizes, fsave);
    } else {
        consumer(rank, n_proc, nlocal);
    }


    if( rank == 0 ) {
        simulation_time = read_timer() - simulation_time;
        printf("n = %d, n_procs = %d, simulation time = %g s\n", n, n_proc+1, simulation_time);
    }

    if( fsave )
        fclose( fsave );

    MPI_Finalize( );

    return 0;
}


