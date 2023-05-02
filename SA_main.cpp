#include <iostream>
#include <mpi.h>
#include "annealer.h"

int main(int argc, char** argv) {

    // create points representing the unit square
    // double x[4] = {1., 1., 0., 0.};
    // double y[4] = {1., 0., 1., 0.};

    // double* x11 = (double*) malloc(4 * sizeof(double));
    // double* x12 = (double*) malloc(4 * sizeof(double));
    // for (int i = 0; i < 4; i++) {
    //     x11[i] = x[i];
    //     x12[i] = x[i];
    // }

    // double* x21 = (double*) malloc(4 * sizeof(double));
    // double* x22 = (double*) malloc(4 * sizeof(double));
    // for (int i = 0; i < 4; i++) {
    //     x21[i] = y[i];
    //     x22[i] = y[i];
    // }

    // long N = 4;

    // // create states and a transition
    // TSP2DState state1 = TSP2DState(N, x11, x21);
    // TSP2DState state2 = TSP2DState(N, x12, x22);
    // TSP2DTransition transition = TSP2DTransition(2, 3);

    // printf("TRANSITION:\n");

    // // print the current path and the path length (the objective value)
    // state1.display_state();
    // state1.display_coords();
    // printf("obj = %f\n", state1.objective());

    // // use the transition to evolve the state
    // state1.step(&transition);

    // state1.display_state();
    // printf("obj = %f\n", state1.objective());

    // printf("ENERGY:\n");

    // printf("energy (1->2) = %f\n", state1.energy(&state2));
    // printf("energy (2->1) = %f\n", state2.energy(&state1));

    // printf("energy local test (1->2) = %f\n", state2.energy_local(transition));

    // printf("FROM FILE:\n");

    // TSP2DState from_file = TSP2DState::from_text_file("example_tsp_in.txt");
    // TSP2DTransition transition1 = TSP2DTransition(0, 1);
    // TSP2DTransition transition2 = TSP2DTransition(71, 49);

    // printf("obj = %f\n", from_file.objective());
    // from_file.step(&transition1);
    // from_file.step(&transition2);
    // printf("obj = %f\n", from_file.objective());

    // Annealer TSP = Annealer(from_file.stops(), LOG);
    // TSP.generate_t_matrix(&from_file);
    // TSP.display_params();
    // TSP.anneal(&from_file);

    // from_file.write_txt("example_tsp_out.txt");

    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    int mpirank, mpisize;
    MPI_Comm_rank(comm, &mpirank);
    MPI_Comm_size(comm, &mpisize);

    // srand(mpirank);
    srand48(mpirank);

    MPI_Status status;

    TSP2DState parallel_state = TSP2DState::from_text_file("tsp_examples/spread=1.0/100.txt");
    Annealer parallel_annealer = Annealer(parallel_state.stops(), LOG);

    double min_objective = parallel_annealer.anneal(&parallel_state, 1000);
    std::vector<long> min_state = parallel_annealer.get_min_state();
    long size = parallel_state.stops();

    MPI_Barrier(comm);

    if (mpirank != 0) {
        MPI_Send(&min_objective, 1, MPI_DOUBLE, 0, 999, comm);
        MPI_Send(min_state.data(), size, MPI_LONG, 0, 999, comm);
    }

    if (mpirank == 0) {
        double recv_min_objective;
        long* recv_min_state = (long*) malloc(size * sizeof(long));
        int min_idx;
        for (int i = 1; i < mpisize; ++i) {
            MPI_Recv(&recv_min_objective, 1, MPI_DOUBLE, i, 999, comm, &status);
            MPI_Recv(recv_min_state, size, MPI_LONG, i, 999, comm, &status);
            // Rank 0 min_objective already computed
            if (recv_min_objective < min_objective) {
                min_objective = recv_min_objective;
                for (long j = 0; j < size; j++)
                    min_state[j] = recv_min_state[j];
                min_idx = i;
            }
        }
        // Rank 0 broadcasts the min state
        MPI_Bcast(min_state.data(), size, MPI_LONG, 0, comm);
        free(recv_min_state);
    }

    MPI_Barrier(comm);

    //Everyone's min_state should be the network minimum;
    parallel_state.set_idxs(min_state);

    std::cout<< "Hello from rank " << mpirank << ", we got annealed." << std::endl;

    MPI_Finalize();
}
