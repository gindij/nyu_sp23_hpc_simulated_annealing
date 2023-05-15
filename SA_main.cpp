#include <iostream>
#include <mpi.h>
#include <unistd.h>
#include <string>
#include "annealer.h"


// command line args:
//    -x : choose file from tsp_examples/spread=1.0
//    -n : set maximum number of iterations in while loop
//    -i : set maximum number of iterations within annealer
//    -j : set number of steps per iterations in annealer
//    -t : set tolerance 
int main(int argc, char** argv) {

    int c;
    std::string filepath = "tsp_examples/spread=1.0/";
    long MAX_ITERATIONS = 100;
    long MAX_ANNEALER_ITERATIONS = 10000;
    long ANNEALING_STEPS_PER_ITERATION = 100;
    double TOLERANCE = 1e-6;

    while((c = getopt(argc, argv, "x:n:i:j:t:")) != -1) {
        switch(c) {
            case 'x':
                filepath = filepath + optarg;
                break;
            case 'n':
                MAX_ITERATIONS = atol(optarg);
                break;
            case 'i':
                MAX_ANNEALER_ITERATIONS = atol(optarg);
                break;
            case 'j':
                ANNEALING_STEPS_PER_ITERATION = atol(optarg);
                break;
            case 't':
                TOLERANCE = std::stod(optarg);
                break;
            case '?':
                if (optopt == 'x' || optopt == 'n' || optopt == 'i' || optopt == 'j' || optopt == 't') {
                    std::cout << "Option " << char(optopt) << " requires an argument" << std::endl;
                }
                else {
                    std::cout << "Unknown arg" << std::endl;
                }
            default:
                return 1;
        }
    }

    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    int mpirank, mpisize;
    MPI_Comm_rank(comm, &mpirank);
    MPI_Comm_size(comm, &mpisize);

    srand48(mpirank);

    MPI_Status status;

    TSP2DState parallel_state = TSP2DState::from_text_file(filepath);
    Annealer parallel_annealer = Annealer(parallel_state.num_stops(), LOG);

    double min_objective;
    double residual;
    std::vector<long> min_state;
    long size = parallel_state.num_stops();
    double curr_objective = parallel_state.objective();

    long iters = 0;

    // while(parallel_annealer.get_iteration() < MAX_ANNEALER_ITERATIONS && iters < MAX_ITERATIONS) {
    // do {
        // each process searches for a next state
        min_objective = parallel_annealer.anneal(&parallel_state, ANNEALING_STEPS_PER_ITERATION, MAX_ANNEALER_ITERATIONS);
        min_state = parallel_annealer.get_min_state();

        if (mpirank != 0) {
            MPI_Send(&min_objective, 1, MPI_DOUBLE, 0, 999, comm);
            MPI_Send(min_state.data(), size, MPI_LONG, 0, 999, comm);
        }

        if (mpirank == 0) {
            long* recv_min_state = (long*) malloc(size * sizeof(long));
            double recv_min_objective;
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

            std::cout << iters+1 << ": Best tour length = " << min_objective << std::endl;
        }

        //Everyone's min_state should be the network minimum;
        residual = curr_objective - min_objective;
        parallel_state.set_idxs(min_state);
        curr_objective = min_objective;
        iters++;
    } 
    // while (residual > TOLERANCE);

    if (mpirank == 0) {
        std::cout << "Final objective: " << min_objective << std::endl;
    }

    if (mpirank == mpisize - 1) {
        std::cout << "Final objective from last node: " << min_objective << std::endl;
    // }

    MPI_Finalize();
}
