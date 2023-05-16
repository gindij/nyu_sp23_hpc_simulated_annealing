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
    std::ofstream vizstream;
    std::string vizpath = "vizdefault.txt";
    long MAX_ITERATIONS = 100;
    long MAX_ANNEALER_ITERATIONS = 10000;
    long ANNEALING_STEPS_PER_ITERATION = 100;
    long TOLERANCE = 20;

    while((c = getopt(argc, argv, "x:n:i:j:t:v:")) != -1) {
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
                TOLERANCE = atol(optarg);
                break;
            case 'v':
                vizpath = optarg;
                break;
            case '?':
                if (optopt == 'x' || optopt == 'n' || optopt == 'i' || optopt == 'j' || optopt == 't' || optopt == 'v') {
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
    long timer = 0;
    std::vector<long> min_state;
    long size = parallel_state.num_stops();
    double global_min = parallel_state.objective();
    std::vector<long> global_state = parallel_state.get_idxs();

    long iters = 0;

    // while(parallel_annealer.get_iteration() < MAX_ANNEALER_ITERATIONS && iters < MAX_ITERATIONS) {
    while (timer < TOLERANCE) {
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
            free(recv_min_state);
            if (min_objective < global_min) {
                global_min = min_objective;
                global_state = min_state;
                timer = 0;
            }

            std::cout << iters+1 << ": Best tour length = " << global_min << std::endl;
        }
        std::cout << "Rank " << mpirank << " minimum = " << min_state << std::endl;

        MPI_Bcast(global_state.data(), size, MPI_LONG, 0, comm);
        MPI_Bcast(&timer, 1, MPI_LONG, 0, comm);


        parallel_state.set_idxs(global_state);

        if (iters % 10 == 0 && mpirank == 0) {
            parallel_state.write_txt(vizstream, vizpath);
        }
        iters++;
        timer++;
    } 

    if (mpirank == 0) {
        std::cout << "Final objective: " << global_min << std::endl;
    }

    if (mpirank == mpisize - 1) {
        std::cout << "Final objective from last node: " << min_objective << std::endl;
    }

    MPI_Finalize();
}
