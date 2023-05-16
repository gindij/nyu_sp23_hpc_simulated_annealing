#include <iostream>
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


        if (min_objective < global_min) {
            global_min = min_objective;
            global_state = min_state;
            timer = 0;
        }

        parallel_state.set_idxs(global_state);

        iters++;
        timer++;
    } 
    std::cout << "Final objective: " << global_min << std::endl;
}
