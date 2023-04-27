#include <iostream>
#include "annealer.cpp"

int main() {

    // create points representing the unit square
    double x[4] = {1., 1., 0., 0.};
    double y[4] = {1., 0., 1., 0.};

    double* x11 = (double*) malloc(4 * sizeof(double));
    double* x12 = (double*) malloc(4 * sizeof(double));
    for (int i = 0; i < 4; i++) {
        x11[i] = x[i];
        x12[i] = x[i];
    }

    double* x21 = (double*) malloc(4 * sizeof(double));
    double* x22 = (double*) malloc(4 * sizeof(double));
    for (int i = 0; i < 4; i++) {
        x21[i] = y[i];
        x22[i] = y[i];
    }

    long N = 4;

    // create states and a transition
    TSP2DState state1 = TSP2DState(N, x11, x21);
    TSP2DState state2 = TSP2DState(N, x12, x22);
    TSP2DTransition transition = TSP2DTransition(2, 3);

    printf("TRANSITION:\n");

    // print the current path and the path length (the objective value)
    state1.display_state();
    state1.display_coords();
    printf("obj = %f\n", state1.objective());

    // use the transition to evolve the state
    state1.step(&transition);

    state1.display_state();
    printf("obj = %f\n", state1.objective());

    printf("ENERGY:\n");

    printf("energy (1->2) = %f\n", state1.energy(&state2));
    printf("energy (2->1) = %f\n", state2.energy(&state1));

    printf("energy local test (1->2) = %f\n", state2.energy_local(transition));

    // Annealer TSP = Annealer(state2.stops(), LOG);
    // TSP.generate_t_matrix(&state2);
    // TSP.display_params();
    // TSP.display_t_matrix();
    // TSP.anneal(&state2);

    printf("FROM FILE:\n");

    TSP2DState from_file = TSP2DState::from_text_file("example_tsp_in.txt");
    TSP2DTransition transition1 = TSP2DTransition(0, 1);
    TSP2DTransition transition2 = TSP2DTransition(71, 49);

    printf("obj = %f\n", from_file.objective());
    from_file.step(&transition1);
    from_file.step(&transition2);
    printf("obj = %f\n", from_file.objective());

    Annealer TSP = Annealer(from_file.stops(), LOG);
    TSP.generate_t_matrix(&from_file);
    TSP.display_params();
    TSP.anneal(&from_file);

    from_file.write_txt("example_tsp_out.txt");
}