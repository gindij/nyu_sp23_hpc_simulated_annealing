// Simulated Annealing - HPC Spring 2023, Stadler
// Brady Edwards & Jack Gindi

#include <vector>
#include <math.h>
#include <stdlib.h>
#include "tsp2d.cpp"

enum BetaScaler {
    LOG,
    LOGLOG,
    EXP
};

class Annealer {
    private:
        int iteration;
        double beta;
        std::Vector<std::Vector<double>> *t_matrix;
    public:
        Annealer(int N, BetaScaler b) : iteration(1) {
            this->tmatrix->resize(N, std::Vector<double>(N));
            switch(b) {
                case LOG:
                    this->beta = log(1.01);
                    break;
                case LOGLOG:
                    this->beta = log(log(1.01));
                    break;
                case EXP:
                    this->beta = 1.001;
                    break;
            }
        }

        void generate_t_matrix(TSP2DState curr_state) {
            t = curr_state.N;
            for (int i = 0; i < t; ++i) {
                for (int j = 0; j < t; ++j) {
                    TSP2DTransition proposal = TSP2DTransition(i,j);
                    dE = curr_state->energy_local(proposal);
                    if (dE > 0) {
                        t_matrix[i][j] = exp(-dE * beta)/(t*t); // mulitply by some power of 10 for large matrices?
                    } else {
                        t_matrix[i][j] = 1/(t*t);
                    }
                }
            }
        }

        // Future: change to a more general transition class?
        TSP2DTransition select_transition(TSP2DState curr_state) {
            TSP2DTransition selection;
            double prob = rand();
            double prob_sum = 0.0;
            t = curr_state.N;
            for (int i = 0; i < t; ++i) {
                for (int j = 0; j < t; ++j) {
                    prob_sum += t_matrix[i][j];
                    if (prob < prob_sum) {
                        selection = TSP2DTransition(i,j);
                        return selection;
                    }
                }
            }
        }

        // Want ability to continue with current beta? or kill if we get stuck in a minima
        void anneal(TSP2DState curr_state) {
            double curr_objective = curr_state.objective();
            double min_objective = curr_objective;
            double residual;
            double tol = 1e-5;
            do {

            } while () // condition where minimum doesn't change for x number of iterations?
        }

};