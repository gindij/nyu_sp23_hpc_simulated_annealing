// Simulated Annealing - HPC Spring 2023, Stadler
// Brady Edwards & Jack Gindi

#include <vector>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include "tsp2d.h"

enum BetaScaler {
    LOG,
    LOGLOG,
    EXP
};

class Annealer {
    private:
        long iteration;
        double beta;
        std::vector<std::vector<double>> t_matrix;
        std::vector<long> min_state;
    public:
        Annealer(long N, BetaScaler b) : iteration(1), t_matrix(N, std::vector<double>(N, 0)), min_state(N, 0) {
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

        long get_iteration() {
          return this->iteration;
        }

        void display_params() {
            std::cout << "Beta = " << this->beta << "\n" <<
                         "Iteration = " << this->iteration << "\n" <<
                         std::endl;
        }

        void generate_t_matrix(TSP2DState* curr_state) {
            long t = curr_state->num_stops();
            TSP2DTransition proposal = TSP2DTransition(0,0);
            for (long i = 0; i < t; ++i) {
                for (long j = 0; j < t; ++j) {
                    proposal = TSP2DTransition(i,j);
                    double dE = curr_state->energy_local(proposal);
                    if (dE > 0) {
                        this->t_matrix[i][j] = 100.*exp(-dE * beta)/(t*t); // mulitply by some power of 10 for large matrices?
                    } else {
                        this->t_matrix[i][j] = 100.*1./(t*t);
                    }
                }
            }
        }

        void display_t_matrix() {
            long t = this->t_matrix.size();
            for (long i = 0; i < t; ++i) {
                for (long j = 0; j < t; ++j) {
                    std::cout << t_matrix[i][j] << " ";
                }
                std::cout << std::endl;
            }
        }

        // Future: change to a more general transition class?
        TSP2DTransition* select_transition(TSP2DState* curr_state) {
            double prob = drand48() * 100.;
            double prob_sum = 0.0;
            long t = curr_state->num_stops();
            for (int i = 0; i < t; ++i) {
                for (int j = 0; j < t; ++j) {
                    prob_sum += t_matrix[i][j];
                    if (prob < prob_sum) {
                        return new TSP2DTransition(i,j);
                    }
                }
            }
            return new TSP2DTransition(0,0);
        }

        // Want ability to continue with current beta? or kill if we get stuck in a minima
        double anneal(TSP2DState* curr_state, long iters_to_run, long max_iters) {
            double curr_objective = curr_state->objective();
            double min_objective = curr_objective;
            double residual;
            double tol = 1e-5;
            long iters = 0;
            while (iters < iters_to_run && this->iteration + iters < max_iters) { //**** temporary cap on iterations ****
                long curr_it = this->iteration + iters;
                if (curr_it % 50000 == 0) {
                    int cont;
                    std::cout << "Current minimum = " << min_objective << std::endl;
                    std::cout << "Continue?" << std::endl;
                    std::cin >> cont;
                    if (!cont) break;
                }
                if (curr_it % 3 == 0) {
                    this->beta = log(exp(this->beta) + 1);
                }
                this->generate_t_matrix(curr_state);
                TSP2DTransition* trans = this->select_transition(curr_state);
                // std::cout << "Selected transition: (" << trans->swap_first << ", " << trans->swap_second << ")" << std::endl; 
                curr_state->step(trans);
                curr_objective = curr_state->objective();
                // std::cout << "Current energy = " << curr_objective << std::endl;
                // std::cout << "Current beta = " << this->beta << std::endl;
                if (curr_objective < min_objective) {
                    min_objective = curr_objective;
                    this->min_state = curr_state->get_idxs();
                }
                iters++;
            }
            this->iteration += iters;
            // std::cout << "iterations: " << iters << " - Final minimum found: "<<min_objective<<std::endl;
            return min_objective;
        }

        std::vector<long> get_min_state() {
            return this->min_state;
        }
};
