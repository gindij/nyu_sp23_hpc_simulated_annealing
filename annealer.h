// Simulated Annealing - HPC Spring 2023, Stadler
// Brady Edwards & Jack Gindi

#include <vector>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <mpi.h>
#include <omp.h>
#include "tsp2d.h"

enum BetaScaler {
    LOG,
    LOGLOG,
    EXP
};

class Annealer {
    private:
        long iteration;
        long size;
        double beta;
        std::vector<std::vector<double>> t_matrix;
        std::vector<long> min_state;
    public:
        Annealer(long N, BetaScaler b) : iteration(1), size(N), t_matrix(N, std::vector<double>(N, 0)), min_state(N, 0) {
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
            long t = this->size;
            TSP2DTransition proposal = TSP2DTransition(0,0);
        #pragma omp parallel for collapse(2)
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

        void update_t_matrix(TSP2DTransition* trans, TSP2DState* curr_state) {
            long t = this->size;
            TSP2DTransition proposal_row = TSP2DTransition(0,0);
            TSP2DTransition proposal_column = TSP2DTransition(0,0);
            long ii = trans->swap_first;
            long jj = trans->swap_second;
            for (long j = 0; j < t; ++j) {
                proposal_row = TSP2DTransition(ii,j);
                proposal_column = TSP2DTransition(j,ii);
                double dE_row = curr_state->energy_local(proposal_row);
                double dE_column = curr_state->energy_local(proposal_column);
                if (dE_row > 0) {
                    this->t_matrix[ii][j] = 100.*exp(-dE_row * beta)/(t*t); // mulitply by some power of 10 for large matrices?
                } else {
                    this->t_matrix[ii][j] = 100.*1./(t*t);
                    
                }
                if (dE_column > 0) {
                    this->t_matrix[j][ii] = 100.*exp(-dE_column * beta)/(t*t); // mulitply by some power of 10 for large matrices?
                } else {
                    this->t_matrix[j][ii] = 100.*1./(t*t);
                    
                }
            }
            for (long i = 0; i < t; ++i) {
                proposal_row = TSP2DTransition(jj,i);
                proposal_column = TSP2DTransition(i,jj);
                double dE_row = curr_state->energy_local(proposal_row);
                double dE_column = curr_state->energy_local(proposal_column);
                if (dE_row > 0) {
                    this->t_matrix[jj][i] = 100.*exp(-dE_row * beta)/(t*t); // mulitply by some power of 10 for large matrices?
                } else {
                    this->t_matrix[jj][i] = 100.*1./(t*t);
                    
                }
                if (dE_column > 0) {
                    this->t_matrix[i][jj] = 100.*exp(-dE_column * beta)/(t*t); // mulitply by some power of 10 for large matrices?
                } else {
                    this->t_matrix[i][jj] = 100.*1./(t*t);
                    
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
        void select_transition(TSP2DState* curr_state, TSP2DTransition* trans) {
            double prob = drand48() * 100.;
            double prob_sum = 0.0;
            long t = curr_state->num_stops();
            for (int i = 0; i < t; ++i) {
                for (int j = 0; j < t; ++j) {
                    prob_sum += t_matrix[i][j];
                    if (prob < prob_sum) {
                        trans->swap_first = i;
                        trans->swap_second = j;
                        return;
                    }
                }
            }
            trans->swap_first = 0;
            trans->swap_second = 0;
            return;
        }

        // Want ability to continue with current beta? or kill if we get stuck in a minima
        double anneal(TSP2DState* curr_state, long iters_to_run, long max_iters) {
            double curr_objective = curr_state->objective();
            double min_objective = curr_objective;
            double b, e;
            long iters = 0;
            TSP2DTransition trans = TSP2DTransition(0,0);
            while (iters < iters_to_run && this->iteration + iters < max_iters) { //**** temporary cap on iterations ****
                long curr_it = this->iteration + iters;
                if (curr_it % 3 == 0) {
                    this->beta = log(exp(this->beta) + 1);
                }
                if (this->iteration == 0) {
                    this->generate_t_matrix(curr_state);
                } else {
                    this->update_t_matrix(&trans, curr_state);
                }
                this->select_transition(curr_state, &trans);
                curr_state->step(&trans);
                curr_objective = curr_state->objective();
                if (curr_objective < min_objective) {
                    min_objective = curr_objective;
                    for (long i = 0; i < this->size; i++) {
                        this->min_state[i] = curr_state->get_single_idx(i);
                    }
                }
                iters++;
            }
            this->iteration += iters;
            return min_objective;
        }

        std::vector<long> get_min_state() {
            return this->min_state;
        }
};