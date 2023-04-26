#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>

std::vector<std::string> split_line(std::string line) {

    std::istringstream iss(line);
    std::vector<std::string> words;
    std::string word;
    while (iss >> word) {
        words.push_back(word);
    }
    return words;
}

class TSP2DTransition {
    public:
        TSP2DTransition(long f, long s): swap_first(f), swap_second(s) {};
        long swap_first;
        long swap_second;
};

class TSP2DState {
    private:
        long N;
        double* x1;
        double* x2;
        long* idxs;
    public:

        TSP2DState(long N, double* x1, double* x2) : N(N), x1(x1), x2(x2) {
            this->idxs = (long*) malloc(sizeof(long) * N);
            for (long i = 0; i < N; i++) this->idxs[i] = i;
        }

        TSP2DState(long N) : TSP2DState(N, nullptr, nullptr) {}

        /**
         * Reads a state from a text file. A state consists of
         * two arrays of points and an ordering of nodes to visit. When
         * the state is updated, only the indices are moved around.
         */
        static TSP2DState from_text_file(std::string path) {

            std::ifstream file(path);
            if (file.is_open()) {
                std::string line;
                std::getline(file, line);
                std::vector<std::string> nums = split_line(line);

                long N = std::stol(nums[0]);

                double* x1 = (double*) malloc(N * sizeof(double));
                double* x2 = (double*) malloc(N * sizeof(double));

                for (long i = 0; i < N; i++) {
                    std::getline(file, line);
                    std::vector<std::string> pts = split_line(line);
                    x1[i] = std::stof(pts[0]);
                    x2[i] = std::stof(pts[1]);
                }

                file.close();

                return TSP2DState(N, x1, x2);
            } else {
                std::cout << "Unable to open file" << std::endl;
                abort();
            }
        }

        /**
         * Destroys allocated memory.
         */
        ~TSP2DState() {
            free(this->idxs);
            free(this->x1);
            free(this->x2);
        }

        /**
         * Returns the current tour structure.
         */
        long* get_idxs() {
            return this->idxs;
        }

        /**
         * Steps the state forward given two indices (cities) to swap.
         */
        void step(TSP2DTransition* t) {
            long tmp = this->idxs[t->swap_first];
            this->idxs[t->swap_first] = this->idxs[t->swap_second];
            this->idxs[t->swap_second] = tmp;
        }

        /**
         * The energy for this problem is the objective of the calling state minus the objective
         * of the other state. (For this problem, the objective is total Euclidean distance.)
         */
        double energy(TSP2DState* other) {
            return this->objective() - other->objective();
        }

        /**
         * changing the energy calculator to only calculate the local change.  saves having to rebuild 
         * a bunch of state objects while building transition matrix
        */
        double energy_local(TSP2DTransition* p) {
            double curr_energy = this->objective();
            double new_energy = curr_energy;
            long first_curr = p->swap_first;
            long first_next = this->idxs[(first_curr + 1)%N];
            long first_prev = this->idxs[(first_curr + N - 1)%N];
            long sec_curr = p->swap_second;
            long sec_next = this->idxs[(sec_curr + 1)%N];
            long sec_prev = this->idxs[(sec_curr + N - 1)%N];
            new_energy += (
                sqrt((x1[first_next]-x1[sec_curr])*(x1[first_next]-x1[sec_curr])) +
                sqrt((x1[first_prev]-x1[sec_curr])*(x1[first_prev]-x1[sec_curr])) +
                sqrt((x1[sec_next]-x1[first_curr])*(x1[sec_next]-x1[first_curr])) + 
                sqrt((x1[sec_prev]-x1[first_curr])*(x1[sec_prev]-x1[first_curr])) -
                sqrt((x1[first_next]-x1[first_curr])*(x1[first_next]-x1[first_curr])) -
                sqrt((x1[first_prev]-x1[first_curr])*(x1[first_prev]-x1[first_curr])) -
                sqrt((x1[sec_next]-x1[sec_curr])*(x1[sec_next]-x1[sec_curr])) -
                sqrt((x1[sec_prev]-x1[sec_curr])*(x1[sec_prev]-x1[sec_curr]))
            );
            return curr_energy - new_energy;
        }

        /**
         * Calculates the total Euclidean distance of a tour.
         */
        double objective() {
            double total_dist = 0.0;
            for (long i = 0; i < N; i++) {
                long i_curr = this->idxs[i];
                // wrap around to include the transition back to where the tour started
                long i_next = this->idxs[(i+1)%N];
                double x1_max = std::max(this->x1[i_curr], this->x1[i_next]);
                double x1_min = std::min(this->x1[i_curr], this->x1[i_next]);
                double x2_max = std::max(this->x2[i_curr], this->x2[i_next]);
                double x2_min = std::min(this->x2[i_curr], this->x2[i_next]);

                total_dist += sqrt(
                    (x1_max - x1_min) * (x1_max - x1_min) +
                    (x2_max - x2_min) * (x2_max - x2_min)
                );
            }
            return total_dist;
        }

        /**
         * return number of stops
        */
       long stops() {
            return this->N;
       }

        /**
         * Display a state by printing to the console.
         */
        void display_state() {
            for (int i = 0; i < this->N + 1; i++)
                printf("%ld ", this->get_idxs()[i%N]);
            printf("\n");
        }

        /**
         * Write tour to file pairs of points.
         */
        void write_txt(std::string path) {
            std::ofstream file(path);
            if (file.is_open()) {
                for (long i = 0; i < this->N; i++) {
                    file << std::setprecision(5)
                         << this->x1[i]
                         << " "
                         << this->x2[i]
                         << "\n";
                }
                file.close();
            } else {
                std::cout << "Unable to open file." << std::endl;
                abort();
            }
        }
};

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
    printf("obj = %f\n", state1.objective());

    // use the transition to evolve the state
    state1.step(&transition);

    state1.display_state();
    printf("obj = %f\n", state1.objective());

    printf("ENERGY:\n");

    printf("energy (1->2) = %f\n", state1.energy(&state2));
    printf("energy (2->1) = %f\n", state2.energy(&state1));

    printf("FROM FILE:\n");

    TSP2DState from_file = TSP2DState::from_text_file("example_tsp_in.txt");
    TSP2DTransition transition1 = TSP2DTransition(0, 1);
    TSP2DTransition transition2 = TSP2DTransition(71, 49);

    printf("obj = %f\n", from_file.objective());
    from_file.step(&transition1);
    from_file.step(&transition2);
    printf("obj = %f\n", from_file.objective());

    from_file.write_txt("example_tsp_out.txt");
}