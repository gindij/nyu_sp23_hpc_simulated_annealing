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
        std::vector<long> get_idxs() {
            std::vector<long> copy = std::vector<long>(this->N, 0);
            for (long i = 0; i < this->N; i++) copy[i] = this->idxs[i];
            return copy;
        }

        long get_single_idx(long i) {
            return this->idxs[i];
        }

        void set_idxs(std::vector<long> new_idxs) {
            for (long i = 0; i < this->N; i++)
                this->idxs[i] = new_idxs[i];
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
        double energy_local(TSP2DTransition p) {
            if (p.swap_first == p.swap_second) return 0.0;
            double curr_energy = this->objective();
            double new_energy = curr_energy;
            long first = std::min(p.swap_first,p.swap_second);
            long second = std::max(p.swap_first,p.swap_second);
            long first_curr = this->idxs[first];
            long first_next = this->idxs[(first + 1)%N];
            long first_prev = this->idxs[(first + N - 1)%N];
            long sec_curr = this->idxs[second];
            long sec_next = this->idxs[(second + 1)%N];
            long sec_prev = this->idxs[(second + N - 1)%N];

            if (second - first > 1) {
                new_energy += (
                    sqrt((x1[first_next]-x1[sec_curr])*(x1[first_next]-x1[sec_curr]) + (x2[first_next]-x2[sec_curr])*(x2[first_next]-x2[sec_curr])) +
                    sqrt((x1[first_prev]-x1[sec_curr])*(x1[first_prev]-x1[sec_curr]) + (x2[first_prev]-x2[sec_curr])*(x2[first_prev]-x2[sec_curr])) +
                    sqrt((x1[sec_next]-x1[first_curr])*(x1[sec_next]-x1[first_curr]) + (x2[sec_next]-x2[first_curr])*(x2[sec_next]-x2[first_curr])) +
                    sqrt((x1[sec_prev]-x1[first_curr])*(x1[sec_prev]-x1[first_curr]) + (x2[sec_prev]-x2[first_curr])*(x2[sec_prev]-x2[first_curr])) -
                    sqrt((x1[first_next]-x1[first_curr])*(x1[first_next]-x1[first_curr]) + (x2[first_next]-x2[first_curr])*(x2[first_next]-x2[first_curr])) -
                    sqrt((x1[first_prev]-x1[first_curr])*(x1[first_prev]-x1[first_curr]) + (x2[first_prev]-x2[first_curr])*(x2[first_prev]-x2[first_curr])) -
                    sqrt((x1[sec_next]-x1[sec_curr])*(x1[sec_next]-x1[sec_curr]) + (x2[sec_next]-x2[sec_curr])*(x2[sec_next]-x2[sec_curr])) -
                    sqrt((x1[sec_prev]-x1[sec_curr])*(x1[sec_prev]-x1[sec_curr]) + (x2[sec_prev]-x2[sec_curr])*(x2[sec_prev]-x2[sec_curr]))
                );
            } else {
                new_energy += (
                    sqrt((x1[first_prev]-x1[sec_curr])*(x1[first_prev]-x1[sec_curr]) + (x2[first_prev]-x2[sec_curr])*(x2[first_prev]-x2[sec_curr])) +
                    sqrt((x1[sec_next]-x1[first_curr])*(x1[sec_next]-x1[first_curr]) + (x2[sec_next]-x2[first_curr])*(x2[sec_next]-x2[first_curr])) -
                    sqrt((x1[first_prev]-x1[first_curr])*(x1[first_prev]-x1[first_curr]) + (x2[first_prev]-x2[first_curr])*(x2[first_prev]-x2[first_curr])) -
                    sqrt((x1[sec_next]-x1[sec_curr])*(x1[sec_next]-x1[sec_curr]) + (x2[sec_next]-x2[sec_curr])*(x2[sec_next]-x2[sec_curr]))
                );
            }
            return new_energy-curr_energy;
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
       long num_stops() {
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
         * Display coordinates of the state.
        */
        void display_coords() {
            for (int i = 0; i < this->N + 1; i++)
                printf("(%f, %f) ", this->x1[this->get_idxs()[i%N]], this->x2[this->get_idxs()[i%N]]);
            printf("\n");
        }

        /**
         * Write tour to file pairs of points.
         */
        void write_txt(std::ofstream &stream, std::string path) {
            stream.open(path, std::ios::app);
            if (stream.is_open()) {
                for (long i = 0; i < this->N; i++) {
                    stream << std::setprecision(5)
                         << this->x1[this->idxs[i]]
                         << " "
                         << this->x2[this->idxs[i]]
                         << "\n";
                }
                stream << std::endl;
                stream.close();
            } else {
                std::cout << "Unable to open file." << std::endl;
                abort();
            }
        }
};
