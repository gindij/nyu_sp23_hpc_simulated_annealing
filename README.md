# HPC Final Project: Simulated Annealing
Authors: Jack Gindi and Brady Edwards

## 2D TSP

### Generating data
You can generate data for the traveling salesman problem using `generate_tsp_data.py`. An example command would look like `python generate_tsp_data.py --n-samples 100 --dim 2 --spread 2.5 --plot`. To see the other available command line arguments, run `python generate_tsp_data.py --help`.

### States and transitions
The state and transition objects for the 2d TSP problem can be found in `tsp2d.cpp`. The code includes the ability to read the output described in the previous section. For examples of how to work with it in C++, see `main()` in `tsp2d.cpp`.