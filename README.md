# ElegansIndividuality
## input data files
Input data includes 
 - `well-list.csv`: list of individual worms/wells and their strain and stress condition (Days Starvation, *DS*)
 - `roam-frac-10.csv`: roaming fraction for each worm in each of 10 bins per development stage
 
## generating output
Run `julia ./run.jl` in the project's home dir. Output is generated in the `output` directory.

Tested with julia 1.7.2.
