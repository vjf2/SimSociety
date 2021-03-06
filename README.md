# SimSociety
Simulate societies to compare null models

**sims.R** creates the simulated dataset 

**sim_obs.R** subsamples a simulated dataset, runs a null model, and compares the output 

**figures.R** plots the figures found on the SimSociety vignette at <http://people.duke.edu/~vjf2/Sim_Society.html>

**animate.R** creates gif of simulation

Input files:

**decent_pp2.RData**: optimized preferred pair network used in example

*SimSociety now uses the Sah et al. 2014 random modular graph generator found at <https://github.com/prathasah/random-modular-network-generator>

Output files:

**sim_results_2000steps.csv**: output of an example 2000 step model with individual locations and groups at each step

**sim_categories_2000steps.csv**: key with assigned categories (random, preference, avoidance) for individuals in results file 
