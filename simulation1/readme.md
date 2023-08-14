### Reproduce simulation data and results

To generate simulation data in each simulation scenario, run the data generation R script in each folder. For example, in folder "scenario1_1", you can run "DataGenerator1.R". The 50 simulation data sets would be generated in the data folder. In each data folder, we gave one example. 

To obtain the simulation results  of different methods, just change the command `data_folder` in "model1.R" and run "model1.R". The R script requires parallel computation with 10 cpu cores. 

Since the computation is intensive, we keep all the results in the results folder in each simulation scenario. 