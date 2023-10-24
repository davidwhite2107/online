# Online learning algorithms for portfolio selection
## Description

This repository contains MATLAB code which implements the Universal Portfolios algorithm by Cover, 
the ANTICOR algorithm by Borodin et al. and the nearest neighbour algorithm by Györfi et al. 
The data on which the algorithms are implemented is from the NYSE and can be found at http://www.cs.bme.hu/~oti/portfolio/data.html. 
The original data is included in the data folder as well as the files nyse.mat and nysemerged.mat which have the data in a usable format.
The code also includes measures of backtest overfitting introduced by Bailey et al.

## Dependencies
The files nyse_all_algorithms.m, nyse_all_algorithms_02.m, nysemerged_all_algorithms_01.m, nysemerged_all_algorithms_02.m
in the scripts folder
all depend on results from other script files in the scripts folder, which must be run first.
## Authors
David White and Tim Gebbie 

