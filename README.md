# Solver for small stochastic differential equations

This is a small library for integration of particle diffusion models based on stochastic differential equations (SDE). 
In particular, the implementation allows Langevin equations with adaptable interfaces for the determinisitic and stochastic forces. 
The `System.h/.cpp` class is the base for the definition of the SDE of Langevin type and handles the numerical integration. 
Particular force terms can be implemented by the user via class inheritance.
For instance, the Force class in `Force.h/.cpp` can be used to implement external and pair forces.
The integration algorithm is the Ermak-McCammon algorithm which also allows implementation of hydrodynamic interactions
through the `HIcouple.h/cpp`.

## Copyright

Copyright (c) 2016-2020, R. Gregor Weiß

## Citations

The code has been used and modified for SDE modeling in these scientific publications:

1. Solvent Fluctuations Induce Non-Markovian Kinetics in Hydrophobic Pocket-Ligand Binding  
R.G. Weiß, P. Setny, and J. Dzubiella J. Phys. Chem. B 120 (33), 8127-8136 (2016)
2. Volume-Scaled Common Nearest Neighbor Clustering Algorithm with Free-Energy Hierarchy  
R.G. Weiß, B. Ries, S. Wang, S. Riniker arXiv:2009.08757 (2020)