# PCASL-simulation
These scripts simulate a blood proton experiencing adiabatic inversion as in pseudo-continuous arterial spin labeling. The simulations are based on [Maccotta et al., 1999](https://doi.org/10.1002/(SICI)1099-1492(199706/08)10:4/5%3C216::AID-NBM468%3E3.0.CO;2-U) and [Dai et al., 2008](https://doi.org/10.1002/mrm.21790).

* The main simulation script is sim_inversion.m which is a function that takes a gradient shape, RF shape, and a sturcture of simulation parameters as input. This script is a functionalized version of pcasl_simulation.m
* sim_inversion_call.m calls sim_inversion() for a few blood velocities and plots them. You can see that slower spins generally experience less inversion.
* maccotta_simulation.m is the CASL simulation.
