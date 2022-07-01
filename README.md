# ShellGrowth
A repository for the various matlab scripts that run the Shell Growth model published as Green et al. 2021 (Published in JGR-Planets)

All files in this repository must be saved to the same directory in order for the model to run.

convectstagnant.m is the hub script that defines parameters and runs the ODE, and is what should be edited and run to collect additional data. Inside convectstagnant.m are several code modules used to produce the results and figures seen in Green et al. 2021. Also included in this repository (but ultimately not used in publication) is the ability to toggle pressure-dependent melting temperature (via StagnantLidODEv2) for large icy satellites such as Ganymede. A loose framework for simulating convective initiation was also constructed, but ultimately not implemented.
