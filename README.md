# ShellGrowth
A repository for the various matlab scripts that run the Shell Growth model published as Green et al. 2020 (currently in press)

convectstagnant.m is the hub script that defines parameters and runs the ODE, and is what should be edited and run to collect additional data. Inside convectstagnant.m are several code modules used to produce the results and figures seen in Green et al. 2020 (in press). Also included in this repository (but ultimately not used in publication) is the ability to toggle pressure-dependent melting temperature (via StagnantLidODEv2) for large icy satellites such as Ganymede. A loose framework for simulating convective initiation was also constructed, but ultimately not implemented.
