This repository relates to a set of excersies on solving, first, the heat equation and then the Schrödinger equation - in various cases. 
Actually, in the end, it goes back to a version of the heat equation - motivated by the time-independent Schrödinger equation.
All examples are 1D ones, and involves simulations in which the solution is question is updated and displayed dynamically as the simulation is running.
To this end, you may have to tweek the settings in your IDE. In Spyder, ror instance, you should set your Graphics backend to "Authomatic", not "Inline".
You can do so this way
Tools -> Preferences -> IPython console -> Graphics tab (2nd tab from the left in the window) -> Graphics backend: Authomatic (scroll down menu) 

The exersies are included here as a pdf file.

The files are the following:
SolveHeatEq.py     -  Solves the heat equation using the Forward-Euler method
SolveHeatEqCN.py   -  Does the same but also including a solution using the Crank-Nicolson method
SolveTDSEnoPot     -  Solves the time-dependent Scrödinger equation for a free particle with a Gaussian initial wave
SolveTDSEwithPot   -  The same as above, but no the particle is exposed to a potential with a smooth rectangular shape
SolveTDSE_ImTime   -  Uses the same implementation as above to find the ground state via propagation in imaginary time

