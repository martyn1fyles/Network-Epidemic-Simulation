## Sellke-Simulator
This repository contains a Sellke simulation for epidemics in general populaiton. Currently we only consider the scenario where contact occurs between every member of the population.

The duration of infectious periods may be drawn from arbitrary distributions. We have implemented the Sellke formulation of epidemics. This is equivalent to saying:
Time between contacts is exponentially distributed
If an infectious person contacts a susceptible person, the probability of infection conditional upon the hazard rate of the infectious individual
As such, the time until infection is a thinned process with time varying thinning.

##Project Structure

# Notebooks
We have several notebooks demostrating usage of the simulations.

# Testing
Code is tested using pytets and has a high coverage.