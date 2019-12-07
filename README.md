## Network-Epidemic-Simulation
Network-Epidemic-Simulation is a project which aims to allow for easy modelling of complex infectious epidemics on networks and the effectiveness of interventions.

By complex, we aim to provide support for the following features:
1) Non-markovian epidemics, where the time-until infection may follow any other distribution
2) Static and Dynamic Networks. For the research which motivated the development of this software, a key feature was that the network must evolve over time, so that is supported and some basic visualisation tools provided.
3) User-defined treatment scenarios. As the key goal of epidemiological modelling is to evaluate the effectiveness of different treatment scenarios, the system for implementing these is quick and customisable. This system can also be dependent on events occurring in the dynamic network. For example, a node may enter a hospital and receive treatment, which updating it's contact network to include those in the hospital
4) Support for many different categories of epidemics. Currently we support SI, SIS and SIR. However, we have plans to extend to a versatile class system, allowing for any possible epidemic model.

[![Build Status](https://travis-ci.org/martyn1fyles/Sellke-Networks.svg?branch=master)](https://travis-ci.org/martyn1fyles/Sellke-Networks)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/martyn1fyles/Sellke-Networks/master)

# Current Features
1) Time-varying infection profiles, leading to non-markovian epidemics.
2) Static networks and dynamic networks
3) A dynamic networks event system, where events in the dynamic network can effect the epidemic
4) The most common epidemic classes
5) User customisable treatment

# In progress
1) Support for birth and death processes on networks
2) Optimisation and implementing parallelisation
3) User-friendly tweaks
4) Better visualisation
5) More tutorials and documentation
6) A class system which allows for user-defined states and behaviours, e.g: "In Treatment", "In Recovery", "Exposed".

# Notebooks
We have several notebooks demonstrating the use of the package, these can be viewed through Binder so that you may get an idea of the use of the package before downloading yourself.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/martyn1fyles/Sellke-Networks/master)

# Testing
Our code is unit tested with a high coverage so that you can be assured of the validity of the results of the simulations.
