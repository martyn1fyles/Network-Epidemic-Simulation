from simulation_code import SIR_Selke
import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt

#Initiate a simulation class with
#A population of size 200
#A hazard rate of 0.008, or equivalently waiting times between contacts distributed with exponential(beta) waiting times
#The length of the infectious period is exponential distributed around 1.5
#The initial number of infected is 5
simulation_class = SIR_Selke(N = 200, beta = 0.006, I_parameters = 1.5, infected_started = 5)
final_size = simulation_class.sim_final_size(10000)
print(f"The mean final size of the epidemic was {np.mean(final_size)}")


plt.figure()


#The following code produces a histogram of the final size of the epidemic
plt.subplot(121)
simulation_class.plot_hist()


#Trying this put again but with a geometrically distributed period of infection
simulation_class = SIR_Selke(N = 200, beta = 0.006, I_parameters = 0.9, infected_started = 5,infection_period_distribution=npr.geometric)

final_size = simulation_class.sim_final_size(10000)
print(f"The mean final size of the epidemic was {np.mean(final_size)}")

#The following code produces a histogram of the final size of the epidemic
plt.subplot(122)
simulation_class.plot_hist()
plt.show()