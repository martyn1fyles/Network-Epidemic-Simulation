from simulation_code import SIR_Selke
import numpy.random as npr
import numpy as np
import pytest

def test_data_gen_simple():
    #We do not specify the distribution, we check that the function draws the correct vector from the exponential distribution
    #the parameter for the infectious periods is set to 1.5
    npr.seed(1)
    simulation = SIR_Selke(200, 0.008, 1.5, 5)
    test_var = simulation.generate_infection_periods()

    npr.seed(1)
    test_var2 = npr.exponential(1.5, 200)
    assert all(test_var2 == test_var)
    

def test_data_gen_non_markovian_1_par():
    npr.seed(1)
    simulation = SIR_Selke(200, 0.008, 0.9, 5, infection_period_distribution = npr.geometric)
    test_var = simulation.generate_infection_periods()

    npr.seed(1)
    test_var2 = npr.geometric(0.9, 200)
    assert all(test_var2 == test_var)

def test_data_gen_non_markovian_int_parameter():
    npr.seed(1)
    simulation = SIR_Selke(200, 0.008, 1, 5, infection_period_distribution = npr.exponential)
    test_var = simulation.generate_infection_periods()

    npr.seed(1)
    test_var2 = npr.exponential(1, 200)
    assert all(test_var2 == test_var)

#def test_error_unless_pos_distribution():
    #simulation = SIR_Selke(200, 0.008, 1, 5, non_markovian_dist = npr.normal)

def test_final_size_calculation():
    #Tests the calculation of the final size epidemic against the manually calculated value
    simulation = SIR_Selke(200, 0.008, 1, 5)

    #Generated the infection periods
    npr.seed(1)
    T = simulation.generate_infection_periods()
    #Generate 200 infection resistances
    #Set the first 5 to 0 to represent the initially infected
    #Sort to determine order of infection
    #Determine cumulative FOI
    res_to_inf = npr.exponential(1,200)
    res_to_inf[:5] = 0
    res_to_inf = np.sort(res_to_inf)
    total_force_of_infection = 0.008*np.cumsum(T)

    #Calcualte the first time that the total force of infection is exceeded by resistance

    counter = 0
    for i in range(200):
        if res_to_inf[i] < total_force_of_infection[i]:
            counter = counter + 1
        else:
            break
    counter

    npr.seed(1)
    test_var = simulation.compute_final_size()
    assert test_var == counter

def test_repeated_sim():
    #Tests that the results from performing repeated iterations of the final size calculation are as expected
    
    npr.seed(1)
    simulation = SIR_Selke(200, 0.008, 1, 5)
    simulation.sim_final_size(2)
