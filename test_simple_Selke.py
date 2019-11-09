from simulation_code import SIR_Selke
from simulation_code import hazard_class
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
    

    npr.seed(1)
    test_var = simulation.compute_final_size()
    assert test_var == counter

def test_repeated_sim():
    #Tests that the results from performing repeated iterations of the final size calculation are as expected
    
    npr.seed(1)
    simulation = SIR_Selke(200, 0.008, 1, 5)
    simulation.sim_final_size(2)

def test_hazard_well_behaved():
    '''
    Tests for out of bounds behaviour
    negative time value returns
    valid time value returns correct value
    time value after the end time point returns 0 
    '''
    def my_hazard_fn(t): return 4*t
    my_hazard = hazard_class(hazard_function= my_hazard_fn)
    t_end = 10
    assert my_hazard.hazard(-1,t_end) == 0
    assert my_hazard.hazard(1,t_end) == 4
    assert my_hazard.hazard(11,t_end) == 0

def test_negative_hazard_func():
    '''
    Checks that the hazard class maps negative values to 0.
    Negative hazard does not make sense
    '''
    def my_hazard_fn(t): return -t
    my_hazard = hazard_class(hazard_function= my_hazard_fn)
    assert my_hazard.hazard(1,10) == 0

def test_total_hazard_function():
    '''
    The total hazard function is the sum of several hazard functions of individuals
    with all times starting at 0.
    The required inputs are:
    list of lengths of infection periods

    '''
    
    infection_lengths = [1,2]
    def my_hazard_fn(t): return 4*t
    my_hazard = hazard_class(hazard_function = my_hazard_fn)

    #tests that it returns 0
    assert my_hazard.total_hazard_function(-1, infection_lengths) == 0

    #tests that both hazards are active
    t_in = 0.5
    assert my_hazard.total_hazard_function(0.5, infection_lengths) == 4

    #Tests the first cut off has happened
    assert my_hazard.total_hazard_function(1.5, infection_lengths) == 6

    #Tests that the cut off has happened at t = 2
    assert my_hazard.total_hazard_function(2.5, infection_lengths) == 0

def test_hazard_integral():
    '''tests the hazard integrator against a predetermined value'''

    def my_hazard(t): return 4*t
    my_hazard = hazard_class(hazard_function = my_hazard)
    assert my_hazard.integrate_hazard(10) == 200

def test_final_size_using_hazard_function():
    '''
    Tests whether the final size of the epidemic is the same when using the
    hazard rate integration for exponential rate, as it is when we precaculate the result
    '''
    npr.seed(1)

    def my_hazard(t): return 1
    simulation = SIR_Selke(200, 0.008, 1, 5, hazard_rate = my_hazard)
    test_var_1 = simulation.compute_final_size()

    npr.seed(1)
    simulation = SIR_Selke(200, 0.008, 1, 5)
    test_var_2 = simulation.compute_final_size()

    assert test_var_1 == test_var_2


