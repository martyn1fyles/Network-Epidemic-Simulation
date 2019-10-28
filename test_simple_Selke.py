from simulation_code import SIR_Selke
import numpy.random as npr
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
    simulation = SIR_Selke(200, 0.008, 0.9, 5, non_markovian_dist = npr.geometric)
    test_var = simulation.generate_infection_periods()

    npr.seed(1)
    test_var2 = npr.geometric(0.9, 200)
    assert all(test_var2 == test_var)

def test_data_gen_non_markovian_int_parameter():
    npr.seed(1)
    simulation = SIR_Selke(200, 0.008, 1, 5, non_markovian_dist = npr.exponential)
    test_var = simulation.generate_infection_periods()

    npr.seed(1)
    test_var2 = npr.exponential(1, 200)
    assert all(test_var2 == test_var)

def test_error_unless_pos_distribution():
    simulation = SIR_Selke(200, 0.008, 1, 5, non_markovian_dist = npr.normal)

def test_final_size_calculation():
    #Tests the the final size of the epidemic against a pre-determined value
    simulation = SIR_Selke(200, 0.008, 1, 5)
    npr.seed(1)
    assert simulation.gen_final_size() == 19 
