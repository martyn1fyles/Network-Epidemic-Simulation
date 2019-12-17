#Testing script for the hazard class

import networkx as nx
import numpy as np
import numpy.random as npr
from NetworkEpidemicSimulation.Simulation import hazard_class
from pytest import approx

def test_hazard_increment():
    """We test the hazard increment works with a simple function, x^2 between the values of 2 and 3.
    It should return 1/3(27 - 8)
    """
    hazard = lambda x: x**2
    t_0 = 2
    t_1 = 3
    expected_out = (3**3 - 2**3)/3
    my_hazard = hazard_class(hazard)
    out = my_hazard.increment_hazard(t_0,t_1, end_of_infection_time = 10)
    assert out == approx(expected_out)

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

def test_hazard_integral():
    '''tests the hazard integrator against a predetermined value'''

    def my_hazard(t): return 4*t
    my_hazard = hazard_class(hazard_function = my_hazard)
    assert my_hazard.increment_hazard(0,10,10) == 200