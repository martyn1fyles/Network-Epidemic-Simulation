#Testing script for the hazard class

import networkx as nx
import numpy as np
import numpy.random as npr
from SellkeSimulation.Simulation import hazard_class
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