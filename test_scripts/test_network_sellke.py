#Testing script for the network Sellke Simulation
import networkx as nx
import numpy as np
import numpy.random as npr
from SellkeSimulation import simulation_code as s


def test_the_testers():
    """
    Does the python package even work? Can we find a single goddamn test? Am I going crazy?
    """
    test_var = 2 + 4
    assert test_var == 6

G_test = nx.complete_graph(200)


def test_initialise_infection():
    #Test that the infection is initialised correctly in the network
    npr.seed(1)
    simulation_class = s.sir_network_sellke_simple(G = G_test,
        beta = 0.008,
        I_parameters = 1.5,
        infected_started = 5)
    simulation_class.initialise_infection()

    npr.seed(1)
    test_list = np.random.choice(range(200), size = 5, replace = False)
    assert all(simulation_class.infected_nodes == test_list)

    test_list = [1,2,3,4,5]
    simulation_class = s.sir_network_sellke_simple(G = G_test,
        beta = 0.008,
        I_parameters = 1.5,
        infected_started = test_list)
    simulation_class.initialise_infection()
    #assert all(simulation_class.infected_nodes == test_list)
