# Testing script for the network Sellke Simulation
import networkx as nx
import numpy as np
import numpy.random as npr
from SellkeSimulation.simulation_code import complex_epidemic_simulation

G_test = nx.complete_graph(200)


def fixed_length(para, n): return np.array([5]*n)


test_sim_2 = complex_epidemic_simulation(G=G_test,
                                         beta=0.008,
                                         I_parameters=1.5,
                                         initial_infected=[0, 1, 2, 3, 4],
                                         infection_period_distribution=fixed_length,
                                         time_increment=0.1,
                                         max_iterations=1000)

test_copy = test_sim_2

test_copy_2 = test_sim_2


"""

This test is no longer viable as we now initialise the infection differently now.

def test_initialise_infection():
    #Test that the infection is initialized correctly in the network
    
    epidemic_simulation = complex_epidemic_simulation(G = G_test,
        beta = 0.008,
        I_parameters = 1.5,
        initial_infected = 5)
    npr.seed(1)
    epidemic_simulation.initialise_infection()
    npr.seed(1)
    test_list = np.random.choice(range(200), size = 5, replace = False)
    test_list = list(test_list)
    test_list = sorted(test_list)
    infected = epidemic_simulation.infected_nodes
    assert test_list == infected

We should just check that the numbers are correct and within bounds.
"""


def test_initialise_infection_list():
    test_list = [1, 2, 3, 4, 5]
    epidemic_simulation = complex_epidemic_simulation(G=G_test,
                                                      beta=0.008,
                                                      I_parameters=1.5,
                                                      initial_infected=test_list,
                                                      time_increment=0.1,
                                                      max_iterations=1000)
    assert epidemic_simulation.infected_nodes == test_list


"""
Same with this test. Anything to do with seeds has broken.

def test_infectious_period_generation():
    #We test the generate_infection_periods method to check it is behaving as expected.
    
    npr.seed(1)
    epidemic_simulation = complex_epidemic_simulation(G = G_test,
        beta = 0.008,
        I_parameters = 1.5,
        initial_infected = [1],
        infection_period_distribution = npr.exponential)
    npr.seed(1)
    test_output = npr.exponential(1.5, 200)
    assert all(test_output == epidemic_simulation.infectious_periods)
"""


def test_cumulative_hazards():
    """
    We test the the cumulative hazard for each node is computed correctly.
    Inputs are for a simulation where the infectious period distribution is of fixed length 5.
    Therefore, the expected output is 5 * beta * hazard rate, where beta = 0.008 and the hazard rate is 1 by default.
    """
    sim = complex_epidemic_simulation(G=G_test,
                                      beta=0.008,
                                      I_parameters=1.5,
                                      initial_infected=5,
                                      infection_period_distribution=fixed_length,
                                      time_increment=0.1,
                                      max_iterations=1000)
    sim.calculate_total_emitted_hazard()
    expected_answer = [5*0.008] * 200
    assert expected_answer == sim.lifetime_emitted_hazard


def test_cumulative_exposure():
    """
    Tests that exposure level is being computed correctly.
    To do this we run the update exposure level method until the infection runs their course.
    We do not actually update infection status.
    The infection period is fixed at 5.
    """
    sim = complex_epidemic_simulation(G=G_test,
                                      beta=0.008,
                                      I_parameters=1.5,
                                      initial_infected=[0, 1, 2, 3, 4],
                                      infection_period_distribution=fixed_length,
                                      time_increment=0.1,
                                      max_iterations=1000)

    # We loop the update exposure level method for a bit.
    for _ in range(100):
        sim.updates_exposure_levels()

    # The fixed_length dist returns an infectious period of length 5
    expected_answer = [4*0.008*5]*5
    # There for the uninfected,
    expected_answer.extend([5*0.008*5]*195)
    assert sim.exposure_level == expected_answer


def test_determine_new_infections():
    """
    Tests the update infection status method is correctly updating.
    We set the resistance low for the first 10 nodes (including the first 5 which are initially infected)
    We expect to see the first 10 to be infected, and no more.
    """

    test_sim_2 = complex_epidemic_simulation(G=G_test,
                                             beta=0.008,
                                             I_parameters=1.5,
                                             initial_infected=[0, 1, 2, 3, 4],
                                             infection_period_distribution=fixed_length,
                                             time_increment=0.1,
                                             max_iterations=1000)

    test_copy = test_sim_2
    [test_copy.epi_data[node].update({"Resistance": 0.039})
     for node in test_copy.node_keys[0:10]]
    [test_copy.epi_data[node].update({"Resistance": 100})
     for node in test_copy.node_keys[10:200]]
    test_copy.determine_new_infections()
    assert test_copy.infected_nodes == list(range(10))


def test_iterate_epidemic_successfully():
    """Tests the control structure for iterating epidemics."""
    test_copy = test_sim_2
    [test_copy.epi_data[node].update({"Resistance": 0.039})
     for node in test_copy.node_keys[0:10]]
    [test_copy.epi_data[node].update({"Resistance": 1000})
     for node in test_copy.node_keys[10:200]]
    test_copy.determine_new_infections()
    test_copy.iterate_epidemic()
    assert test_copy.final_size == 10
    assert test_copy.iterations == 61


def test_iterate_epidemic_real_epidemic():
    npr.seed(1)
    [test_copy.epi_data[node].update({"Resistance": 0.039})
     for node in test_copy.node_keys[0:10]]
    [test_copy.epi_data[node].update({"Resistance": 0.21})
     for node in test_copy.node_keys[10:20]]
    test_copy.determine_new_infections()
    test_copy.iterate_epidemic()
    assert test_copy.time > 1


def test_node_list_complex():
    """For some graphs, i.e; a lattice where nodes are named using lists. We cannot refer to them numerically, so we have more advanced logic to handle this.
    """
    G_test_lattice = nx.grid_2d_graph(1, 1)
    my_epidemic = complex_epidemic_simulation(G_test_lattice,
                                              beta=0.008,
                                              I_parameters=1.5,
                                              initial_infected=1,
                                              time_increment=0.1,
                                              max_iterations=1000)
    my_epidemic.iterate_epidemic()
    assert my_epidemic.final_size == 1


def test_complex_node_list_iteration():
    """We check that the iteration can be successfully performed over a complex node list
    """
    G_test_lattice = nx.grid_2d_graph(2, 2)
    my_epidemic = complex_epidemic_simulation(G_test_lattice,
                                              beta=100,
                                              I_parameters=1,
                                              initial_infected=1,
                                              time_increment=0.1,
                                              max_iterations=1000)
    my_epidemic.iterate_epidemic()
    assert my_epidemic.final_size == 4


def test_complex_node_list_iteration_larger_network():
    """We check that the iteration can be successfully performed over a complex node list
    """
    G_test_lattice = nx.grid_2d_graph(10, 10)
    my_epidemic = complex_epidemic_simulation(G_test_lattice,
                                              beta=100,
                                              I_parameters=1,
                                              initial_infected=1,
                                              time_increment=0.1,
                                              max_iterations=1000)
    my_epidemic.iterate_epidemic()
    assert my_epidemic.final_size > 4
    assert my_epidemic.time > 3

# def test_initialise_with_list_node():
#    """We check that the iteration can be successfully performed over a complex node list
#    """
#    G_test_lattice = nx.grid_2d_graph(10,10)
#    my_epidemic = complex_epidemic_simulation(G_test_lattice,beta = 100, I_parameters = 1, initial_infected = (1,1))
#    my_epidemic.iterate_epidemic()
#    assert my_epidemic.final_size > 4
#    assert my_epidemic.time > 3
