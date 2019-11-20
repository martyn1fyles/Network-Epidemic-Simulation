
#Testing script for epidemic_data class methods
import networkx as nx
import numpy as np
import numpy.random as npr
from SellkeSimulation.EpidemicSimulation import epidemic_data

def test_initialise_data_structure():
    """Tests that the output is as expected for two types of networks.
    """
    G_test = nx.complete_graph(10)
    my_data = epidemic_data(G_test, initial_infected = [0], pre_gen_data = 100)
    assert 0 in my_data.epi_data
    assert my_data.epi_data[1]["Infection Stage"] == "Susceptible"

def test_update_infection_stage():
    """We test that the dictionary updater updates current values correctly, we are not checking the history yet
    """
    G_test = nx.complete_graph(10)
    my_data = epidemic_data(G_test, initial_infected = [1], pre_gen_data = 100)
    my_data.update_infection_stage([0], "Infected", 10)
    assert my_data.epi_data[0]["Infection Stage"] == "Infected"
    assert my_data.epi_data[0]["Infection Stage Started"] == 10
    assert my_data.epi_data[2]["Infection Stage"] == "Susceptible"
    assert my_data.epi_data[2]["Infection Stage Started"] == 0


def test_update_infection_stage_list_input():
    """We test that the dictionary updater updates the correct values
    and only the correct values.
    """
    G_test = nx.complete_graph(10)
    my_data = epidemic_data(G_test, initial_infected = [0], pre_gen_data = 100)
    my_data.update_infection_stage([1,2], "Infected", 10)
    assert my_data.epi_data[1]["Infection Stage"] == "Infected"
    assert my_data.epi_data[1]["Infection Stage Started"] == 10
    assert my_data.epi_data[2]["Infection Stage"] == "Infected"
    assert my_data.epi_data[2]["Infection Stage Started"] == 10

def test_update_infection_stage_2_updates():
    """We test that the dictionary updater updates the correct values
    and only the correct values.
    """
    G_test = nx.complete_graph(10)
    my_data = epidemic_data(G_test, initial_infected = [1], pre_gen_data = 100)
    my_data.update_infection_stage([0], "Infected", 10)
    my_data.update_infection_stage([0], "Susceptible", 20)
    assert my_data.epi_data[0]["Infection Stage"] == "Susceptible"
    assert my_data.epi_data[0]["History"]["Infection Stage Times"] == [0,10,20]
    assert my_data.epi_data[2]["Infection Stage"] == "Susceptible"
    assert my_data.epi_data[2]["History"]["Infection Stage Times"] == [0]



def test_initialise_data_structure_tuple_nodes():
    """Tests that the output is as expected for two types of networks.
    """
    G_test = nx.grid_2d_graph(2,2)
    my_data = epidemic_data(G_test, initial_infected = [(0,0)], pre_gen_data = 100)
    assert (1,1) in my_data.epi_data
    assert my_data.epi_data[(1,1)]["Infection Stage"] == "Susceptible"

def test_update_infection_stage_tuple_nodes():
    """We test that the dictionary updater updates current values correctly, we are not checking the history yet
    """
    G_test = nx.grid_2d_graph(2,2)
    my_data = epidemic_data(G_test, initial_infected = [(1,1)], pre_gen_data = 100)
    my_data.update_infection_stage([(0,0)], "Infected", 10)
    assert my_data.epi_data[(0,0)]["Infection Stage"] == "Infected"
    assert my_data.epi_data[(0,0)]["Infection Stage Started"] == 10
    assert my_data.epi_data[(0,1)]["Infection Stage"] == "Susceptible"
    assert my_data.epi_data[(0,1)]["Infection Stage Started"] == 0


def test_update_infection_stage_list_input_tuple_nodes():
    """We test that the dictionary updater updates the correct values
    and only the correct values.
    """
    G_test = nx.grid_2d_graph(2,2)
    my_data = epidemic_data(G_test, initial_infected = [(1,1)], pre_gen_data = 100)
    my_data.update_infection_stage([(0,0),(0,1)], "Infected", 10)
    assert my_data.epi_data[(0,0)]["Infection Stage"] == "Infected"
    assert my_data.epi_data[(0,0)]["Infection Stage Started"] == 10
    assert my_data.epi_data[(0,1)]["Infection Stage"] == "Infected"
    assert my_data.epi_data[(0,1)]["Infection Stage Started"] == 10

def test_update_infection_stage_2_updates_tuple_nodes():
    """We test that the dictionary updater updates the correct values
    and only the correct values.
    """
    G_test = nx.grid_2d_graph(2,2)
    my_data = epidemic_data(G_test, initial_infected = [(1,1)], pre_gen_data = 100)
    my_data.update_infection_stage([(0,0)], "Infected", 10)
    my_data.update_infection_stage([(0,0)], "Susceptible", 20)
    assert my_data.epi_data[(0,0)]["Infection Stage"] == "Susceptible"
    assert my_data.epi_data[(0,0)]["History"]["Infection Stage Times"] == [0,10,20]
    assert my_data.epi_data[(0,1)]["Infection Stage"] == "Susceptible"
    assert my_data.epi_data[(0,1)]["History"]["Infection Stage Times"] == [0]


def test_initialise_data_structure_infected_as_int():
    """Tests that the output is as expected for two types of networks.
    """
    G_test = nx.grid_2d_graph(2,2)
    my_data = epidemic_data(G_test, initial_infected = 4, pre_gen_data = 100)
    
    assert my_data.epi_data[(1,1)]["Infection Stage"] == "Infected"
    assert my_data.epi_data[(0,0)]["Infection Stage"] == "Infected"

def test_resistances_not_same_element():
    """Tests that when we generate the node resistance we get a different quantity for each node instead of the same value for each node.
    """
    G_test = nx.grid_2d_graph(2,2)
    my_data = epidemic_data(G_test, initial_infected = 1, pre_gen_data = 100)
    assert my_data.epi_data[(0,0)]["Resistance"] != my_data.epi_data[(0,1)]["Resistance"] 



