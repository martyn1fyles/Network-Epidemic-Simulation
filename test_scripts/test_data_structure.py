#Testing script for the network Sellke Simulation
import networkx as nx
import numpy as np
import numpy.random as npr
from SellkeSimulation.EpidemicSimulation import epidemic_data

G_complete = nx.complete_graph(10)
G_lattice = nx.grid_2d_graph(5,5)

def test_update_exposure_level():
    """We call the method which increments a nodes exposure and check that the correct entry in the dictionary is updated"""
    my_data_structure =  epidemic_data(G_complete, [1])
    my_data_structure.update_exposure_level(2, 3)
    assert my_data_structure.epi_data[2]["Exposure Level"] == 3