# Testing script for the network Sellke Simulation
import networkx as nx
import numpy as np
import numpy.random as npr
from SellkeSimulation.EpidemicSimulation import epidemic_data

G_complete = nx.complete_graph(10)
G_lattice = nx.grid_2d_graph(5, 5)


def test_update_exposure_level():
    """We call the method which increments a nodes exposure and check that the correct entry in the dictionary is updated"""
    my_data_structure = epidemic_data(G_complete,
                                     initial_infected=[1],
                                     pre_gen_data=100)
    my_data_structure.update_exposure_level(2, 3)
    assert my_data_structure.epi_data[2]["Exposure Level"] == 3
    assert my_data_structure.epi_data[1]["Exposure Level"] == 0
    assert my_data_structure.epi_data[3]["Exposure Level"] == 0


def test_1_par_distribution():

    def dist_1_par(par_1, n): return np.array([5]*n)*par_1
    my_data = epidemic_data(G_complete,
                            initial_infected=[1],
                            pre_gen_data=100,
                            infection_period_distribution=dist_1_par,
                            infection_period_parameters=4)
    assert my_data.epi_data[1]["Infection Period"] == 20


def test_2_par_distribution():

    def dist_2_par(par_1, par_2, n): return np.array([5]*n)*par_1*par_2
    my_data = epidemic_data(G_complete,
                            initial_infected=[1],
                            pre_gen_data=100,
                            infection_period_distribution=dist_2_par,
                            infection_period_parameters=[1, 2])
    assert my_data.epi_data[1]["Infection Period"] == 10


def test_3_par_distribution():

    def dist_3_par(par_1, par_2, par_3, n): return np.array(
        [5]*n)*par_1*par_2*par_3
    my_data = epidemic_data(G_complete,
                            initial_infected=[1],
                            pre_gen_data=100,
                            infection_period_distribution=dist_3_par,
                            infection_period_parameters=[1, 2, 3])
    assert my_data.epi_data[1]["Infection Period"] == 30


def test_initialise_resistances():
    """Cant think of a better way to check that the resistances are exponential(1) w/out resorting to monte carlo.
    We also check that it's not the same in memory.
    """
    my_data = epidemic_data(G_complete,
                            initial_infected=[0],
                            pre_gen_data=100)
    assert all([my_data.epi_data[node]["Resistance"]
                > 0 for node in my_data.epi_data if my_data.epi_data[node]["Infection Stage"]== "Susceptible"])
    assert my_data.epi_data[1]["Resistance"] != my_data.epi_data[2]["Resistance"]


def test_update_infection_stage():
    """We test that when we update an infection stage that:
    1) The update goes through correctly
    2) Only the intended node is getting updated
    3) The status change is correctly added to the log """
    my_data = epidemic_data(G_complete,
                            initial_infected=[1],
                            pre_gen_data=100)
    my_data.update_infection_stage([5, 6, 7, 8, 9], "Stage 1", 1)
    my_data.update_infection_stage([8, 9], "Stage 2", 2)

    assert my_data.epi_data[5]["Infection Stage"] == "Stage 1"
    assert my_data.epi_data[5]["Infection Stage Started"] == 1

    assert my_data.epi_data[8]["Infection Stage"] == "Stage 2"
    assert my_data.epi_data[8]["Infection Stage Started"] == 2

    assert my_data.epi_data[5]["History"]["Infection Stage Log"] == ["Susceptible", "Stage 1"]
    assert my_data.epi_data[8]["History"]["Infection Stage Log"] == ["Susceptible", "Stage 1", "Stage 2"]

    assert my_data.epi_data[5]["History"]["Infection Stage Times"] == [0, 1]
    assert my_data.epi_data[8]["History"]["Infection Stage Times"] == [0, 1, 2]


def test_update_exposure_level_tuple():
    """We call the method which increments a nodes exposure and check that the correct entry in the dictionary is updated and only the correct entry"""
    my_data_structure = epidemic_data(G_lattice,
                                     initial_infected=[(1, 1)],
                                     pre_gen_data=100)
    my_data_structure.update_exposure_level((2, 2), 3)
    assert my_data_structure.epi_data[(2, 2)]["Exposure Level"] == 3
    assert my_data_structure.epi_data[(3, 3)]["Exposure Level"] == 0


def test_1_par_distribution_tuple():

    def dist_1_par(par_1, n): return np.array([5]*n)*par_1
    my_data = epidemic_data(G_lattice,
                            initial_infected=[(1, 1)],
                            pre_gen_data=100,
                            infection_period_distribution=dist_1_par,
                            infection_period_parameters=4)
    assert my_data.epi_data[(1, 1)]["Infection Period"] == 20


def test_2_par_distribution_tuple():

    def dist_2_par(par_1, par_2, n): return np.array([5]*n)*par_1*par_2
    my_data = epidemic_data(G_lattice,
                            initial_infected=[(1, 1)],
                            pre_gen_data=100,
                            infection_period_distribution=dist_2_par,
                            infection_period_parameters=[1, 2])
    assert my_data.epi_data[(1, 1)]["Infection Period"] == 10


def test_3_par_distribution_tuple():

    def dist_3_par(par_1, par_2, par_3, n): return np.array(
        [5]*n)*par_1*par_2*par_3
    my_data = epidemic_data(G_lattice,
                            initial_infected=[(1, 1)],
                            pre_gen_data=100,
                            infection_period_distribution=dist_3_par,
                            infection_period_parameters=[1, 2, 3])
    assert my_data.epi_data[(1, 1)]["Infection Period"] == 30


def test_initialise_resistances_tuple():
    """Cant think of a better way to check that the resistances are exponential(1) w/out resorting to monte carlo.
    We also check that it's not the same in memory.
    """
    my_data = epidemic_data(G_lattice,
                            pre_gen_data=100,
                            initial_infected=[(1, 1)])
    assert all([my_data.epi_data[node]["Resistance"]
                > 0 for node in my_data.epi_data if my_data.epi_data[node]["Infection Stage"]== "Susceptible"] )
    assert my_data.epi_data[(0, 0)]["Resistance"] != my_data.epi_data[
        (2, 2)]["Resistance"]


def test_update_infection_stage_tuple():
    """We test that when we update an infection stage that:
    1) The update goes through correctly
    2) Only the intended node is getting updated
    3) The status change is correctly added to the log """
    my_data = epidemic_data(G_lattice,
                            initial_infected=[(0, 0)],
                            pre_gen_data=100)
    my_data.update_infection_stage([(1, 1), (2, 2)], "Stage 1", 1)
    my_data.update_infection_stage([(2, 2)], "Stage 2", 2)

    assert my_data.epi_data[(1, 1)]["Infection Stage"] == "Stage 1"
    assert my_data.epi_data[(1, 1)]["Infection Stage Started"] == 1

    assert my_data.epi_data[(2, 2)]["Infection Stage"] == "Stage 2"
    assert my_data.epi_data[(2, 2)]["Infection Stage Started"] == 2

    assert my_data.epi_data[(1, 1)]["History"]["Infection Stage Log"] == [
        "Susceptible", "Stage 1"]
    assert my_data.epi_data[(2, 2)]["History"]["Infection Stage Log"] == [
        "Susceptible", "Stage 1", "Stage 2"]

    assert my_data.epi_data[(1, 1)]["History"]["Infection Stage Times"] == [
        0, 1]
    assert my_data.epi_data[(2, 2)]["History"]["Infection Stage Times"] == [
        0, 1, 2]

def test_data_pre_gen():

    def dist_1_par(par_1, n): return np.array([5]*n)*par_1

    my_data = epidemic_data(G_lattice,
                            pre_gen_data=100,
                            initial_infected=[(1, 1)],
                            infection_period_distribution=dist_1_par,
                            infection_period_parameters=1)
    
    test_resistance = my_data.epi_data[(1,1)]["Pre-generated Data"]["Resistance"]
    test_infection_periods = my_data.epi_data[(1,1)]["Pre-generated Data"]["Infection Period"]

    assert len(test_resistance) == 100
    assert len(test_infection_periods) == 100

    assert all([test_infection_periods[i] == 5 for i in range(100)]) == True
