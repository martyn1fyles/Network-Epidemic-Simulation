import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
import networkx as nx
from SellkeSimulation.EpidemicSimulation import epidemic_data




class hazard_class:
    '''
    This class handles all the calculations associated with hazard functions.

    Particular duties include:
    hazard rate evaluation
    hazard rate integration
    producing the cumulative hazard vector
    '''
    
    def __init__(self, hazard_function):
        self.hazard_function = hazard_function

    def hazard(self, t, t_end):
        '''
        A housekeeping variant of the hazard function.
        Retuns 0 if t is less than 0, truncates to 0 after t_end
        If no hazard rate was specified, returns 1 so the system defaults to exponential waiting times.
        '''
        self.t_end = t_end
        if t < 0 or t > t_end:
            return 0
        elif self.hazard_function == None:
            return(1)
        else:
            temp = self.hazard_function(t)
            if temp < 0:
                return 0
            else:
                return temp
    
    def total_hazard_function(self, t, T_endpoints):
        '''
        !Not in use!

        Given a predefined hazard function, an input time, and the lengths of the infectious periods
        The function caculates the rate at which the hazard is emitted at time t, assuming they were all intially infected

        Designed as an input function for the BVP Solver, which doesn't currently work
        '''
        self.T_endpoints = T_endpoints
        total_hazard_rate = 0
        for i in range(len(T_endpoints)):
            total_hazard_rate =  total_hazard_rate + self.hazard(t,T_endpoints[i])
        return total_hazard_rate
    
    def integrate_hazard(self, t_end):
        '''
        Integrate a hazard function
        '''
        f = lambda t: self.hazard(t,t_end)
        integral = spi.quad(f, 0, t_end)
        return integral[0]

class complex_epidemic_simulation(epidemic_data):
    '''
    This is the class we use to simulate Sellke on a network

    I think it's possible to run this by calling the existing sellke code

    If you think about it, when we determine whether a susceptible node becomes infected, we are consider the sellke general population problem
    with 1 susceptible and the number of initial infected equal to the number of infected connected via edge
    You know they got infected if he final size of the epidemic is 1 greater than the initial starting size.

    Alternatively, when a node gets infected, we immediately calculate their total hazard contribution.

    This way, when we calculate whether a susceptible becomes infected we can simply sum over the infected neighbours.

    ::Inputs::
    G a network x graph object
    '''

    def __init__(self, G, beta, I_parameters, initial_infected, hazard_rate = None, infection_period_distribution = None):
        self.G = G
        self.node_keys = list(G.nodes())
        self.beta = beta
        self.inf_starting = initial_infected
        self.I_parameters = I_parameters
        self.hazard_rate = hazard_rate
        self.N = nx.number_of_nodes(self.G)
        self.data_structure = epidemic_data(G, initial_infected, infection_period_distribution, I_parameters)
        self.epi_data = self.data_structure.epi_data
        self.new_infections = self.infected_nodes[:]
        self.calculate_total_emitted_hazard()
        self.time = 0
    
    @property
    def infected_nodes(self):
        """Returns the calls the data structure to return the set of infected
        """
        return [nodes for nodes in self.data_structure.epi_data if self.data_structure.epi_data[nodes]["Infection Stage"] == "Infected"]

    @property
    def susceptible_nodes(self):
        """Returns the calls the data structure to return the set of susceptibles
        """
        return [nodes for nodes in self.epi_data if self.epi_data[nodes]["Infection Stage"] == "Susceptible"]

    @property
    def infectious_periods(self):
        #There's a function that generates the infection periods as it is shared between several class objects
        return [self.epi_data[node]["Infection Period"] for node in self.epi_data]

    @property
    def exposure_level(self):
        """Call to the data structure to return the current exposure level for every node in the network.
        """
        return [self.epi_data[node]["Exposure Level"] for node in self.epi_data]
    
    def calculate_total_emitted_hazard(self):
        """For every node in the network, we calculate the total hazard emitted if they were infected.
        """

        haz = hazard_class(self.hazard_rate)

        self.lifetime_emitted_hazard = [self.beta * haz.integrate_hazard(length) for length in self.infectious_periods]

    def updates_exposure_levels(self):
        """Updates the exposure levels of nodes that are connected to newly infected nodes.
        """

        for node in self.new_infections:
            connected_nodes = self.G.neighbors(node)
            node_index = self.node_keys.index(node)
            emitted_hazard = self.lifetime_emitted_hazard[node_index]

            for connected_node in connected_nodes:
                self.data_structure.update_exposure_level(connected_node, emitted_hazard)
        
    def determine_new_infections(self):
        """The infected list gets updated based upon their exposure level.
        """
        #Loops over all nodes which isn't terribly efficient.
        self.updates_exposure_levels()
        self.new_infections = [susceptible for susceptible in self.susceptible_nodes if (self.epi_data[susceptible]["Resistance"] < self.epi_data[susceptible]["Exposure Level"]) ]
        self.update_infection_stage(self.new_infections, "Infected", self.time)
    

    def iterate_epidemic(self):
        """Control structure for looping the epidemic until it completes. Only useful for estimating the final size of an SIR epidemic.
        """
        #The set of infected at the previous step of the iteration
        previously_infected = self.infected_nodes[:]
        #The nodes whose neighbours exposure levels will be updated.
        self.new_infections = self.infected_nodes[:]
        epidemic_ended = False
        while epidemic_ended == False:
            self.time += 1
            self.determine_new_infections()
            if previously_infected == self.infected_nodes:
                epidemic_ended = True
            previously_infected = self.infected_nodes[:]
        self.final_size = len(self.infected_nodes)
        
