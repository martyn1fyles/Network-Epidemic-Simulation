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
        Returns 0 if t is less than 0, truncates to 0 after t_end
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
        The function calculates the rate at which the hazard is emitted at time t, assuming they were all intially infected

        Designed as an input function for the BVP Solver, which doesn't currently work
        '''
        self.T_endpoints = T_endpoints
        total_hazard_rate = 0
        for i in range(len(T_endpoints)):
            total_hazard_rate = total_hazard_rate + \
                self.hazard(t, T_endpoints[i])
        return total_hazard_rate

    def integrate_hazard(self, t_end):
        '''
        Integrate a hazard function
        '''
        def f(t): return self.hazard(t, t_end)
        integral = spi.quad(f, 0, t_end)
        return integral[0]

    def increment_hazard(self, t_0, t_1, end_of_infection_time):
        """Integrates the hazard between the times of t_0 and t_1, with t_0 < t_1

        Arguments:
            t_0 {float} -- The first timepoint
            t_1 {float} -- The seconde timepoint
        """
        def f(t): return self.hazard(t, end_of_infection_time)
        hazard_emitted = spi.quad(f, t_0, t_1)
        return hazard_emitted[0]


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

    def __init__(self, G, beta, I_parameters, initial_infected, time_increment, max_iterations, hazard_rate=None, infection_period_distribution=None, SIS = None):
        self.G = G
        self.node_keys = list(G.nodes())
        self.beta = beta
        self.time_increment = time_increment
        self.inf_starting = initial_infected
        self.I_parameters = I_parameters
        self.hazard_rate = hazard_rate
        self.N = nx.number_of_nodes(self.G)
        self.data_structure = epidemic_data(
            G, initial_infected, 100, infection_period_distribution, I_parameters)
        self.SIS = SIS
        if self.SIS == None:
            self.SIS = False
        self.epi_data = self.data_structure.epi_data
        self.new_infections = self.infected_nodes[:]
        self.calculate_total_emitted_hazard()
        self.time = 0
        self.hazard = hazard_class(self.hazard_rate)
        self.max_iterations = max_iterations

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
        # There's a function that generates the infection periods as it is shared between several class objects
        return [self.epi_data[node]["Infection Period"] for node in self.epi_data]

    @property
    def recovered_nodes(self):
        """Returns the calls the data structure to return the set of susceptibles
        """
        return [nodes for nodes in self.epi_data if self.epi_data[nodes]["Infection Stage"] == "Recovered"]

    @property
    def exposure_level(self):
        """Call to the data structure to return the current exposure level for every node in the network.
        """
        return [self.epi_data[node]["Exposure Level"] for node in self.epi_data]

    def calculate_total_emitted_hazard(self):
        """For every node in the network, we calculate the total hazard emitted if they were infected.

        DEPRECATED
        """

        haz = hazard_class(self.hazard_rate)

        self.lifetime_emitted_hazard = [
            self.beta * haz.integrate_hazard(length) for length in self.infectious_periods]

    def updates_exposure_levels(self):
        """Updates the exposure levels of nodes that are connected to newly infected nodes. This should only affect susceptible nodes.
        """

        for node in self.infected_nodes:

            infection_started = self.epi_data[node]["Infection Stage Started"]
            end_of_infection = infection_started +\
                self.epi_data[node]["Infection Period"]
            time_since_infected = self.time - infection_started
            emitted_hazard = self.beta * self.hazard.increment_hazard(
                time_since_infected, time_since_infected + self.time_increment, end_of_infection)

            connected_susceptibles = list(
                set(self.G.neighbors(node)) & set(self.susceptible_nodes))

            for exposed_node in connected_susceptibles:
                self.data_structure.update_exposure_level(
                    exposed_node, emitted_hazard)

    def determine_new_infections(self):
        """The infected list gets updated based upon their exposure level.
        """
        self.new_infections = [susceptible for susceptible in self.susceptible_nodes if (
            self.epi_data[susceptible]["Resistance"] < self.epi_data[susceptible]["Exposure Level"])]
        self.update_infection_stage(self.new_infections, "Infected", self.time)

    def determine_recoveries(self):
        """Nodes if the amount of time since the infected started is greater than the length of the infection period, the nodes status changes to Recovered.
        """
        recoveries = [infected for infected in self.infected_nodes
                    if self.epi_data[infected]["Infection Stage Started"] + self.epi_data[infected]["Infection Period"] < self.time]

        if self.SIS == False: 
            self.update_infection_stage(recoveries, "Recovered", self.time)
        elif self.SIS == True:
            self.update_infection_stage(recoveries, "Susceptible", self.time)
        else:
            raise ValueError("SIS parameter not set to true or false.")


    def iterate_epidemic(self):
        """Control structure for looping the epidemic until it completes. Only useful for estimating the final size of an SIR epidemic.
        """
        # The set of infected at the previous step of the iteration
        # The nodes whose neighbours exposure levels will be updated.

        iteration = 0
        epidemic_ended = False
        max_iterations_reached = False

        #We will be recording data into these lists.
        time = [self.time]
        susceptible_count = [len(self.susceptible_nodes)]
        infected_count = [len(self.infected_nodes)]
        recovered_count = [len(self.recovered_nodes)]

        infected_set = [self.infected_nodes]
        suscuptible_set = [self.susceptible_nodes]
        recovered_set = [self.recovered_nodes]

        while (epidemic_ended == False) and (max_iterations_reached == False):

            self.determine_recoveries()
            self.updates_exposure_levels()
            self.determine_new_infections()

            iteration += 1
            self.time += self.time_increment

            time.append(self.time)
            infected_count.append(len(self.infected_nodes))
            susceptible_count.append(len(self.susceptible_nodes))
            recovered_count.append(len(self.recovered_nodes))

            suscuptible_set.append(self.susceptible_nodes)
            infected_set.append(self.infected_nodes)
            recovered_set.append(self.recovered_nodes)

            if self.infected_nodes == []:
                epidemic_ended = True

            if iteration == self.max_iterations:
                max_iterations_reached = True

        self.final_size = len(self.recovered_nodes)
        self.iterations = iteration
        self.data_time = time

        self.data_susceptible_counts = susceptible_count
        self.data_infected_counts = infected_count
        self.data_recovered_counts = recovered_count

        self.data_susceptible_nodes = suscuptible_set
        self.data_infected_nodes = infected_set
        self.data_recovered_nodes = recovered_set

        if epidemic_ended == True:
            self.stop_reason = f"The epidemic died out at time = {self.time} ({iteration} iterations)"
        else:
            self.stop_reason = f"The simulation stopped because the max number of iteration was reached (max = {iteration} iterations)."
