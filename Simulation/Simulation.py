import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
import networkx as nx
from Simulation.EpidemicSimulation import epidemic_data


class hazard_class:
    """For a specified hazard function, this class manages calculations of useful quantities"""

    def __init__(self, hazard_function):
        """Initialises the class
        
        Arguments:
            hazard_function {function} -- A function of the form f(t)
        """
        self.hazard_function = hazard_function

    def hazard(self, t, t_end):
        """Returns a variant of the hazard rate function which rounds negative values up to 0.

        If t > t_end then it also returns 0.
        
        Arguments:
            t {int, float} -- evalutate the hazard rate at this t
            t_end {int, float} -- Cut-off value to return 0, typically set to the when the infection ends
        
        Returns:
            f(t)
        """

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

    def increment_hazard(self, t_0, t_1, end_of_infection_time):
        """Integrates the hazard function over the domain [t_0, t_1]
        
        Arguments:
            t_0 {float} -- The first timepoint
            t_1 {float} -- The second timepoint
            end_of_infection_time {float} -- The time at which a nodes infection will end. This is required so that values after this time are returned as 0
        """
        def f(t): return self.hazard(t, end_of_infection_time)
        hazard_emitted = spi.quad(f, t_0, t_1)
        return hazard_emitted[0]



class complex_epidemic_simulation(epidemic_data):
    """This class manages the simulation of the epidemic and dynamic network behavior."""

    def __init__(self, G, beta, infection_period_parameters, initial_infected, time_increment, max_iterations, hazard_rate=None,
                 infection_period_distribution=None, SIS = False, increment_network = None, custom_behaviour = None):
        """This class manages the simulation of the epidemic and the simulation of the dynamic network (if the network is dynamic).
        If the network is static, then
        
        Arguments:
            G {Networkx.graph} -- A NetworkX graph object
            beta {float} -- The thinning parameter 
            infection_period_parameters {list} -- A list of parameters for the infection period distribution
            initial_infected {int, list} -- [description]
            time_increment {float} -- The length of the time step of the simulation
            max_iterations {int} -- The maximum number of iterations that will be performed
        
        Keyword Arguments:
            hazard_rate {function} -- A function of the form f(x) (default: f(x) = 1)
            infection_period_distribution {function} -- A numpy random number distributon (default: {None})
            SIS {bool} -- Boolean on whether the epidemic is SIS, if not it will be treated as SIR (default: False)
            increment_network {method} -- A method of the form increment_network(increment_length). This method will be called during the simulation to move the network forward by the network_increment.
            custom_behaviour {function} -- Allows users to execute custom behaviour during the simulation. This is useful for customising the simulation to your own purposes, such as treatment scenarios. (default: {None})

        TODO: Remove the beta parameter, too confusing
        """

        self.G = G
        self.beta = beta
        self.infection_period_parameters = infection_period_parameters
        self.inf_starting = initial_infected
        self.time_increment = time_increment
        self.max_iterations = max_iterations
        self.hazard_rate = hazard_rate
        self.SIS = SIS
        self.increment_network = increment_network
        self.custom_behaviour = custom_behaviour

        self.N = nx.number_of_nodes(self.G)
        self.data_structure = epidemic_data(
            G, initial_infected, 100, infection_period_distribution, infection_period_parameters)
        self.epi_data = self.data_structure.epi_data
        self.hazard = hazard_class(self.hazard_rate)

    @property
    def infected_nodes(self):
        """Returns a list of dictionary keys for the nodes who are currently infected.
        
        Returns:
            [list] -- List of infected nodes
        """
        return [nodes for nodes in self.data_structure.epi_data if self.data_structure.epi_data[nodes]["Infection Stage"] == "Infected"]

    @property
    def susceptible_nodes(self):
        """Returns a list of dictionary keys for the nodes who are currently susceptible.
        
        Returns:
            [list] -- List of susceptible nodes
        """
        return [nodes for nodes in self.epi_data if self.epi_data[nodes]["Infection Stage"] == "Susceptible"]

    @property
    def infectious_periods(self):
        """Returns a list of the infectious periods for all nodes.
        
        This will include susceptible nodes, for whom this will be the length of the infectious period the next time they are infected.
        
        Returns:
            [list] -- A list of infectious periods
        """
        # There's a function that generates the infection periods as it is shared between several class objects
        return [self.epi_data[node]["Infection Period"] for node in self.epi_data]

    @property
    def recovered_nodes(self):
        """Returns a list of dictionary keys for the nodes who are currently recovered.
        
        Returns:
            [list] -- List of recovered nodes
        """
        return [nodes for nodes in self.epi_data if self.epi_data[nodes]["Infection Stage"] == "Recovered"]

    @property
    def exposure_level(self):
        """Returns a list of the current exposure level of every node in the network.
        
        Returns:
            [list] -- A list of node exposure levels
        """
        return [self.epi_data[node]["Exposure Level"] for node in self.epi_data]

    def updates_exposure_levels(self):
        """Loops over all infected nodes and updates the exposure levels of connected susceptible nodes.
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
        """Compares a nodes exposure level to it's resistance and determines which nodes have been infecetd during this step of the iteration.
        """
        self.new_infections = [susceptible for susceptible in self.susceptible_nodes if (
            self.epi_data[susceptible]["Resistance"] < self.epi_data[susceptible]["Exposure Level"])]
        self.update_infection_stage(self.new_infections, "Infected", self.time)

    def determine_recoveries(self):
        """For nodes whose infections have ended, this method updates to the appropiate status.
        
        Raises:
            ValueError: Raises an error if the SIS is not a boolean
        """
        recoveries = [infected for infected in self.infected_nodes
                    if self.epi_data[infected]["Infection Stage Started"] + self.epi_data[infected]["Infection Period"] < self.time]

        if self.SIS == False: 
            self.update_infection_stage(recoveries, "Recovered", self.time)
        elif self.SIS == True:
            self.update_infection_stage(recoveries, "Susceptible", self.time)
        else:
            raise ValueError("SIS parameter not set to true or false.")

    def perform_iteration(self):
        """Executes one step of the simulation in the following order:
        1) Update the network structure
        2) Determine which infections have ended
        3) Update node exposure levels
        4) Determine new infection
        5) Perform custom behaviours
        """
        
        #Computation Steps
        if self.increment_network != None:
            self.increment_network(self.time_increment)
        self.determine_recoveries()
        self.updates_exposure_levels()
        self.determine_new_infections()

        self.iteration += 1
        self.time += self.time_increment

        if self.custom_behaviour != None:
            self.custom_behaviour(self)


        #Recording data from here onwards
        self.data_time.append(self.time)

        self.data_susceptible_counts.append(len(self.susceptible_nodes))
        self.data_infected_counts.append(len(self.infected_nodes))
        self.data_recovered_counts.append(len(self.recovered_nodes))

        self.data_susceptible_nodes.append(self.susceptible_nodes)
        self.data_infected_nodes.append(self.infected_nodes)
        self.data_recovered_nodes.append(self.recovered_nodes)

        if self.infected_nodes == []:
            self.epidemic_ended = True

        if self.iteration == self.max_iterations:
            self.max_iterations_reached = True


    def iterate_epidemic(self):
        """Performs iterations of the simulation until either there is epidemic die out, or the maximum number of iterations is reached.
        """
        # The set of infected at the previous step of the iteration
        # The nodes whose neighbours exposure levels will be updated.
        #We will be recording data into these lists.

        #variables for controlling the iteration
        self.time = 0
        self.iteration = 0
        self.epidemic_ended = False
        self.max_iterations_reached = False

                #Data for the results
        self.data_time = [0]
        self.data_susceptible_counts = [len(self.susceptible_nodes)]
        self.data_infected_counts = [len(self.infected_nodes)]
        self.data_recovered_counts = [len(self.recovered_nodes)]

        self.data_susceptible_nodes = [self.susceptible_nodes]
        self.data_infected_nodes = [self.infected_nodes]
        self.data_recovered_nodes = [self.recovered_nodes]

        while (self.epidemic_ended == False) and (self.max_iterations_reached == False):
            self.perform_iteration()

        self.final_size = len(self.recovered_nodes)

        if self.epidemic_ended == True:
            self.stop_reason = f"The epidemic died out at time = {self.time} ({self.iteration} iterations)"
        else:
            self.stop_reason = f"The simulation stopped because the max number of iteration was reached (max = {self.iteration} iterations)."
