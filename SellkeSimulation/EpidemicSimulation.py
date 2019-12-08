import numpy as np
import copy as c

class infection_period_handler:

    def __init__(self, N, infection_period_distribution = None, infection_period_parameters = None):
        """The class calls numpy.random distributions differently, depending on how many parameters are being passed to it.

        Examples:

        If parameters = [parameter_1] then the class will attempt to return:
        infection_period_distribution(parameter_1, N)

        If parameters = [parameter_1, parameter_2] then the class will attempt to return:
        infection_period_distribution(parameter_1, parameter_2, N)
        
        Arguments:
            N {int} -- The number of observations to draw from the random distribution
        
        Keyword Arguments:
            infection_period_distribution {function} -- A numpy.random distribution (default: {exponential})
            infection_period_parameters {list} -- A list of parameters (default: {None})
        """
        self.infection_period_distribution = infection_period_distribution
        self.infection_period_parameters = infection_period_parameters
        self.N = N
        if infection_period_parameters == None: self.infection_period_parameters = 1 
    
    def generate(self):
        if self.infection_period_distribution is None:
            if type(self.infection_period_parameters) == int or type(self.infection_period_parameters) == float:
                self.infection_periods = np.random.exponential(self.infection_period_parameters,self.N)
            else:
                print("Put something here to stop everything else going ahead because the infection_period_parameters are not correct.")
        else:
            if type(self.infection_period_parameters) == int or type(self.infection_period_parameters) == float:
                self.infection_periods = self.infection_period_distribution(self.infection_period_parameters, self.N)
            elif len(self.infection_period_parameters) == 2:
                self.infection_periods = self.infection_period_distribution(self.infection_period_parameters[0]
                                                                    ,self.infection_period_parameters[1]
                                                                    ,self.N)
            elif len(self.infection_period_parameters) == 3:
                self.infection_periods = self.infection_period_distribution(self.infection_period_parameters[0]
                                                                    ,self.infection_period_parameters[1]
                                                                    ,self.infection_period_parameters[2]
                                                                    ,self.N)
            else:
                print("There is something incorrect with the infection_period_parameters.")
        #if any(infection_periods < 0):
        #    raise ValueError("Negative values for length of infectious period detected. Ensure use of a positive distribution.")
        return(self.infection_periods)



class epidemic_data(infection_period_handler):
    def __init__(self, G, initial_infected, pre_gen_data, infection_period_distribution = None, infection_period_parameters = None, treatment_class = False, treatment_dist = None):
        """A class used to store the data about the epidemic. Includes a number of methods to easily update the data, and return useful data sets.

        Note:
        We generate all the variables about the epidemic first. This is so that you can run exactly the same epidemic, without other random variables incrementing the random number generator.
        For example, suppose you want to test a treatment, and this treatment requires generating random numbers this would increment the random number generator seed.
        If we were generating random number on the fly for the epidemic, you would end up with different epidemics and therefore a much larger number of simulations required to generate good statistical values.

        In effect, by generating the epidemic random numbers first, you are able to "rewind time" and test your treatment on exactly the same epidemic.
        
        Arguments:
            infection_period_handler {[type]} -- [description]
            G {NetworkX graph} -- The network that will be used to initialise the node data
            initial_infected {int, list} -- Either a number of nodes to be randomly infected at time 0, or a list of nodes who will be infected at time 0
            pre_gen_data {int} -- The number of times we draw a variable for the data generation. If pre-gen-data = 10, then a node will have enough pre-generated data to be infected and recover 10 times.
        
        Keyword Arguments:
            infection_period_distribution {function} -- The distribution that will be used to generate the length of an infection period (default: exponential)
            infection_period_parameters {list} -- A list of parameters to be passed to the infection period distribution (default: 1)
        """
        self.G = G
        self.node_keys = list(G.nodes())
        self.N = len(self.node_keys)
        self.pre_gen_data = pre_gen_data
        self.initial_infected = initial_infected
        self.infection_period_distribution = infection_period_distribution
        self.infection_period_parameters = infection_period_parameters
        self.infection_period_handler = infection_period_handler(self.pre_gen_data, self.infection_period_distribution, self.infection_period_parameters)
        self.initialise_data_structure()
        self.pre_generate_data()
        self.initialise_infection()


    def initialise_data_structure(self):
        """Creates an initial dictionary of parameters for each node with default values.
        """
        
        #Create a dictionary where the keys are the node name.
        epi_data = dict(self.G.nodes())

        #This is a basic set of information we need to know for each node.
        node_info = {
            #The current status of the node is stored at the top level.
            "Infection Stage": None,
            "Infection Stage Started": None,
            "Resistance": 0,
            "Infection Period": 0,
            "Exposure Level": 0,
            "Times Infected": 0,
            "Times Susceptible": 0,

            #The History sub dictionary is updated whenever the status of a node changes
            "History": {
                #Records the times at which events occur
                "Node Created": 0,
                "Infection Stage Log": [],
                "Infection Stage Times": []
            },

            #Pre-generated Data is stored here, whenever the nodes status is updated, the new value is taken form the sequence defined in this sub-dictionary.
            "Pre-generated Data": {
                "Resistance": [],
                "Infection Period": []
            }
        }

        for node in epi_data:
            #I think the deepcopy can now just be changed to a copy.
            epi_data[node] = c.deepcopy(node_info)
        
        self.epi_data = epi_data

    def initialise_infection(self):
        """This method initialises the infection, by updating the node status to 0.
        If initial infected is an integer, then the nodes chosen to be infected at time 0 are chosen at random.
        If a list of NetworkX dictionary keys for the nodes is supplied, then the specified nodes are chosen to be infected.
        """

        if type(self.initial_infected) == int:

            #Randomly choose the initial infected
            key_index = len(self.node_keys)
            self.initial_infected = np.random.choice(key_index, replace = False, size = self.initial_infected)

            #Set the initial infected stage'
            self.initial_infected = [self.node_keys[index] for index in self.initial_infected]
            self.update_infection_stage(self.initial_infected, "Infected", 0)

        self.update_infection_stage(self.initial_infected, "Infected", 0)
        initial_susceptibles = [node for node in self.epi_data if self.epi_data[node]["Infection Stage"] != "Infected"]
        self.update_infection_stage(initial_susceptibles, "Susceptible", 0)
        
    def update_infection_stage(self,node_list, new_stage, timepoint):
        """Allows you to update the status of a node, and records the times at which this occurs.
        
        Arguments:
            node {list} -- list of node dictionary keys, the specified nodes will be updated
            new_stage {str} -- The new infection stage the nodes will be updated to
            timepoint {float, int} -- The time at which the change occurs

        TODO: Input is not list should work
        """

        for node in node_list:
            
            #update the infection stage
            self.epi_data[node].update({"Infection Stage": str(new_stage)})

            #If the new stage is susceptible, we give them a new resistance value and set their exposure to 0.
            if new_stage == "Susceptible":
                #How many times have they been in the susceptible state
                times_susceptible = self.epi_data[node]["Times Susceptible"]
                #Get the resistance for that susceptible state
                new_resistance = self.epi_data[node]["Pre-generated Data"]["Resistance"][times_susceptible]
                #Update their current resistance to the new resistance
                self.epi_data[node].update({"Resistance": new_resistance})
                #Add one to the number of times they've been in the susceptible state
                self.epi_data[node].update({"Times Susceptible": times_susceptible + 1})
                #Set their exposure level to 0
                self.epi_data[node].update({"Exposure Level": 0})
            
            #If the new stage is infected, we give them a new infection period value and set their exposure to 0.
            if new_stage == "Infected":
                #How many times have they been in the susceptible state
                times_infected = self.epi_data[node]["Times Infected"]
                #Get the resistance for that susceptible state
                new_infection_period = self.epi_data[node]["Pre-generated Data"]["Infection Period"][times_infected]
                #Update their current resistance to the new resistance
                self.epi_data[node].update({"Infection Period": new_infection_period})
                #Add one to the number of times they've been in the susceptible state
                self.epi_data[node].update({"Times Infected": times_infected + 1})

            #Update the infection stage history
            new_log = self.epi_data[node]["History"]["Infection Stage Log"]
            new_log.append(new_stage)
            self.epi_data[node]["History"].update({"Infection Stage Log": new_log})

            #Update the the timepoints.
            self.epi_data[node].update({"Infection Stage Started": timepoint})
            new_times = self.epi_data[node]["History"]["Infection Stage Times"]
            new_times.append(timepoint)
            self.epi_data[node]["History"].update({"Infection Stage Times": new_times})

    def update_exposure_level(self, node, exposure_increment):
        """Increases a nodes exposure level by an amount equal to the exposure_increment
        
        Arguments:
            node {int, tuple} -- The dictionary key of the node who is receiving the exposure
            exposure_increment {int, float} -- The amount that the nodes exposure level will be increased by
        """

        new_exposure = self.epi_data[node]["Exposure Level"] + exposure_increment
        self.epi_data[node].update({"Exposure Level": new_exposure})
        
    def pre_generate_data(self):
        """This method pre-generates the infection periods and resistances of a node. This is important, because we want to be able to re-run epidemics with an intervention to see how effective it is.
        
        Arguments:
            count {int} -- How many infection periods to generate. If the number of infections for a node exceeds the pre generated data, an error will be thrown.

        TODO Throw error if infection periods run out
        """
        for node in self.node_keys:
            resistances = list(np.random.exponential(1,self.pre_gen_data))
            self.epi_data[node]["Pre-generated Data"].update({"Resistance": resistances})

            inf_periods = list(self.infection_period_handler.generate())
            self.epi_data[node]["Pre-generated Data"].update({"Infection Period": inf_periods})
