#Epidemic Simulation

"""We store the objects used specifically in performing the Sellke Simulation of epidemics on networks
"""
import numpy as np
import copy as c

class infection_period_handler:
    """[summary]
    
    Raises:
        ValueError: [description]
    """
    def __init__(self, N, infection_period_distribution = None, infection_period_parameters = None):
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
    def __init__(self, G, initial_infected, infection_period_distribution = None, infection_period_parameters = None, debug = False):
        self.G = G
        #create a list of the keys that are used in the data structure.
        self.node_keys = list(G.nodes())
        self.N = len(self.node_keys)
        self.initial_infected = initial_infected
        self.infection_period_distribution = infection_period_distribution
        self.infection_period_parameters = infection_period_parameters
        self.infection_period_handler = infection_period_handler(self.N, self.infection_period_distribution, self.infection_period_parameters)
        self.initialise_data_structure()
        self.initialise_resistances()
        self.initialise_infection()
        self.initialise_infection_periods()
        if debug == True: print("epidemic_data.__init__ has been called")


    def initialise_data_structure(self):
        
        #Create a dictionary where the keys are the node name.
        epi_data = dict(self.G.nodes())

        #We assign each node a basic dictionary of things that we want to track
        node_info = {
            #The latest set of information for a node is stored here
            "Infection Stage": None,
            "Infection Stage Started": None,
            "Resistance": 0,
            "Infection Period": 0,
            "Exposure Level": 0,
            "History": {
                #Records the times at which events occur
                "Node Created": 0,
                "Infection Stage Log": [],
                "Infection Stage Times": []
            },
        }

        for node in epi_data:
            epi_data[node] = c.deepcopy(node_info)
        
        self.epi_data = epi_data

    def initialise_infection(self):

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
    
    def initialise_infection_periods(self):
        infection_periods = self.infection_period_handler.generate()
        i=0
        for node in self.node_keys:
            self.epi_data[node].update({"Infection Period": infection_periods[i]})
            i += 1
        
    def initialise_resistances(self):
        resistances = np.random.exponential(size = self.N)
        i=0
        for node in self.node_keys:
            self.epi_data[node].update({"Resistance": resistances[i]})
            i += 1

        
    def update_infection_stage(self,node_list, new_stage, timepoint):
        """Changes the stage of a node and records the time at which this occurs.
        
        Arguments:
            node {list} -- list of nodes to be updated. Each entry must match an element in G.nodes() where G is a networkx graph.
            new_stage {str} -- The new stage to be recorded
            timepoint {float, int} -- The time at which the change occurs
        """
        
        for node in node_list:

            #update the infection stage
            self.epi_data[node].update({"Infection Stage": str(new_stage)})
            #append the new stage to the log
            new_log = self.epi_data[node]["History"]["Infection Stage Log"]
            new_log.append(new_stage)
            #update the log
            self.epi_data[node]["History"].update({"Infection Stage Log": new_log})


            self.epi_data[node].update({"Infection Stage Started": timepoint})
            new_times = self.epi_data[node]["History"]["Infection Stage Times"]
            new_times.append(timepoint)
            self.epi_data[node]["History"].update({"Infection Stage Times": new_times})

    def update_exposure_level(self, node, exposure_increment):
        """Increases a nodes exposure level by the exposure_increment
        
        Arguments:
            node {int, tuple} -- The node whose exposure level will be incremented.
            exposure_increment {int, float} -- The amount that the nodes exposure level will be increased by.
        """

        new_exposure = self.epi_data[node]["Exposure Level"] + exposure_increment
        self.epi_data[node].update({"Exposure Level": new_exposure})
        