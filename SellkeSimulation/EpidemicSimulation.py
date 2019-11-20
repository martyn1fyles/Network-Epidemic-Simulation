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
    def __init__(self, G, initial_infected, pre_gen_data, infection_period_distribution = None, infection_period_parameters = None):
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
            "Times Infected": 0,
            "Times Susceptible": 0,
            "History": {
                #Records the times at which events occur
                "Node Created": 0,
                "Infection Stage Log": [],
                "Infection Stage Times": []
            },
            "Pre-generated Data": {
                "Resistance": [],
                "Infection Period": []
            }
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
        """Increases a nodes exposure level by the exposure_increment
        
        Arguments:
            node {int, tuple} -- The node whose exposure level will be incremented.
            exposure_increment {int, float} -- The amount that the nodes exposure level will be increased by.
        """

        new_exposure = self.epi_data[node]["Exposure Level"] + exposure_increment
        self.epi_data[node].update({"Exposure Level": new_exposure})
        
    def pre_generate_data(self):
        """This method pre-generates the infection periods. This is important, because we want to be able to re-run epidemics with an intervention to see how effective it is.
        
        Arguments:
            count {int} -- How many infection periods to generate. If the number of infections for a node exceeds the pre generated data, an error will be thrown.
        """
        for node in self.node_keys:
            resistances = list(np.random.exponential(1,self.pre_gen_data))
            self.epi_data[node]["Pre-generated Data"].update({"Resistance": resistances})

            inf_periods = list(self.infection_period_handler.generate())
            self.epi_data[node]["Pre-generated Data"].update({"Infection Period": inf_periods})
