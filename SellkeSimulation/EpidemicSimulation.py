#Epidemic Simulation

"""We store the objects used specifically in performing the Sellke Simulation of epidemics on networks
"""
import numpy as np
import copy as c

class epidemic_data:
    def __init__(self, G, initial_infected,debug = False):
        self.G = G
        #create a list of the keys that are used in the data structure.
        self.node_keys = list(G.nodes())
        self.initial_infected = initial_infected
        self.initialise_data_structure()
        self.initialise_infection()
        if debug == True: print("epidemic_data.__init__ has been called")


    def initialise_data_structure(self):
        
        #Create a dictionary where the keys are the node name.
        epi_data = dict(self.G.nodes())

        #We assign each node a basic dictionary
        node_info = {
            #The latest set of information for a node is stored here
            "Infection Stage": None,
            "Infection Stage Started": None,
            "Resistance": None,
            "Infection Length": 0,
            
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
    
    def update_infection_stage(self,node_list, new_stage, timepoint):
        """Changes the stage of a node and records the time at which this occurs.
        
        Arguments:
            node {list} -- list of nodes to be updated. Each entry must match an element in G.nodes() where G is a networkx graph.
            new_stage {str} -- The new stage to be recorded
            timepoint {float, int} -- The time at which the change occurs
        """
        
        for node in node_list:

            #update the infection stage
            self.epi_data[node].update({"Infection Stage": new_stage})
            #append the new stage to the log
            new_log = self.epi_data[node]["History"]["Infection Stage Log"]
            new_log.append(new_stage)
            #update the log
            self.epi_data[node]["History"].update({"Infection Stage Log": new_log})


            self.epi_data[node].update({"Infection Stage Started": timepoint})
            new_times = self.epi_data[node]["History"]["Infection Stage Times"]
            new_times.append(timepoint)
            self.epi_data[node]["History"].update({"Infection Stage Times": new_times})

    def initialise_infection(self):

        if type(self.initial_infected) == int:

            #Randomly choose the initial infected
            key_index = len(self.node_keys)
            self.initial_infected = np.random.choice(key_index, replace = False, size = self.initial_infected)

            #Set the initial infected stage'
            self.initial_infected = self.node_keys[self.initial_infected]
            self.update_infection_stage(self.initial_infected, "Infected", 0)

        self.update_infection_stage(self.initial_infected, "Infected", 0)
        initial_susceptibles = [node for node in self.epi_data if self.epi_data[node]["Infection Stage"] != "Infected"]
        self.update_infection_stage(initial_susceptibles, "Susceptible", 0)
