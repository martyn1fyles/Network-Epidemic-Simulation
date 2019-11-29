#This module contains code we use to construct dynamic networks
import copy as c
import networkx as nx
import numpy as np

class dynamic_stochastic_block_model:
    """This class enables dynamics for Stochastic Block Models (SBM) in the form of Birth and Death Processes and migration, where nodes are allowed to move between groups at random times.
    
    We require that an end-time is specified, as the data cannot be generated on the fly without impacting the reproducibility of an experiment, in the future, we can relax this restraint."""
    
    def __init__(self, sizes, p, m, waiting_time_par, end_time, node_list = None, birth_rate = 0, custom_attribute = None, custom_migration_behaviour = None):
        #I have dropped the directed parameter since I cannot think of a simple way to implement it.
        self.G = nx.generators.community.stochastic_block_model(sizes, p, node_list)
        self.end_time = end_time
        self.waiting_time_par = waiting_time_par
        self.m = m
        self.p = p
        self.time = 0
        self.custom_migration_behaviour = custom_migration_behaviour
        self.custom_attribute = custom_attribute
        [self.assign_membership_data(node) for node in self.G.nodes]
    
    def generate_migration_times(self, node, birth_time = 0):
        """For a single node, this method returns the times at which the node moves between groups.

        In the future, we will make it so that the length of time a node stays in a group is conditional upon the group.

        For now, we will only allow exponentially distributed waiting times.
        """

        current_block = self.G.nodes[node]["block"]
        time = birth_time
        memberships = [(current_block, 0)] #(which block are they in, what time they started being in that block)
        while time < self.end_time:
            length_of_stay = np.random.exponential(self.waiting_time_par)
            time = time + length_of_stay

            # We draw a sample from a multinomial distribution to determine the new group. The migration matrix m is used to sample the migration probabilities. We then work out the index to find out which group to move to
            draw_sample = list(np.random.multinomial(1,self.m[current_block]))
            new_block = draw_sample.index(1)

            memberships.append((new_block, time))

            current_block = new_block
        return(memberships)

    def assign_membership_data(self, node, birth_time = 0):
        #This method calls the method for data generation or migration times, and assigns the data to the node in the dictionary
        #It also sets up other variabls necessary for fast updates
        #The current membership index simply stores the index of the tuple which represents the current membership of the node
        #This allows to quickly look at the next tuple, to work out migration times.

        membership_data = self.generate_migration_times(node)
        self.G.nodes[node].update({"Membership Data": membership_data})
        self.G.nodes[node].update({"Current Membership Index": 0})

        next_time = self.get_node_migration_times(node)[1]
        self.G.nodes[node].update({"Next Migration Time": next_time})

        if self.custom_attribute != None:
            self.G.nodes[node].update(self.custom_attribute)


    def get_current_node_membership_index(self, node):
        """Returns the index for the tuple in Membership data that represents the nodes current block membership
        
        Arguments:
            node {str, int, tuple} -- the NetworkX name for a node in the network
        
        Returns:
            int -- The index for the membership tuple that represents the nodes current block membership
        """
        return self.G.nodes[node]["Current Membership Index"]

    def get_node_current_block(self, node):
        return self.G.nodes[node]["block"]
    
    def get_node_migration_times(self, node):
        """Provides a list of times at which the node will move between blocks.
        
        Arguments:
            node {str, int, tuple} -- the NetworkX name for a node in the network
        
        Returns:
            list -- The list of times at which a node will move between blocks.
        """
        return [data[1] for data in self.G.nodes[node]["Membership Data"]]

    def get_node_memberships(self, node):
        """Provides the list of block that a node will belong to
        
        Arguments:
            node {str, int, tuple} -- the NetworkX name for a node in the network
        
        Returns:
            list -- The list of times at which a node will move between blocks.
        """
        return [data[0] for data in self.G.nodes[node]["Membership Data"]]

    def get_node_next_migration_time(self, node):
        return self.G.nodes[node]["Next Migration Time"]

    def get_next_migration_times(self):
        return [self.G.nodes[node]["Next Migration Time"] for node in self.G.nodes()]
        
    def determine_nodes_to_migrate(self, time):
        return [node for node in self.G.nodes if self.get_node_next_migration_time(node) < time]
        
    def perform_migration_event(self, node):
        #Get the current index, add one to it, and update
        index = self.G.nodes[node]["Current Membership Index"]
        new_index = index + 1
        self.G.nodes[node].update({"Current Membership Index": new_index})

        # Using the new index, update the next migration time variable
        migration_times = self.get_node_migration_times(node)
        self.G.nodes[node].update({"Next Migration Time": migration_times[new_index + 1]})

        # Using the new index, update the blocks membership
        memberships = self.get_node_memberships(node)
        self.G.nodes[node].update({"block": memberships[new_index]})

        self.update_edges(node)

        if self.custom_migration_behaviour != None:
            self.custom_migration_behaviour(self, node)
    
    def update_edges(self, node):
        #Remove all the edges from the node
        neighbours = list(self.G.neighbors(node))
        [self.G.remove_edge(node, connected_node) for connected_node in neighbours]

        node_membership = self.get_node_current_block(node)
        for potential_neighbour in self.G.nodes:
            if potential_neighbour != node:

                potential_neighbour_membership = self.get_node_current_block(potential_neighbour)

                edge_forming_prob = self.p[node_membership][potential_neighbour_membership]

                if np.random.binomial(1, edge_forming_prob) == 1:
                    self.G.add_edge(node, potential_neighbour)
            
    def increment_network(self, increment_length):
        self.time += increment_length

        #Which nodes need to have a migration?
        to_be_migrated = self.determine_nodes_to_migrate(self.time)
        #This while loop is important, incase multiple migration events happen during a time increment.
        while to_be_migrated != []:
            [self.perform_migration_event(node) for node in to_be_migrated]
            to_be_migrated = self.determine_nodes_to_migrate(self.time)
        

