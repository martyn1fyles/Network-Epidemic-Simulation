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
        """ For a given node, generate the times at which they will migrate between the blocks.
        
        Arguments:
            node {str, int, tuple} -- The dictionary key of the node whose times will be generated
        
        Keyword Arguments:
            birth_time {int, float} -- The time at which the node was added to the network (default: {0})
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
        """Adds variables to a nodes dictionary that are required for quick calculation in the future.

        These variables are:
        Membership Data - The sequence of blocks that the node will belong to
        Current Membership Index - The index of "Membership Data" which gives the current location in the sequence of block memberships
        Next Migration Time - The time at which the node will change membership
        
        Arguments:
            node {str, int, tuple} -- The dictionary key of the node whose times will be generated
        
        Keyword Arguments:
            birth_time {int, float} -- The time at which the node was added to the network (default: {0})
        """

        membership_data = self.generate_migration_times(node)
        self.G.nodes[node].update({"Membership Data": membership_data})
        self.G.nodes[node].update({"Current Membership Index": 0})

        next_time = self.get_node_migration_times(node)[1]
        self.G.nodes[node].update({"Next Migration Time": next_time})

        if self.custom_attribute != None:
            self.G.nodes[node].update(self.custom_attribute)

    def get_node_current_block(self, node):
        """Returns the current block membership of a node
        
        Arguments:
            node {str, int, float} -- The dictionary key of the node whose current block membership will be returned
        
        Returns:
            {str, int, float} -- The dictionary key representing of the block which the node currently belongs to
        """
        return self.G.nodes[node]["block"]
    
    def get_node_migration_times(self, node):
        """Returns the time at which a node will migrate between block in the network
        
        Arguments:
            node {str, int, tuple} -- The dictionary key of a node in the network
        
        Returns:
            list -- The list of times at which a node will move between blocks.
        """
        return [data[1] for data in self.G.nodes[node]["Membership Data"]]

    def get_node_memberships(self, node):
        """Returns the sequence of blocks that the node will belong to
        
        Arguments:
            node {str, int, tuple} -- The dictionary key of a node in the network
        
        Returns:
            list -- The sequence of blocks that the node will belong to
        """
        return [data[0] for data in self.G.nodes[node]["Membership Data"]]

    def get_node_next_migration_time(self, node):
        """Returns the next time that a node will migrate to another block
        
        Arguments:
            node {str, int, tuple} -- The dictionary key of a node in the network
        
        Returns:
            float -- The time next time at which a migration will occur for the chosen node
        """
        return self.G.nodes[node]["Next Migration Time"]

    def get_next_migration_times(self):
        """Returns a list of the times of next migration for every node in the network
        
        Returns:
            list -- The list of times at which the next migration will occur for every node in the network
        """
        return [self.G.nodes[node]["Next Migration Time"] for node in self.G.nodes()]
        
    def determine_nodes_to_migrate(self, time_limit):
        """For a specified time_limit, this returns a list of nodes whose next migration event occurs before the specified time_limit
        
        Arguments:
            time_limit {int, float} -- If a nodes next migration < time_limit, then the nodes dictionary key will be added to the returned list
        
        Returns:
            list -- A list of nodes keys
        """
        return [node for node in self.G.nodes if self.get_node_next_migration_time(node) < time_limit]
        
    def perform_migration_event(self, node):
        """For a specified node, perform their next migration by updating the nodes dictionary entry, removing all edges and adding new edges based upon the specified probabilities
        
        Arguments:
            node {str, int, tuple} -- The dictionary key of the node
        """
        #Get the current membership index, add one to it, and update the dictionary parameter
        index = self.G.nodes[node]["Current Membership Index"]
        new_index = index + 1
        self.G.nodes[node].update({"Current Membership Index": new_index})

        # Using the new index, update the next migration time variable
        migration_times = self.get_node_migration_times(node)
        self.G.nodes[node].update({"Next Migration Time": migration_times[new_index + 1]})

        # Using the new index, update the blocks membership
        memberships = self.get_node_memberships(node)
        self.G.nodes[node].update({"block": memberships[new_index]})

        # The node now has it's new block membership, and the edges will be added based upon the network parameters
        self.update_edges(node)

        # If the user has specified a command to be executed upon migrations, then it will be executed
        if self.custom_migration_behaviour != None:
            self.custom_migration_behaviour(self, node)
    
    def update_edges(self, node):
        """Removes all edges from a nodes and added new edges based upon the networks edge probabilities
        
        Arguments:
            node {str, int, float} -- The dictionary key of the node
        """
        #Remove all the edges from the node
        neighbours = list(self.G.neighbors(node))
        [self.G.remove_edge(node, connected_node) for connected_node in neighbours]

        # Get the block membership of the node
        node_membership = self.get_node_current_block(node)

        # Loop over every other node in the network
        for potential_neighbour in self.G.nodes:
            if potential_neighbour != node:
                
                # Get the block membership of the other nodes in the network
                potential_neighbour_membership = self.get_node_current_block(potential_neighbour)

                #Based upon the returned block memberships, we perform a bernoulli trial and add an edge
                edge_forming_prob = self.p[node_membership][potential_neighbour_membership]

                if np.random.binomial(1, edge_forming_prob) == 1:
                    self.G.add_edge(node, potential_neighbour)

    def increment_network(self, increment_length):
        """Increment the network foraward in time
        
        Arguments:
            increment_length {int, float} -- The length of time to move the network forward
        """
        self.time += increment_length

        #Which nodes need to have a migration?
        to_be_migrated = self.determine_nodes_to_migrate(self.time)
        #This while loop is important, incase multiple migration events happen during a time increment.
        while to_be_migrated != []:
            [self.perform_migration_event(node) for node in to_be_migrated]
            to_be_migrated = self.determine_nodes_to_migrate(self.time)
        

