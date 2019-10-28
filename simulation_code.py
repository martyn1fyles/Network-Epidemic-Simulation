import numpy as np
import matplotlib.pyplot as plt

class SIR_Selke:
    '''
    This method is able to simulate the final size of an SIR epidemics with non-MArkovian distribution. Further, it has tools for the investigation of the final size of the epidemic
    such as being able to repeatedly generate observations, create histograms, plot trajectories.
    
    Arguements
    N = total size of population
    beta = force of infection
    I_parameters = either a float, int, or list of parameters that are passed to the numpy distribution function
    
    Choosing a non-markovian distribution:
    If parameter non_markovian_distribution is left blank then we take it to be the exponential distribution
    In order to change the distribution to something non-markovian, pass the name of a numpy.random univariate function to the class
    '''
    
    #We use the init function to assign values to object that are necessary to do when the object is run
    def __init__(self, N, beta, I_parameters, infected_started, non_markovian_dist = None):
        self.N = N
        self.beta = beta
        self.inf_starting = infected_started
        self.I_parameters = I_parameters
        self.non_markovian_dist = non_markovian_dist
        if self.non_markovian_dist == None:
            print("Taking the exponential distribution for the infectious periods.")
        
    def generate_infection_periods(self):
        '''
        This method contains the logic required to handle change of distributions.
        Defaults to choosing the exponential distribution if no other distribution is specified.
        Different distributions have a different number of parameters passed to them 
        
        To Do:
        We need something to make sure that the distributions that are passed to this function are non-negative.
        Rename the self.datas variables to something more intelligible
        Add some proper break variables to this whole business
        Make the error messages better
        '''
        if self.non_markovian_dist is None:
            if type(self.I_parameters) == int or type(self.I_parameters) == float:
                self.datas = np.random.exponential(self.I_parameters,self.N)
            else:
                print("Put something here to stop everything else going ahead because te parametes are not correct.")
        else:
            if type(self.I_parameters) == int or type(self.I_parameters) == float:
                self.datas = self.non_markovian_dist(self.I_parameters, self.N)
            elif len(self.I_parameters) == 2:
                self.datas = self.non_markovian_dist(self.I_parameters[0]
                                                    ,self.I_parameters[1]
                                                    ,self.N)
            elif len(self.I_parameters) == 3:
                self.datas = self.non_markovian_dist(self.I_parameters[0]
                                                    ,self.I_parameters[1]
                                                    ,self.I_parameters[2]
                                                    ,self.N)
            else:
                print("nah ya fucked it matey")
            
        return(self.datas)
    
    
    def gen_final_size(self):
        '''Generates 1 observation of the final size of an epidemic from the parameters defined earlier'''
        #Infectious Period variables, function as defined in the __init__ section
        T = self.generate_infection_periods()
        assert all(T) > 0

        #The resistance to infection variables, always an exponential 1 rv
        Q = np.random.exponential(scale = 1,
                                  size = self.N)
        
        #We take the formulation found in the Thomas House paper and set the intial infecteds to zero
        Q[ :self.inf_starting] = 0
        
        #The ordering
        ordering = np.argsort(Q)  #The vector that gives us the correct ordering, incase we want to use it later
        Q_sorted = Q[ordering]    #The sorted vector
        T_cumulative = self.beta * np.cumsum(T) #The cumulative sum of the tolerance 
        
        counter = 0
        for i in range(self.N-1):
            if Q_sorted[i] < T_cumulative[i]:
                counter = counter + 1
            else:
                break
        return(counter+1)
    
    def sim_final_size(self, n_sim):
        '''Generates multiple observations of the final size of the epidemic.'''
        self.n_sim = n_sim
        
        self.observations = []
        
        print("Performing", self.n_sim, "iterations!")
        for i in range(self.n_sim):
            
            new_obs = self.gen_final_size()
            self.observations.append(new_obs)
        
        #return(self.observations)
        
    def plot_hist(self):
        '''Calls the method for generating the observations and then creates the '''
        self.plot = plt.hist(self.observations, bins = range(self.N))
        plt.hist(self.observations, bins = range(self.N), density = True)
        plt.title(f'Final epidemic size of {self.n_sim} observations')
        
