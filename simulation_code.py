import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi

class SIR_Selke:
    '''
    This method is able to simulate the final size of an SIR epidemics with non-MArkovian distribution. Further, it has tools for the investigation of the final size of the epidemic
    such as being able to repeatedly generate observations, create histograms, plot trajectories.
    
    Arguments
    N = total size of population
    beta = force of infection
    I_parameters = either a float, int, or list of parameters that are passed to a numpy distribution function
    

    Choosing a non-markovian distribution:
    If parameter infectious_period_distribution is left blank then we take it to be the exponential distribution
    In order to change the distribution to something non-markovian, pass the name of a numpy.random univariate function to the class
    '''
    
    #We use the init function to assign values to object that are necessary to do when the object is run
    def __init__(self, N, beta, I_parameters, infected_started, hazard_rate = None, infection_period_distribution = None):
        self.N = N
        self.beta = beta
        self.inf_starting = infected_started
        self.I_parameters = I_parameters
        self.inf_period_dist = infection_period_distribution
        self.hazard_rate = hazard_rate
        if self.inf_period_dist == None:
            print("Taking the exponential distribution of the length of the infectious periods.")
        
    def generate_infection_periods(self):
        '''
        This method contains the logic required to handle change of infectious period distributions.
        Defaults to choosing the exponential distribution if no other distribution is specified.
        Different distributions have a different number of parameters passed to them 
        
        To Do:
        Add some proper break variables to this whole business if we reach an error
        Make the error messages better
        '''
        if self.inf_period_dist is None:
            if type(self.I_parameters) == int or type(self.I_parameters) == float:
                self.inf_periods = np.random.exponential(self.I_parameters,self.N)
            else:
                print("Put something here to stop everything else going ahead because te parametes are not correct.")
        else:
            if type(self.I_parameters) == int or type(self.I_parameters) == float:
                self.inf_periods = self.inf_period_dist(self.I_parameters, self.N)
            elif len(self.I_parameters) == 2:
                self.inf_periods = self.inf_period_dist(self.I_parameters[0]
                                                    ,self.I_parameters[1]
                                                    ,self.N)
            elif len(self.I_parameters) == 3:
                self.inf_periods = self.inf_period_dist(self.I_parameters[0]
                                                    ,self.I_parameters[1]
                                                    ,self.I_parameters[2]
                                                    ,self.N)
            else:
                print("There is something incorrect with the parameters.")
            
        return(self.inf_periods)
    
    def hazard(self, t, t_end):
        '''A housekeeping variant of the hazard function.
        Retuns 0 if t is less than 0, truncates to 0 after t_end
        If no hazard rate was specified, returns 1 so the system defaults to exponential waiting times.
        '''
        
        if t < 0 or t > t_end:
            return 0
        elif self.hazard_rate == None:
            return(1)
        else:
            temp = self.hazard_rate(t)
            if temp < 0:
                return 0
            else:
                return temp

    def total_hazard_function(self, t, T_ends):
        '''
        Given a predefined hazard function, an input time, and the lengths of the infectious periods
        The function caculates the rate at which the hazard is emitted at time t, assuming they were all intially infected

        Designed as an input function for the BVP Solver, which doesn't currently work
        '''
        total_hazard_rate = 0
        for i in range(len(T_ends)):
            total_hazard_rate =  total_hazard_rate + self.hazard(t,T_ends[i])
        return total_hazard_rate

    def integrate_hazard(self, t_end):
        '''
        Integrate a hazard function
        '''
        f = lambda t: self.hazard(t,t_end)
        integral = spi.quad(f, 0, t_end)
        return integral[0]
    
    def compute_final_size(self):
        '''Generates 1 observation of the final size of an epidemic from the parameters defined earlier'''
        #Infectious Period variables, function as defined in the __init__ section
        T = self.generate_infection_periods()
        assert all(T) > 0

        #The resistance to infection variables, always an exponential 1 rv
        self.Q = np.random.exponential(scale = 1,
                                  size = self.N)
        
        #We take the formulation found in the Thomas House paper and set the intial infecteds to zero
        self.Q[ :self.inf_starting] = 0
        
        #The ordering
        ordering = np.argsort(self.Q)  #The vector that gives us the correct ordering, incase we want to use it later
        Q_sorted = self.Q[ordering]    #The sorted vector

        if self.hazard_rate == None:
            #Then we assume they emit hazard at a constant rate
            #As such, infection times are exponentially distributed
            cumulative_hazard = self.beta * np.cumsum(T)
        else:
            #The hazard rate is more complicated
            cumulative_hazard = np.empty(self.N)
            for i in range(self.N):
                cumulative_hazard[i] = self.integrate_hazard(T[i])
            cumulative_hazard = self.beta * np.cumsum(cumulative_hazard)
        
        counter = 0
        for i in range(self.N):
            if Q_sorted[i] < cumulative_hazard[i]:
                counter = counter + 1
            else:
                break
        
        return(counter)

    
    def sim_final_size(self, n_sim):
        '''Generates multiple observations of the final size of the epidemic.'''
        self.n_sim = n_sim
        
        self.observations = []
        
        print("Performing", self.n_sim, "iterations!")
        for _ in range(self.n_sim):
            
            new_obs = self.compute_final_size()
            self.observations.append(new_obs)
        
        return self.observations
        
    def plot_hist(self):
        '''Calls the method for generating the observations and then creates the '''
        self.plot = plt.hist(self.observations, bins = range(self.N))
        plt.hist(self.observations, bins = range(self.N), density = True)
        plt.title(f'Final epidemic size of {self.n_sim} observations')
        
