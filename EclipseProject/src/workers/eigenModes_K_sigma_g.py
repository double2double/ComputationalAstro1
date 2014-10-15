'''
Created on 14 Oct 2014

@author: bob
'''
from workers.worker import Worker as worker
import system.waveSystem as wave
import system.function as func
import integrators.rungeKutta as rK
import scipy
from scipy.linalg import norm
from scipy.signal import argrelextrema


# Global parameters to tweak the algoritm
DIST = 1
GUESS_W = 0.1
UPPERZERO = 0.05



class rho0(func.Function):
    def __init__(self,Ksqr,sigma,g,wsqr):
        func.Function.__init__(self)
        self.Ksqr = Ksqr
        self.sigma = sigma
        self.g = g
        self.wsqr = wsqr
    def evaluate(self, x):
        return (1+self.sigma*x)
# Create the two objects to represent the functions P and Q
class P(func.Function):
    def __init__(self,Ksqr,sigma,g,wsqr):
        func.Function.__init__(self)
        self.Ksqr = Ksqr
        self.sigma = sigma
        self.g = g
        self.wsqr = wsqr
    def evaluate(self, x):
        return self.wsqr*rho0(self.Ksqr,self.sigma,self.g,self.wsqr).evaluate(x)
class Q(func.Function):
    def __init__(self,Ksqr,sigma,g,wsqr):
        func.Function.__init__(self)
        self.Ksqr = Ksqr
        self.sigma = sigma
        self.g = g
        self.wsqr = wsqr
    def evaluate(self, x):
        return -self.Ksqr*(rho0(self.Ksqr,self.sigma,self.g,self.wsqr).evaluate(x)*self.wsqr+
                           rho0(self.Ksqr,self.sigma,self.g,self.wsqr).derivative(x)*self.g)



class EigenModes_K_sigma_g(worker):
    '''
    A class to represent a worker to find the eigenmodes as 
    a function of K^2, sigma and g
    '''
    def __init__(self, Ksqr=[1],sigma=[1],g=[1],y0=[0.,1.],n=10,t0=0,tend=1,h=0.01):
        '''
        The constructor to set up the right parameters and to create
        the ode's
        '''
        self.Ksqr = Ksqr
        self.sigma = sigma
        self.g = g
        self.n = n
        self.y0 = y0
        self.t0 = t0
        self.tend = tend
        self.h = h
    def task(self):
        '''
        This is the main task. It calculates the first n values for the 
        eigenmodes of w.
        These eigen modes will be placed in a 4 dimensional array so that
        Spectrum[ksqr,sigma,g,:] are the first n spectral lines for the 
        assososieted parameters.
        '''
        # Create a 4 dimentional array to store the values in
        spectrum = scipy.zeros((len(self.Ksqr),len(self.sigma),len(self.g),self.n))
        # Starting the masive three dimtional itteration...
        for i in range(len(self.Ksqr)):
            for j in range(len(self.sigma)):
                for k in range(len(self.g)):
                    Ksqr = self.Ksqr[i]
                    sigma = self.sigma[j]
                    g = self.sigma[k]
                    eigen_nodes = scipy.zeros(self.n)
                    # No we are going to exploit the continuity of the solution
                    # and look in the vicinity of the previous points for solutions.
                    # But of course only if there are previous solutions.
                    if (i==0 or j==0 or k == 0):
                        eigen_nodes = self.search_hard(Ksqr,sigma,g)
                    if (i !=0 and j != 0 and k != 0):
                        # If there are previous solutions than the distance between
                        # these will first be calculated to get an idea about the
                        # sensibility of these solution to a change in the variables.
                        # The weighted average will be used depending on how far away
                        # the vectors are in param space.
                        dist,average = self._distance(spectrum[i-1,j,k], spectrum[i,j-1,k], spectrum[i,j,k-1],
                                              1./scipy.sqrt(2)*spectrum[i-1,j-1,k],1./scipy.sqrt(2)*spectrum[i-1,j,k-1],1./scipy.sqrt(2)*spectrum[i,j-1,k-1],
                                              1./scipy.sqrt(3)*spectrum[i-1,j-1,k-1])
                        if (dist>self._DIST):
                            eigen_nodes = self.search_hard(Ksqr,sigma,g)
                        else:
                            eigen_nodes= self.search_soft(Ksqr,sigma,g,average)
                    spectrum[i,j,k,:] = eigen_nodes
                    
    def _distance(self,v1,v2,v3,v4,v5,v6,v7):
        """
        A method to calculate the distance between the 7  vectors.
        It returns the distance and the average of the 7 vectors.
        """
        average = 1./7*(v1+v2+v3+v4+v5+v6+v7)
        dist = norm(average-v1, 2)
        dist = dist + norm(average-v2, 2)
        dist = dist + norm(average-v3, 2)
        dist = dist + norm(average-v4, 2)
        dist = dist + norm(average-v5, 2)
        dist = dist + norm(average-v6, 2)
        dist = dist + norm(average-v7, 2)
        dist = 1./7*dist
        return dist,average
    def search_hard(self,Ksqrnum,sigmanum,gnum):
        '''
        A method to search for the first n eigen modes, the hard way.
        This method will first guess the value of w. And based on the number of 
        modes and n will search for the first n eigen modes.
        '''
        eigenModes = scipy.zeros(self.n)
        guess = GUESS_W
        wsqrnum = guess
        nb_of_zero , index_last_min , end_point = self.zero_point_info(Ksqrnum, sigmanum, gnum, wsqrnum)
        '''
        nb_of_zeros will be eighter of three following options
            - 0 => Decrese the guess
            - nb_of_zero < self.n => the guess is in the wright region
            - nb_of_zero > self.n => Increase the guess.
        The first to criteria are good, if it happens to be the 3 we will
        place an other guess till it is one of the two first.
        '''
        while (nb_of_zero > self.n):
            guess = 2 * guess
            nb_of_zero , index_last_min , end_point = self.zero_point_info(Ksqrnum, sigmanum, gnum, wsqrnum)
        # At this point we know we have a good guess and can start finding 
        # the eigen modes.
        while (nb_of_zero==0):
            # Empirisch is er een snelle convergentie naar het nulpunt
            # als we de eindwaarde van w aftrekken.
            guess = scipy.absolute(guess - end_point)
            nb_of_zero , index_last_min , end_point = self.zero_point_info(Ksqrnum, sigmanum, gnum, wsqrnum)
            pass
        # At this point there will be at least one zero point in the function.
        # Lets now start looking for all the values of w for which the numbur
        # of modes is smaller than nb_of_zero.
        nb_of_nodes_to_find = nb_of_zero
        while (nb_of_nodes_to_find != 0):
            
            pass
        
    def zero_point_info(self,Ksqrnum,sigmanum,gnum,wsqrnum):
        '''
        This method will give some information about the zero points of the 
        function.
        Output:
            - nb_of_zero: Holds the number of 0 points
            - index_last_min: Holds the index of the last 0 point
            - value_end_point: Holds the value of the endpoint of the equation
        '''
        #Create the ode
        funcP = P(Ksqrnum,sigmanum,gnum,wsqrnum)
        funcQ = Q(Ksqrnum,sigmanum,gnum,wsqrnum)
        vgl = wave.WaveSystem(funcP,funcQ)
        fe = rK.RungeKutta(vgl)
        t_runge,soln_runge = fe.integrate(self.y0, self.t0, self.tend, self.h)
        # Now we are going to calculate the local minima of the absolute value of the solution.
        solution_runge = [soln_runge[i][0] for i in range(len(soln_runge))]
        index_local = argrelextrema(scipy.absolute(solution_runge), scipy.less)[0]
        # This will in theory give all the points wther the data is zero,
        # plus the starting point and possible also the end point.
        # But due to the possible unsmoothness of the data some points could
        # appear several times.
        # Seen that we are looking for the first n values of we can safely
        # assume that two 0 points should lie at a minimum distance of say
        # ceil(tend/h) + 1)/20 = N/20
        presision = scipy.ceil(self.tend/self.h)/20
        nb_of_zero = 0
        # Count the number real of local 0 points.
        end_point = solution_runge[len(solution_runge)-1]
        if (len(index_local) == 0):
            return 0,0,end_point,len(solution_runge)
        if index_local[0]>presision:
            nb_of_zero = 1
        for i in range(1,len(index_local)):
            if ((index_local[i]-index_local[i-1])<presision):
                # if this is the case the two points are to close to
                # each other and are the same minima.
                pass
            else:
                # In this case we check the number of actual 0
                if(solution_runge[index_local[i]]< UPPERZERO):
                    nb_of_zero = nb_of_zero + 1
        # nb_of_zero should now be the number of times the function was 0  
        index_last_min = index_local[len(index_local)-1]
        return nb_of_zero , index_last_min , end_point , len(solution_runge)
    def search_soft(self,Ksqr,sigma,g,initial_guess):
        '''
        A method to search for the first n eigen modes, the 'easy' way.
        This method will base its searching on a given set of guessed
        values for w
        '''
        pass
        
        
        
        
        
        
        
        
        