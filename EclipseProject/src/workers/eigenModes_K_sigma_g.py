'''
Created on 14 Oct 2014

@author: bob
'''
from workers.worker import Worker as worker
import system.waveSystem as wave
import function.function as func
import integrators.rungeKutta as rK
import scipy
from scipy.linalg import norm
from scipy.signal import argrelextrema


# Global parameters to tweak the algoritm
DIST = 1
GUESS_W = 0.0001
UPPERZERO = 0.05




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
        super(EigenModes_K_sigma_g, self).__init__(Ksqr, sigma, g, y0, n, t0, tend, h)
    def task(self):
        '''
        This is the main task. It calculates the first n values for the 
        eigenmodes of w.
        These eigen modes will be placed in a 4 dimensional array so that
        Spectrum[ksqr,sigma,g,:] are the first n spectral lines for the 
        assososieted parameters.
        '''
        print 'started the task'
        # Create a 4 dimentional array to store the values in
        spectrum = scipy.zeros((len(self.Ksqr),len(self.sigma),len(self.g),self.n))
        # Starting the masive three dimtional itteration...
        for i in range(len(self.Ksqr)):
            for j in range(len(self.sigma)):
                for k in range(len(self.g)):
                    print 'Loop loop loop'
                    Ksqr = self.Ksqr[i]
                    sigma = self.sigma[j]
                    g = self.sigma[k]
                    eigen_nodes = scipy.zeros(self.n)
                    # No we are going to exploit the continuity of the solution
                    # and look in the vicinity of the previous points for solutions.
                    # But of course only if there are previous solutions.
                    eigen_nodes = self.search_hard(Ksqr,sigma,g)
                    print ('Eigen Nodes for (%i,%i,%i) = %s')%(i,j,k,eigen_nodes)
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
    def find_previous_eigen_mode(self,Ksqrnum,sigmanum,gnum,wguess):
        '''
        A method that calculates the solution for the eigen modes with w smaller than wguess
        '''
        nb_of_zero , index_last_min , end_point , length_Set = self.zero_point_info(Ksqrnum, sigmanum, gnum, wguess)
        newW = wguess
        while (nb_of_zero==0):
            newW = (newW+0.0)/2
            nb_of_zero , index_last_min , end_point , length_Set = self.zero_point_info(Ksqrnum, sigmanum, gnum, newW)
        previousW = newW
        counter = 0
        while(abs(end_point)>0.001):
            counter = counter +1
            nb_of_zero_new , index_last_min_new , end_point_new , length_Set_new = self.zero_point_info(Ksqrnum, sigmanum, gnum, newW)
            print counter, newW , nb_of_zero_new , index_last_min_new , end_point_new , length_Set_new
            if (abs(end_point_new)<0.001):
                break
            if (nb_of_zero_new>=nb_of_zero):
                # Increase w
                previousW = newW
                newW = newW*(1+((length_Set-index_last_min_new)+0.0)/index_last_min_new)
            if (nb_of_zero_new < nb_of_zero):
                # Nu weten we dat het vorige punt wel nog achter het nulpunt lag dus nemen we het 
                # gemmidelde tussen het slechte punt en het vorige goede punt.
                newW = ((newW+previousW)+0.0)/2
        return newW , nb_of_zero
    def find_next_eigen_mode(self,Ksqrnum,sigmanum,gnum,wguess): 
        'Deze methode gaat hier niet lukken op de trein...'
        '''
        A method that calculates the solution for the eigen modes with w greater than wguess
        '''
        nb_of_zero , index_last_min , end_point , length_Set = self.zero_point_info(Ksqrnum, sigmanum, gnum, wguess)
        newW = wguess
        counter = 0
        while(abs(end_point)>0.001):
            counter = counter +1
            nb_of_zero_new , index_last_min_new , end_point_new , length_Set_new = self.zero_point_info(Ksqrnum, sigmanum, gnum, newW)
            print counter, newW , nb_of_zero_new , index_last_min_new , end_point_new , length_Set_new
            if (abs(end_point_new)<0.001):
                break
            if (nb_of_zero_new == nb_of_zero):
                # Decrease w
                previousW = newW
                newW = newW*(1+((-length_Set+index_last_min_new)+0.0)/index_last_min_new)
            if (nb_of_zero_new < nb_of_zero):
                # Nu weten we dat het vorige punt wel nog achter het nulpunt lag dus nemen we het 
                # gemmidelde tussen het slechte punt en het vorige goede punt.
                newW = ((newW+previousW)+0.0)/2
        return newW , nb_of_zero
    
    
    def search_hard(self,Ksqrnum,sigmanum,gnum):
        '''
        A method to search for the first n eigen modes, the hard way.
        This method will first guess the value of w. And based on the number of 
        modes it will find the value for tune thise guess.
        If it finds a good value for w it will search for the neigbour values 
        for w who are eigenmodes.
        '''
        eigenModes = scipy.zeros(self.n)
        newW = GUESS_W
        while (scipy.count_nonzero(eigenModes)!=self.n):
            newW , nb_of_zero = self.find_previous_eigen_mode(Ksqrnum,sigmanum,gnum,newW)
            print newW , nb_of_zero
            if(nb_of_zero<=(self.n)):
                eigenModes[nb_of_zero-1] = newW
            newW = newW*1.1          
        return eigenModes
        
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
        funcP = func.P(Ksqrnum,sigmanum,gnum,wsqrnum)
        funcQ = func.Q(Ksqrnum,sigmanum,gnum,wsqrnum)
        vgl = wave.WaveSystem(funcP,funcQ)
        fe = rK.RungeKutta(vgl)
        t_runge,soln_runge = fe.integrate(self.y0, self.t0, self.tend, self.h)
        # Now we are going to calculate the local minima of the absolute value of the solution.
        solution_runge = [soln_runge[i][0] for i in range(len(soln_runge))]
        index_local = argrelextrema(scipy.absolute(solution_runge), scipy.less)[0]
        # This will in theory give all the points where the data is zero,
        # plus the starting point and possible also the end point.
        # But due to the possible un smoothness of the data some points could
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
        
        
        
        
        
        
        
        
        