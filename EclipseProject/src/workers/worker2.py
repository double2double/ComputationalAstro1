'''
Created on 27 Oct 2014

@author: bob
'''

from workers.worker import Worker as worker
import scipy
from scipy.signal import argrelextrema



# Global parameters to tweak the algorithm
DIST = 1
GUESS_W = 10
UPPERZERO = 0.05
PRECISION = 400



class Worker2(worker):
    '''
    A class to represent a worker to find the eigenmodes as 
    a function of K^2, sigma and g
    '''
    def __init__(self, Ksqr=[1],sigma=[1],g=[1],y0=[0.,1.],n=3,t0=0,tend=1,h=0.01,filename='out.txt',name='HardWorkingStudent'):
        '''
        The constructor to set up the right parameters and to create
        the ode's
        '''
        super(Worker2, self).__init__(Ksqr, sigma, g, y0, n, t0, tend, h)
        self.f = open(filename, 'w')
        self.filename = filename
        self.name = name
                        
    def search(self,Ksqrnum,sigmanum,gnum,n):
        '''
        A method to search for the first n eigen modes, the hard way.
        This method will first guess the value of w. And based on the number of 
        modes it will find the value for tune thise guess.
        If it finds a good value for w it will search for the neigbour values 
        for w who are eigenmodes.
        '''
        self.tempKsqrnum = Ksqrnum
        self.tempsigmanum = sigmanum
        self.tempgnum = gnum
        N = 100      
        #print w
        fx = scipy.zeros(N)
        nb_of_eigen = 0
        count = 0
        solutionW = []
        while nb_of_eigen < n:
            count = count + 1
            w = scipy.logspace(0.1*count,10*count,N,base=0.5)
            count = count+1
            for x in xrange(N):
                fx[x] = self.endPoint(w[x])
            index_local = argrelextrema(scipy.absolute(fx), scipy.less)[0]
            solutW = scipy.zeros(len(index_local))
            #matplotlib.pyplot.show(block=False)
            for i in xrange(len(index_local)):
                solutW[i] = w[index_local[i]] 
            nb_of_eigen = len(index_local) + nb_of_eigen
            solutionW = scipy.append(solutionW, solutW, 1)
            #matplotlib.pyplot.draw()
        #matplotlib.pyplot.show(block=False)
        return solutW[0:n]
        
    def getAnswer(self):
        if (self.spectrum==None):
            raise RuntimeError('No spectrum calculated')
        return self.spectrum
        



