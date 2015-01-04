'''
Created on 15 Oct 2014

@author: bob

This task tests the worker class in more detail.
'''

from workers.worker2 import Worker2 as worker2
import scipy 
from matplotlib import pyplot as plot

if __name__ == '__main__':
    Ksqr = scipy.linspace(0.5, 1.5, 1).tolist()
    sigma = scipy.linspace(0.5, 1.5,10).tolist()
    g = scipy.linspace(0.5, 1.5, 1).tolist()
    n = 3
    y0=[0.,1.];
    t0=0
    tend=1
    h=0.01
    eigenSys = worker2(Ksqr,sigma,g,y0,n,t0,tend,h)
    data = eigenSys.taskNonParallel()
    print data
    plot.plot(data[:,1],data[:,3])
    plot.show()
    
    

