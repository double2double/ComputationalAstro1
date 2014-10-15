'''
Created on 15 Oct 2014

@author: bob

This task tests the worker class in more detail.
'''
from workers.eigenModes_K_sigma_g import EigenModes_K_sigma_g as eigenMode
import system.waveSystem as wave
import system.function as func
import integrators.rungeKutta as rK
import matplotlib.pyplot as plt
import scipy 

if __name__ == '__main__':
    Ksqr = scipy.linspace(0.5, 1.5, 5).tolist()
    sigma = scipy.linspace(0.5, 1.5, 5).tolist()
    g = scipy.linspace(0.5, 1.5, 5).tolist()
    n = 10
    y0 = 0
    t0 = 1
    tend = 1
    h = 0.1
    y0=[0.,1.];n=10;t0=0;tend=1;h=0.01
    eigenSys = eigenMode(Ksqr,sigma,g,y0,n,t0,tend,h)
    eigenSys.task()
    
    
    
    
    
    # Setting up the values where we need to find w
    '''
    Ksqr = scipy.linspace(0.5, 1.5, 5).tolist()
    sigma = scipy.linspace(0.5, 1.5, 5).tolist()
    g = scipy.linspace(0.5, 1.5, 5).tolist()
    n = 10
    Ksqr = [1]
    sigma = [1]
    g = [1]
    n = 10
    # This will in total calculate 1250 values for w
    t0 = 1
    y0=[0.,1.]
    tend = 1
    tend=1
    h=0.01
    '''

