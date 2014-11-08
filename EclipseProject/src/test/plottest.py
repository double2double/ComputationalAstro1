'''
Created on 29 Oct 2014

@author: bob
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import function.function as func
import system.waveSystem as wave
import scipy
import matplotlib.pyplot as plot



class MyClass(object):
    '''
    classdocs
    '''

    def __init__(self, data,fig=None):
        self.fig = fig
        for i in range(len(data[:,0])):
            self.plotRow(data[i,:])
    def f(self,y,t):
        return self.vgl.f(t, y)
    def plotRow(self,row):
        Solution = row[3:]
        Ksqr = row[0]
        sigma = row[1]
        g = row[2]
        for i in Solution:
            funcP = func.P(Ksqr,sigma,g,i)
            funcQ = func.Q(Ksqr,sigma,g,i)
            # create the ODE    
            self.vgl = wave.WaveSystem(funcP,funcQ)
            # solve the DEs
            soln = odeint(self.f, y0, t)
            S =soln[:, 0]
            self.fig.plot(t,S)




y0 = [0.,1.]
t  = np.linspace(0, 1., 1000)


        