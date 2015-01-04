'''
Created on 29 Oct 2014

@author: bob
'''

import scipy 
from scipy.integrate import odeint
import function.function as func
import system.waveSystem as wave



class PlotWave(object):
    '''
    A class to plot the curves for the given parameters.
    '''
    def __init__(self, data,fig=None):
        self.fig = fig
        for i in range(len(data[:,0])):
            self.plotRow(data[i,:])
    def f(self,y,t):
        return self.vgl.f(t, y)
    def plotRow(self,row):
        y0 = [0.,1.]
        t  = scipy.linspace(0, 1., 1000)
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
            # Normalising the solution
            MAX = max(S[:])
            S = scipy.multiply(1/MAX,S)
            self.fig.plot(t,S)




        