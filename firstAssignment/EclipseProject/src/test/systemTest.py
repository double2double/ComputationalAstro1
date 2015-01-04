# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 17:17:48 2014

@author: bob

A class to visually check the solutions of the differential equations.

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import system.waveSystem as wave
import function.function as func

# A dummy function for the runge kutta solver.
def f(y,t):
    return vgl.f(t, y)
# Values for wsqr
Solution = [ 0.1 ]
# Parameters of the differential equations
sigma = 1
Ksqr = 1
g = -1
# initial condition 
y0 = [0.,1.]
t  = np.linspace(0, 1., 1000)
# Start the calculation of the ode for the differert values of wsqr
for i in Solution:
    funcP = func.P(Ksqr,sigma,g,i)
    funcQ = func.Q(Ksqr,sigma,g,i)
    # create the ODE    
    vgl = wave.WaveSystem(funcP,funcQ)
    # solve the DEs
    soln = odeint(f, y0, t)
    S =soln[:, 0]
    plt.plot(t,S)
# plot results
plt.plot(t,S)
plt.show()





