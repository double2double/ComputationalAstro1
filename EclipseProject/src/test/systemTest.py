# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 17:17:48 2014

@author: bob
"""
import system.linearProblem as lnProblem


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import system.function as func
import system.waveSystem as wave

def f(y,t):
    return vgl.f(t, y)



w = 0.01
sigma = 1.
K = 1.
g = 1.

# Create a aid function of the density
class rho0(func.Function):
    def evaluate(self, x):
        return (1+sigma*x)
# Create the two objects to represent the functions P and Q
class P(func.Function):
    def evaluate(self, x):
        return w**2*rho0().evaluate(x)
class Q(func.Function):
    def evaluate(self, x):
        return -K**2*(rho0().evaluate(x)*w**2+rho0().derivative(x)*g)



    
# create the ODE    
vgl = wave.WaveSystem(P(),Q())

# initial condition 
y0 = [0.,1.]
t  = np.linspace(0, 1., 1000)

# solve the DEs
soln = odeint(f, y0, t)
S = soln[:, 0]
Z = soln[:, 1]


 # plot results
plt.plot(t,S)
#plt.plot(t,Z)
#plt.plot(t,R)

plt.show()