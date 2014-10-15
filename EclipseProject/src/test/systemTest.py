# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 17:17:48 2014

@author: bob
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import system.function as func
import system.waveSystem as wave
from scipy import zeros
import scipy

def f(y,t):
    return vgl.f(t, y)

Solution = [ 0.1057867,   0.02557878,  0.01131672 ,
             0.00634341,  0.00407022,  0.0028213,
              0.00207537 , 0.00158514 , 0.00124911 , 0.00101265]



# Create a aid function of the density
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

wsqr = 0.00101265
sigma = 1.
Ksqr = 1.
g = 1.
# initial condition 
y0 = [0.,1.]
t  = np.linspace(0, 1., 1000)

for i in Solution:

    funcP = P(Ksqr,sigma,g,i)
    funcQ = Q(Ksqr,sigma,g,i)
    
    # create the ODE    
    vgl = wave.WaveSystem(funcP,funcQ)
    # solve the DEs
    soln = odeint(f, y0, t)
    S =soln[:, 0]
    plt.plot(t,S)

# plot results
plt.plot(t,S)
#plt.plot(t,Z)
#plt.plot(t,R)

plt.show()