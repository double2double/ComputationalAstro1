# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 17:17:48 2014

@author: bob
"""

import scipy
import integrators.multiStepIntegrator as msi
from scipy.io.matlab.mio5_utils import scipy
from fpformat import sci

class AdamBashford4(msi.MultiStepIntegrator):

    def __init__(self, ode):
        """
        Initializes a forward Euler time integration object for the system of
        ODEs specified by the function rhs.
        Input:
            ode -- an object of the class ODESystem that contain the ODE to be
            integrated
        Output:
            an object of the class ForwardEulerIntegrator that can integrate
            the ODE y'=f(t,y)
        """
        self.ode = ode
        self.s = 4

    def step(self, tn, yn, h,fn):
        """
        A method to calculate the next step using the adams bashford
        method. This method is based on y_(n+s) = y_(n+s-1) + sum(b_i*f_i)
        Input:
            tn -- A vector of length 4 having the previous times 
            yn -- state at times t_n, t_{n+1}, ..., t_{n+s}
            fn -- time derivatives at times t_n ... t_{n+s}
            h - size of time step
        Output:
            y -- state at time t_{h+s+1}
            f -- time derivative at time t_{n+s+1} 
        """
        #Generate the vector for b, this is a fixed vector for each step.
        b = [1 , 1 , 1 , 1 ]
        lincombinatie = 0
        for i in range(self.s):
            prod = scipy.multiply(b[i],fn[i])
            lincombinatie = scipy.add(lincombinatie,prod)
        lincombinatie = scipy.multiply(h,lincombinatie)
        return scipy.add(lincombinatie, yn[self.s-1])
    def y_exact(self,t,t0,y0):
        y = scipy.zeros(self.s)
        for i in range(self.s):
            y[i] = 1 + scipy.exp(-self.b*t[i])*scipy.sin(self.c*t[i]) + scipy.exp(-self.a*t[i])
        print y
        return y

        
    def integrate(self,y0,t0,tend,h):
        """
        Integrates over a number of time steps
        Input:
            t0 -- initial time 
            y0 -- initial condition at time t0
            tend -- time horizon of time integration
            h -- size of time step
        """
        # obtain the number of time steps
        N = int(scipy.ceil(tend/h))
        # create a vector of time instances 
        t = scipy.arange(t0,N*h+h/2.,h)
        # obtain the number of equations
        D = scipy.size(y0)
        # create the matrix that will contain the solutions
        y = scipy.zeros((N+1,D))
        # create a matrix that will contain all time derivatives
        # to avoid re-evaluation of time derivatives
        f = scipy.zeros((N+1,D))
        # perform the time integration
        s = self.s
        # put starting values to exact solution
        print y[0:s,:]
        print self.ode.y_exact(t[0:s],t0,y0)
        y[0:s,:]=self.ode.y_exact(t[0:s],t0,y0)
        for m in range(s):
            f[m,:]=self.ode.f(t[m],y[m,:])
        for n in range(N-s+1):
            # take the next step (and also obtain the time derivative in 
            # the next point)
            print y[n+s,:],f[n+s,:]
            y[n+s,:],f[n+s,:]=self.step(t[n],y[n:n+s,:],h,f[n:n+s,:])
        return t,y

