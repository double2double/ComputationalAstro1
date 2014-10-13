# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 17:17:48 2014

@author: bob
"""



class ODESystem(object):
    
    def __init__(self):
        """
        Constructor for ODE system.  Input can take the form that the 
        subclasses desire. This abstract class is meant to create a common 
        interface that can be used by Integrator objects.
        """
        raise NotImplementedError

    def f(self,t,y):
        """
        Return the righthand side of the ODE
        """
        raise NotImplementedError        
        
    def y_exact(self,t,y0):
        """
        Returns the exact solution for the ODE at times t, 
        starting from y0 at time t=0
        """
        raise NotImplementedError
        
