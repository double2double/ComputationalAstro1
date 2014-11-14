'''
Created on 13 Oct 2014

@author: bob
'''
import ode_system
class WaveSystem(ode_system.ODESystem):
    
    def __init__(self,P=None,Q=None):
        """
        Creates an object to represent a differential equation of the form:
            d/dx[p(x,w) d/dx e(x)] - q(x,w) e(x) = 0
        This equation is internally converted to a system of first order differential equations
            x'_1 = x_2
            x'_2 = (q(x,w) x_1  - p(x,w) x_2)/(p(x,w)')
            
        Input:
            p -- an object of type function to represent p(x,w) in the equation
            q -- an object of type function to represent q(x,w) in the equation
        Output:
            an object of the class WaveSystem 

        """
        self.P = P
        self.Q = Q

    def f(self,x,y):
        """
        Return the righthand side of the ODE
        """
        P = self.P
        Q = self.Q 
        dy_1 = y[1]/P.evaluate(x)
        dy_2 = Q.evaluate(x)*y[0]
        return [dy_1, dy_2]
        
    
    
    
    
