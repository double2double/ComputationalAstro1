# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 17:17:48 2014

@author: bob
"""

import scipy

import ode_system

class linearProblem(ode_system.ODESystem):
    
    def __init__(self,a,b,c):
        """ 
        creates an system of differential equations of the form:
            x_1 = ((b^2+c^2+2ab) x'_1 + (2b+a) x'_2 + x'_3 + a(b^2+c^2))/(a(b^2+c^2))
            x_2 = x'_1
            x_3 = x'_2
          or in the form x' = Ax + c
            x'_1 = x_2
            x'_2 = x_3
            x'_3 = a(b^2+c^2) x_1 - (b^2+c^2+2ab) x_2 - (2b+a) x_3 - a(b^2+c^2)
          or
            
        Input:
            a -- one of the parameters of the differential equations
            b -- one of the parameters of the differential equations
            c -- one of the parameters of the differential equations
        Output:
            an object of the class linearProblem 

        """
        self.a = a
        self.b = b
        self.c = c
    
    def f(self,t,y):
        """ 
        contains the righthand side of a linear scalar ODE 
        of the form y'=f(t,y) as specified above
        Input:
            t -- current time
            y -- current state, this is a 3 by 1 vector
        Output:
            righthand side f(t,y) 
        
        """
        a = self.a
        b = self.b
        c = self.c
        dy_1 = y[1]
        dy_2 = y[2]
        dy_3 = a*(b**2+c**2) * y[0] - (b**2+c**2+2*a*b) * y[1] - (2*b+a) * y[2] - a*(b**2+c**2)
        return [dy_1,dy_2,dy_3]
        
    
    def Jacobian(self,t,y):
        """ 
        Returns the jacobian of the system, in this case it will
        be a 3 by 3 matrix. For this simple linear system it is 
        just the matrix A as described above in the init.
        Input:
            t -- current time
            y -- current state, this is a 3 by 1 vector
        Output:
            Jacobian matrix of the system of differential equations. 
        """
        a = self.a
        b = self.b
        c = self.c
        # Setting up the Jacobian
        jacobian = [[0 * j *i for i in range(3)] for  j in range(3)]
        jacobian[1][2] = 1
        jacobian[2][3] = 1
        jacobian[3][1] = a*(b^2+c^2)
        jacobian[3][2] = -(b^2+c^2+2*a*b) 
        jacobian[3][3] = (2*b+a)
        
        return jacobian
    def jac(self,t,y):
        return self.Jacobian(t, y)
    
    def y_exact(self,t,t0,y0):
        """ 
        contains the exact solution of the logistic equation  
        Input:
            time -- array of values at which time is required
            y0 -- initial value
            t0 -- initial time
        Output:
            exact solution y(time)
        """
        y = scipy.zeros(len(t))
        for i in range(len(t)):
            y[i] = 1 + scipy.exp(-self.b*t[i])*scipy.sin(self.c*t[i]) + scipy.exp(-self.a*t[i])
        print y
        return y

    

