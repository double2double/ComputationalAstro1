# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 17:17:48 2014

@author: bob
"""

import scipy.linalg

CONVERGENCE = 1
FIXED_NB_ITERATIONS = 2

import nonLinearSolver as nls

class NewtonSolver(nls.NonLinearSolver):
            
    def iterate(self,y):
        """
        takes the current guess and updates it using Newton's method
        Input:
            yk -- current guess of solution
        Output:
            an updated guess of the solution                
        """
        res = self.eq.residual(y)
        J = self.eq.Jacobian(y)
        dy = scipy.linalg.solve(J,-res)
        return y + dy        

