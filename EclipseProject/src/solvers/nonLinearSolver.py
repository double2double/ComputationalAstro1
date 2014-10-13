# -*- coding: utf-8 -*-

import scipy

CONVERGENCE = 1
FIXED_NB_ITERATIONS = 2

class NonLinearSolver(object):
        
    def __init__(self, param=None):
        """
        Initializes a nonlinear solver for the equation specified by eq
        Input:
            eq -- an object to which one can ask the information that is 
                required by the specific nonlinear solver (specifically, residual 
                and, if necessary, Jacobian information)
            param -- method parameters (either tolerances or a number of iterations, not both!)
        Output:
            an object of the class NonLinearSolver
        """
        self.eq = None
        if param == None:
            param = NonLinearSolver.getDefaultParameters()
        if param['abs_tol']!=None:
            self.stop_criterion = CONVERGENCE
        elif param['nb_iter']!=None:
            self.stop_criterion = FIXED_NB_ITERATIONS
        else:
            raise ValueError, "Bad combination of input parameters"
        self.param = param
    
    def getDefaultParameters():
        """
        Sets the default parameters of the iterative method
        Output: 
            a structure containing default Parameters for iteration until
            convergence
        """
        param = {}
        param['abs_tol'] = 1e-8
        param['nb_iter'] = None
        return param
    getDefaultParameters = staticmethod(getDefaultParameters)
    
    def iterate(self,yk):
        """
        takes the current guess and updates it
        Beware: this method needs to be overridden for specific methods !
        Input:
            yk -- current guess of solution
        Output:
            an updated guess of the solution                
        """
        raise NotImplementedError
        
    def solve(self,y0):
        """            
        iterate until stopping criterion is satisfied
        Input:
            y0 -- starting value
        Output:
            y -- solution of the nonlinear system             
        """
        # initialise the iterative scheme        
        yk = y0        
        k = 0
        # perform iterations while stopping criterion not satisfied
        while not self.stop(yk,k):
            yk = self.iterate(yk)
            k+=1
        return yk
    
    def stop(self,y,k):
        """
        stopping criterion
        Input:
            y -- current guess of solution
            k -- current iteration number
        Output:
            stop -- Boolean (True if finished, False otherwise)
        """
        if self.stop_criterion==FIXED_NB_ITERATIONS:
            # check if the maximum number of iterations was reached
            cond = (k>=self.param['nb_iter'])
        elif self.stop_criterion==CONVERGENCE:
            # check if the residual is small enough
            cond = (scipy.absolute(self.eq.residual(y))<self.param['abs_tol'])
        return cond

