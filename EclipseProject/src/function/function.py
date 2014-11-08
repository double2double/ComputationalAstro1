'''
Created on 14 Oct 2014

@author: bob
'''

class Function(object):
    '''
    A class to represent function objects, these functions
    must be able to be evaluated and to be derived at some point.
    '''


    def __init__(self, delta=0.001):
        '''
        Constructor
        '''
        self._DELTA = delta
        pass
    def evaluate(self,x):
        '''
        A method to evaluate the function
        '''
        raise NotImplementedError
    def derivative(self,x):
        '''
        A method to calculate the derivative of the function
        '''
        return (self.evaluate(x+self._DELTA/2) - self.evaluate(x-self._DELTA))/(self._DELTA)


class P(Function):
    def __init__(self,Ksqr,sigma,g,wsqr):
        Function.__init__(self)
        self.Ksqr = Ksqr
        self.sigma = sigma
        self.g = g
        self.wsqr = wsqr
    def evaluate(self, x):
        return self.wsqr*rho0(self.Ksqr,self.sigma,self.g,self.wsqr).evaluate(x)
    
class Q(Function):
    def __init__(self,Ksqr,sigma,g,wsqr):
        Function.__init__(self)
        self.Ksqr = Ksqr
        self.sigma = sigma
        self.g = g
        self.wsqr = wsqr
    def evaluate(self, x):
        return -self.Ksqr*(rho0(self.Ksqr,self.sigma,self.g,self.wsqr).evaluate(x)*self.wsqr+
                           rho0(self.Ksqr,self.sigma,self.g,self.wsqr).derivative(x)*self.g)
        
class rho0(Function):
    def __init__(self,Ksqr,sigma,g,wsqr):
        Function.__init__(self)
        self.Ksqr = Ksqr
        self.sigma = sigma
        self.g = g
        self.wsqr = wsqr
    def evaluate(self, x):
        return (1+self.sigma*x)