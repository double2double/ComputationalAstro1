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
        