'''
Created on 08 Nov 2014

@author: bob
'''
from workers.worker import Worker as worker
import scipy
from scipy.optimize import anderson
from scipy.optimize import newton_krylov

GUESS_W = 1

class WorkerNLS(worker):
    '''
    A class to represent a worker to find the eigenmodes as 
    a function of K^2, sigma and g
    '''
    def __init__(self, Ksqr=[1],sigma=[1],g=[1],y0=[0.,1.],n=3,t0=0,tend=1,h=0.01,filename='out.txt',name='HardWorkingStudent'):
        '''
        The constructor to set up the right parameters and to create
        the ode's
        '''
        super(WorkerNLS, self).__init__(Ksqr, sigma, g, y0, n, t0, tend, h)
        self.f = open(filename, 'w')
        self.filename = filename
        self.name = name
    def search(self,Ksqrnum,sigmanum,gnum,n):
        '''
        Search for the first n roots.
        '''
        guess = GUESS_W
        self.tempKsqrnum = Ksqrnum
        self.tempsigmanum = sigmanum
        self.tempgnum = gnum
        Nroots =  scipy.zeros(n)
        allRoots = scipy.zeros(100)
        goal = 0
        if n ==1:
            orde = 1
            while(orde!=0):
                oneOnGuess = 1./guess
                print 'Root guess: %s'%(1./oneOnGuess)
                root = (1./scipy.absolute(anderson(self.aid_f, oneOnGuess,verbose=True,f_tol=6e-5)))
                info = self.zero_point_info(Ksqrnum,sigmanum,gnum,root)
                orde = info[0]
                print 'Root and order: %s , %s'%(root,orde)
                guess = root*(orde + 1)
        else:
            pass
        
        
        return [root]
    def aid_f(self,oneonguess):
        return self.endPoint(scipy.absolute(1/oneonguess))
    
    def _foundAll(self,Nroots):
        for root in Nroots:
            if (root==0):
                return False
        return True
    def _nextGuess(self,Nroots,guess,previous=None):
        '''
        A method to do an educated guess for the next omega value
        '''
        # Find index of lowest 0
        index = 0
        for root in Nroots:
            if (root==0):
                break
            index = index +1
        print '-'*40
        print previous
        print index
        print '-'*40
        if (previous!=0 and previous==index):
            return (guess/2,index)
        newguess = 0
        count = 1
        nonzero = 0
        for root in Nroots:
            newguess = newguess + root/count
            count = 1 + count
            if (root!=0):
                nonzero = nonzero +1
        if nonzero!=0:
            newguess  = (newguess/nonzero)/(index+1)
            return (newguess,index)
        else:
            return (guess/2,index)
        
        
        
        
 
        
        
        
        
        
        
        
        
        
