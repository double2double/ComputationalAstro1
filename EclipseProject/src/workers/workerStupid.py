'''
Created on 08 Nov 2014

@author: bob
'''


'''
Created on 08 Nov 2014

@author: bob
'''
from workers.worker import Worker as worker
import scipy
from scipy.optimize import anderson
from scipy.optimize import newton_krylov

GUESS_W = 10
MAX_TOL = 0.000001

class workerStupid(worker):
    '''
    A class to represent a worker to find the eigenmodes as 
    a function of K^2, sigma and g
    '''
    def __init__(self, Ksqr=[1],sigma=[1],g=[1],y0=[0.,1.],n=3,t0=0,tend=1,h=0.01,filename='out.txt',name='HardWorkingStudent'):
        '''
        The constructor to set up the right parameters and to create
        the ode's
        '''
        super(workerStupid, self).__init__(Ksqr, sigma, g, y0, n, t0, tend, h)
        self.f = open(filename, 'w')
        self.filename = filename
        self.name = name
    def search(self,Ksqrnum,sigmanum,gnum,n):
        guess = GUESS_W/0.9
        self.tempKsqrnum = Ksqrnum
        self.tempsigmanum = sigmanum
        self.tempgnum = gnum
        endP = self.endPoint(guess)
        while endP<0:
            guess = guess*10
            endP = self.endPoint(guess)
        info = self.zero_point_info(Ksqrnum,sigmanum,gnum,guess)
        nb_zero = info[0]
        jumpFactor = 1
        guesses = scipy.zeros(n)
        positiveGuess = guess
        for i in range(n):
            teken = (-1)**i
            print i
            jumpFactor = 1
            endP = self.endPoint(guess)
            info = self.zero_point_info(Ksqrnum,sigmanum,gnum,guess)
            nb_zero = info[0]
            shrinksize = 0.9
            overCount = 1
            while scipy.absolute(endP)>MAX_TOL:
                previous = guess
                guess = previous * shrinksize
                endP = self.endPoint(guess)
                print ('endP,PreviousPos,value, shrink = %s ,%s , %s ,%s' )%(guess, positiveGuess, endP,shrinksize)
                info = self.zero_point_info(Ksqrnum,sigmanum,gnum,guess)
                if (scipy.absolute(endP)<MAX_TOL):
                    break
                if (teken*endP<0):
                    #print 'Ooh now'
                    overCount = overCount +1
                    shrinksize = shrinksize + 9./(10**overCount)
                    guess = positiveGuess
                    guess = guess/shrinksize
                elif (teken*endP>0):
                    positiveGuess = guess
            guesses[i] = guess
            while (scipy.absolute(endP)<MAX_TOL ):
                guess = positiveGuess
                guess = guess/shrinksize
                endP = self.endPoint(guess)
                print 'Stepping away: %s'%(endP)
                print guesses
                pass
                
        return guesses
                
            
            
            
            
            
            
            
            