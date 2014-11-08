'''
Created on 15 Oct 2014

@author: bob


A test file for the worker class.

'''

from workers.eigenModes_K_sigma_g import EigenModes_K_sigma_g as eigenMode
import system.waveSystem as wave
import function.function as func
import integrators.rungeKutta as rK
import matplotlib.pyplot as plt

Ksqr = [1]
sigma = [1]
g = [1]
n = 10
y0 = 0
t0 = 1
tend = 1
h = 0.01

#wsqrnum = 0.3-0.285
wsqrnum = 0.015*(1+23./50)
wsqrnum = wsqrnum*(1+7./50)
wsqrnum = 0.001

class rho0(func.Function):
    def __init__(self,Ksqr,sigma,g,wsqr):
        func.Function.__init__(self)
        self.Ksqr = Ksqr
        self.sigma = sigma
        self.g = g
        self.wsqr = wsqr
    def evaluate(self, x):
        return (1+self.sigma*x)
# Create the two objects to represent the functions P and Q
class P(func.Function):
    def __init__(self,Ksqr,sigma,g,wsqr):
        func.Function.__init__(self)
        self.Ksqr = Ksqr
        self.sigma = sigma
        self.g = g
        self.wsqr = wsqr
    def evaluate(self, x):
        return self.wsqr*rho0(self.Ksqr,self.sigma,self.g,self.wsqr).evaluate(x)
class Q(func.Function):
    def __init__(self,Ksqr,sigma,g,wsqr):
        func.Function.__init__(self)
        self.Ksqr = Ksqr
        self.sigma = sigma
        self.g = g
        self.wsqr = wsqr
    def evaluate(self, x):
        return -self.Ksqr*(rho0(self.Ksqr,self.sigma,self.g,self.wsqr).evaluate(x)*self.wsqr+
                           rho0(self.Ksqr,self.sigma,self.g,self.wsqr).derivative(x)*self.g)


def plot_ode(Ksqr, sigma, g, wsqr):
    funcP = P(Ksqr,sigma,g,wsqr)
    funcQ = Q(Ksqr,sigma,g,wsqr)
    vgl = wave.WaveSystem(funcP,funcQ)
    fe = rK.RungeKutta(vgl)
    t_runge,soln_runge = fe.integrate(y0,t0,tend,h)
    solution_runge = [soln_runge[i][0] for i in range(len(soln_runge))]
    plt.plot(t_runge,solution_runge)
    plt.show()
    pass


def nbOfZerosnumber_of_zeros():
    eigenSys = eigenMode(Ksqr,sigma,g,y0,n,t0,tend,h)
    print eigenSys.zero_point_info(Ksqrnum=1, sigmanum=1, gnum=1, wsqrnum=wsqrnum)
    
def find_eigen_mode():
    w = wsqrnum
    eigenSys = eigenMode(Ksqr,sigma,g,y0,n,t0,tend,h)
    nb_of_zero , index_last_min , end_point , length_Set = eigenSys.zero_point_info(Ksqrnum=1, sigmanum=1, gnum=1, wsqrnum=w)
    newW = w
    while (nb_of_zero==0):
        newW = (newW+0.0)/2
        nb_of_zero , index_last_min , end_point , length_Set = eigenSys.zero_point_info(Ksqrnum=1, sigmanum=1, gnum=1, wsqrnum=newW)
    print nb_of_zero , index_last_min , end_point , length_Set
    
    print newW
    previousW = newW
    counter = 0
    while(abs(end_point)>0.001):
        counter = counter +1
        nb_of_zero_new , index_last_min_new , end_point_new , length_Set_new = eigenSys.zero_point_info(Ksqrnum=1, sigmanum=1, gnum=1, wsqrnum=newW)
        print counter, newW , nb_of_zero_new , index_last_min_new , end_point_new , length_Set_new
        if (abs(end_point_new)<0.001):
            break
        if (nb_of_zero_new>=nb_of_zero):
            # Vergroot w
            previousW = newW
            newW = newW*(1+((length_Set-index_last_min_new)+0.0)/index_last_min_new)
        if (nb_of_zero_new < nb_of_zero):
            # Nu weten we dat het vorige punt wel nog achter het nulpunt lag dus nemen we het 
            # gemmidelde tussen het slechte punt en het vorige goede punt.
            newW = ((newW+previousW)+0.0)/2
    print newW
    plot_ode(Ksqr=1, sigma=1, g=1, wsqr=newW)

def task():
    eigenSys = eigenMode(Ksqr,sigma,g,y0,n,t0,tend,h)
    eigenSys.task()
    pass

if __name__ == '__main__':
    Ksqr = [1]
    sigma = [1]
    g = [1]
    n = 10
    y0 = 0
    t0 = 1
    tend = 1
    h = 0.1
    Ksqr=[1];sigma=[1];g=[1];y0=[0.,1.];n=10;t0=0;tend=1;h=0.01
    #nbOfZerosnumber_of_zeros()
    #find_eigen_mode()
    #plot_ode(Ksqr=1, sigma=1, g=1, wsqr=wsqrnum)
    task()



