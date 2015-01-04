'''
Created on 29 Oct 2014

@author: bob
'''
import multiprocessing
import scipy
from workers.worker2 import Worker2 as worker2
from workers import workerSimple as workerSimple
from workers.worker_nls import WorkerNLS
import matplotlib.pyplot as plot
from test import plottest


if __name__ == '__main__':
    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    jobs = []
    nb_ofLoopsPerProces = 10
    nb_Proces = 8
    Ksqr = [1]
    sigma = [1]
    g = scipy.linspace(-0.1, -500,  nb_ofLoopsPerProces*nb_Proces).tolist()
    n = 10
    y0=[0.,1.];
    t0=0
    tend=1
    h=0.01
    for i in range(nb_Proces):
        name = 'Worker %s'%(i+1)
        filename = 'File%s.txt'%(i+1)
        eigenSys = workerSimple(Ksqr,sigma,g[nb_ofLoopsPerProces*i:nb_ofLoopsPerProces*(i+1)],y0,n,t0,tend,h,filename,name)
        p = multiprocessing.Process(target=eigenSys.task,args=(i+1,return_dict))
        jobs.append(p)
        p.start()
    i = 1; 
    for p in jobs:
        p.join()
        print 'Job %s finished, from %s'% (i,len(jobs))
        i+=1
    result = scipy.zeros((nb_ofLoopsPerProces*nb_Proces,n+3))
    for i in range(nb_Proces):
        result[i*nb_ofLoopsPerProces:(i+1)*nb_ofLoopsPerProces,:] =  return_dict[i+1]
    fig = plot.figure()
    fig.suptitle('Dispersion Relation', fontsize=18)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.85)
    ax.set_xlabel('g', fontsize=16)
    ax.set_ylabel('Omega sqr', fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=14)
    for i in range(n):
        ax.plot(result[:,2],result[:,3+i])
    plot.savefig('../../plot/dispersionG.eps')

    plot.show()

    
    
        

    