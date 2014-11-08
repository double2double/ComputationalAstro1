'''
Created on 29 Oct 2014

@author: bob
'''
import multiprocessing
import scipy
from workers.worker2 import Worker2 as worker2
from workers.workerStupid import workerStupid
from workers.worker_nls import WorkerNLS
import matplotlib.pyplot as plot
from test import plottest


if __name__ == '__main__':
    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    jobs = []
    nb_ofLoopsPerProces = 2
    nb_Proces = 8
    Ksqr = scipy.linspace(0.1, 10,  nb_ofLoopsPerProces*nb_Proces).tolist()
    sigma = scipy.linspace(1, 1.5, 1).tolist()
    g = scipy.linspace(1, 1.5, 1).tolist()
    n = 6
    y0=[0.,1.];
    t0=0
    tend=1
    h=0.01
    for i in range(nb_Proces):
        name = 'Worker %s'%(i+1)
        filename = 'File%s.txt'%(i+1)
        eigenSys = workerStupid(Ksqr[nb_ofLoopsPerProces*i:nb_ofLoopsPerProces*(i+1)],sigma,g,y0,n,t0,tend,h,filename,name)
        p = multiprocessing.Process(target=eigenSys.task,args=(i+1,return_dict))
        jobs.append(p)
        p.start()
    i = 1; 
    for p in jobs:
        p.join()
        print 'Job %s finished, from %s'% (i,len(jobs))
        i+=1
    result = scipy.zeros((nb_ofLoopsPerProces*nb_Proces,n+3))
    print return_dict
    for i in range(nb_Proces):
        result[i*nb_ofLoopsPerProces:(i+1)*nb_ofLoopsPerProces,:] =  return_dict[i+1]
    ax1 = plot.subplot2grid((3,4), (0,0), colspan=1)
    ax2 = plot.subplot2grid((3,4), (1,0), colspan=1)
    ax3 = plot.subplot2grid((3,4), (2, 0), colspan=1)
    ax4 = plot.subplot2grid((3,4), (0, 1), colspan=1)
    ax5 = plot.subplot2grid((3,4), (1, 1), colspan=1)
    ax6 = plot.subplot2grid((3,4), (2, 1), colspan=1)
    ax7 = plot.subplot2grid((3,4), (0, 2), rowspan=3,colspan=2)

    
    ax1.plot(result[:,3],result[:,0])
    ax2.plot(result[:,4],result[:,0])
    ax3.plot(result[:,5],result[:,0])
    ax4.plot(result[:,6],result[:,0])
    ax5.plot(result[:,7],result[:,0])
    ax6.plot(result[:,8],result[:,0])
    plottest.MyClass(data=result,fig = ax7)
    plot.savefig('../../plot/mainResult.eps')
    plot.show()
    
    
    
        

    