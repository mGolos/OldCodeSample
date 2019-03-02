#!/usr/bin/env python
import warnings
warnings.filterwarnings("ignore",category=RuntimeWarning)


import sys
sys.path.append('../evaCure/.')
import multiprocessing
from time import time
from pylab import *
from numpy import save, load
from os import mkdir, path


def whatyouwant1(paramin):  # 2D
    out = 0
    for i in range(10):
        prog = "../evaCure/Main.py"
        sys.argv = [prog, "-m",str(paramin), \
                    '-N','64', \
                    #'-m','1', \
                        '-dt','0.01', \
                        '-dur','100', \
                    '-a','-1.', \
                        '-b','1.', \
                    '-theta','0.', \
                        '-nper','1000', \
                        '-err','0.00001', \
                        '-GUI','0', \
                        '-save','0', \
                        '-test','0', \

                    '-netw','0', \
                    '-patt','0', \
                    '-init','1', \
                        '-who','0', \
                        '-fileIN','../FilesIN/initialization.npy', \
                        '-dens','0.5', \
                    '-noise','0.10', \
                    '-conn','0', \
                        '-upd','0', \
                        '-Dx','0.1', \
                        '-posit','0', \
                        '-ap0e0','0', \
                        '-mT','0', \
                        '-p_var','0.', \
                        '-normW','2.', \
                    '-tau_theta','0', \
                    '-theta_ref','0.', \
                        '-dens_moy','0.5', \

                        '-priT','0', \
                        '-resh','1', \
                        '-nbcol','1', \

                        '-thre','0', \
                        '-gamma','30', \
                        '-derr','0.2', \
                        '-tau','1', \
                        '-intFDx','0']
        execfile(prog)
        out += Network.dCorrelation(Network.states)[0]
    out /= 10.
    print paramin, out
    return out

def whatyouwant2(paramin):  # 3D xfois
    import sys
    n2 = 33 ; para2 = linspace(0.05,0.3,n2) # Dx
    out = zeros(n2)

    for j2 in range(n2):
         sys.argv = ["Main.py", "-Dx",string_(para2[j2]), "-m",string_(paramin)]
         execfile("Main.py")
         out[j2] += Network.dCorrelation(Network.states)[0]
         #sys.argv = ["Main.py", "-Dx",string_(para2[j2]), "-m",string_(paramin), "-mT","1", "-upd","2"]
         #execfile("Main.py")
         #out[j2+1*n2] += Network.dCorrelation(Network.states)[0]
    print paramin, out
    return out

def whatyouwant3(paramin): # 4D variation de Dx,m et dens pour adatron
    import sys
    #n2 = 12 ; para2 = linspace(0.08,0.30,n2) # Dx
    para2 = arange(0.001,0.3011,0.01) # Dx
    n2 = len(para2)
    n3 = 22 ; para3 = int0(linspace(1,n3,n3)) # m/N
    out = zeros((n2,n3))

    file2D = 'result2D_%.3f.npy' %paramin
    if path.exists(file2D):
        return load(file2D)
    else:
        for j2 in range(n2):
            file1D = 'result1D_%.3i_%.3f.npy' %(j2 ,paramin)
            if path.exists(file1D):
                out[j2] = load(file1D)
            else:
                open(file1D, 'a').close()
                for j3 in range(n3):
                    sys.argv = ["Main.py", "-dens",string_(paramin), \
                                           "-Dx",string_(para2[j2]), \
                                           "-m",string_(para3[j3]),\
                                '-N','66', \
                                #'-m','1', \
                                '-dt','0.01', \
                                '-dur','100', \
                                '-a','0.', \
                                '-b','1.', \
                                '-theta','0.5', \
                                '-nper','1000', \
                                '-err','0.00001', \
                                '-netw','1', \
                                '-patt','1', \
                                '-init','1', \
                                '-who','0', \
                                #'-dens','0.5', \
                                '-noise','0.05', \
                                '-conn','2', \
                                '-upd','0', \
                                #'-Dx','0.1', \
                                '-posit','0', \
                                '-ap0e0','1', \
                                '-mT','0', \
                                '-p_var','0.', \
                                '-normW','2.', \
                                '-tau_theta','0', \
                                '-priT','0', \
                                '-savS','0', \
                                '-savW','0', \
                                '-savP','0', \
                                '-resh','1', \
                                '-nbcol','1', \
                                '-derr','0.2', \
                                '-tau','1', \
                                '-thre','0', \
                                '-gamma','10', \
                                '-dens_moy','0.3', \
                                '-GUI','0', \
                                '-intFDx','0']
                    execfile("Main.py")
                    #out[j2,j3] += Network.dCorrelation(Network.states)[0]
                    out[j2,j3] += Network.connectivityDensity()
                save(file1D, out[j2])
        save(file2D, out)
        return out

if __name__ == '__main__':
    t0 = time()
    pool = multiprocessing.Pool()
    #pool = multiprocessing.Pool(33)

    args = int0(linspace(1,124,124)).tolist() #m
    results = pool.map(whatyouwant1, args) # launches multiple processes

    #args = int0(linspace(1,33,33)).tolist() #m
    #results = pool.map(whatyouwant2, args) # launches multiple processes

    #args = arange(0.02,0.981,0.03).tolist() # dens
    #results = pool.map(whatyouwant3, args) # launches multiple processes

    save('results.npy', array(results))
    print "Simulation took :",time()-t0,"s"
    ofi = open('Simulation_took_' + str(time()-t0) + 's', 'a')
    ofi.close()