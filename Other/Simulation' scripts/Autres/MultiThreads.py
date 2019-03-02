#!/usr/bin/env python

import sys
sys.path.append('../evaCure/.')
import multiprocessing
from FromToFiles import *
from scipy.misc import toimage
from pylab import linspace, array, zeros, arange, corrcoef
from numpy import load as np_load, save as np_save
from time import time
import os


def whatyouwant(paramin):
    global dir_pri
    dir_sub = dir_pri + '/Norm_%.2f' %paramin
    os.mkdir(dir_sub)

    multi_dens = arange(0.02,0.981,0.03) # (33)
    nb_patterns = 0
    tendance =[]
    for d in range(len(multi_dens)):
        for i in range(100):
            prog = "../evaCure/Main.py"
            sys.argv = [prog , "-dens",str(multi_dens[d]), \
                               '-normW',str(paramin), \
                        '-N','66', \
                        '-m','1', \
                        '-dt','0.01', \
                        '-dur','100', \
                        '-a','0.', \
                        '-b','1.', \
                        '-theta','0.5', \
                        '-nper','3000', \
                        '-err','0.00001', \
                        '-GUI','0', \
                        '-save','0', \
                        '-test','1', \

                        '-netw','0', \
                        '-patt','1', \
                        '-init','0', \
                        '-who','0', \
                        '-fileIN','../FilesIN/initialization.npy', \
                    #'-dens','0.5', \
                        '-noise','0.', \
                        '-conn','5', \
                        '-upd','0', \
                        '-Dx','0.1', \
                        '-posit','0', \
                        '-ap0e0','0', \
                        '-mT','0', \
                        '-p_var','0.', \
                    #'-normW','2.', \
                        '-tau_theta','1', \
                        '-theta_ref','0.5', \
                    '-dens_moy','0.26', \

                        '-priT','0', \
                        '-resh','1', \
                        '-nbcol','1', \

                        '-thre','0', \
                        '-gamma','30', \
                        '-derr','0.2', \
                        '-tau','1', \
                        '-intFDx','0']
            execfile(prog)

            maybe = True
            for j in range(1,nb_patterns + 1):
                #if corrcoef(np_load(dir_sub+'/pattern_%.3i.npy'%(j-1)), Network.X)[0,1] > 0.90:
                if corrcoef(np_load(dir_sub+'/pattern_%.3i.npy'%(j-1)), Network.states)[0,1] > 0.90:
                    maybe = False
                    tendance[j-1][d] += 1
                    break
            if maybe:
                #saveStates(dir_sub+'/', Network.X, resh=1, opt=nb_patterns)
                saveStates(dir_sub+'/', Network.states, resh=1, opt=nb_patterns)
                tendance.append(zeros(len(multi_dens)))
                tendance[nb_patterns][d] = 1
                nb_patterns += 1
        print 'norm :%.2f ; dens:%.3f'%(paramin,multi_dens[d])

    ofi = open(dir_sub+'/tendance.txt', 'w')
    for i in range(len(multi_dens)):
        ofi.write('%.2f\t' %multi_dens[i])
        for j in range(len(tendance)):
            ofi.write('%.3i\t' %tendance[j][i])
        ofi.write('\n')
    ofi.close()

    return 1

if __name__ == '__main__':
    dir_pri = './Patterns'
    if not os.path.exists(dir_pri):
        os.mkdir(dir_pri)
    t0 = time()
    pool = multiprocessing.Pool()
    #args = arange(1.65,2.61,0.05).tolist() # normW
    args = arange(0.50,4.501,0.2).tolist() # normW
    #args = [1.90] # normW
    results = pool.map(whatyouwant, args)

    print "Simulation took :",time()-t0,"s"