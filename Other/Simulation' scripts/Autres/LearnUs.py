#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore",category=Warning)
from pylab import *
from time import time
from numpy import load as np_load
import multiprocessing
import commands
import os
import sys
sys.path.append('../evaCure/.')
from FromToFiles import *


def saveConnectivity2(dirOUT,vector,f1,f2):
    '''Save 2D vector as "connectivity.dat"'''
    ofi = open(dirOUT + '/connectivity%.2f_%.2f.dat'%(f1,f2), 'w')
    for i in range(vector.shape[0]):
        for j in range(vector.shape[0]):
            ofi.write(str(vector[i,j]) + '\t')
        ofi.write('\n')
    ofi.close()

def loadConnectivity2(dirIN,f1,f2):
    '''Define connectivity from file in FolderIN folder'''
    ofi = open(dirIN + '/connectivity%.2f_%.2f.dat'%(f1,f2), 'r')
    W = array([line.split() for line in ofi]).astype(float)
    ofi.close()
    return W

def whatyouwant(paramin):
    global dir_pri, dir_pri2

    dens_moy = arange(0.02,0.981,0.03) # (33)
    for d in range(len(dens_moy)):
        dir_sub = '/Norm_%.2f_DensMoy_%.2f' %(paramin,dens_moy[d])
        if not os.path.exists(dir_pri + dir_sub):
            print 'Pas de sous-dossier ' + dir_sub
            return 0

        nb_de_candidats = len(os.listdir(dir_pri + dir_sub)) / 2
        dirIN = dir_pri + dir_sub

        prog = "../evaCure/Main.py"
        sys.argv = [prog ,  "-dens_moy",str(dens_moy[d]), \
                            '-normW',str(paramin), \
                            '-direIN', str(dirIN + '/'), \
                            '-m', str(nb_de_candidats), \
                    '-N','66', \
                    '-dt','0.01', \
                '-dur','0.01', \
                    '-a','0.', \
                    '-b','1.', \
                    '-theta','0.5', \
                    '-nper','3000', \
                    '-err','0.00001', \
                    '-GUI','0', \
                    '-save','0', \
                    '-test','2', \

                    '-netw','1', \
                '-patt','3', \
                    '-init','3', \
                    '-who','0', \
                    '-fileIN','../FilesIN/initialization.npy', \
                    '-dens','0.5', \
                    '-noise','0.', \
                '-conn','3', \
                '-upd','1', \
                '-Dx','0.1', \
                '-posit','1', \
                '-ap0e0','1', \
                    '-mT','0', \
                    '-p_var','0.', \
                #'-normW','2.', \
                    '-tau_theta','1', \
                    '-theta_ref','0.5', \
                #'-dens_moy','0.3', \

                    '-priT','0', \
                    '-resh','1', \
                    '-nbcol','1', \

                    '-thre','0', \
                    '-gamma','30', \
                    '-derr','0.2', \
                    '-tau','1', \
                    '-intFDx','0']
        execfile(prog)

        saveConnectivity2(dir_pri2, Network.W, paramin, dens_moy[d])
        print 'Norm %.2f, DensMoy %.2f, maxW %.2f' %(paramin,dens_moy[d],Network.W.max())
    return 1

if __name__ == '__main__':
    dir_pri  = './Ca_x_Gamma_30'
    dir_pri2 = './MatricesG30TVarXAda3'

    if not os.path.exists(dir_pri):
        print 'Pas de dossier ' + dir_pri
        quit()
    if not os.path.exists(dir_pri2):
        os.mkdir(dir_pri2)

    t0 = time()
    pool = multiprocessing.Pool()
    args = arange(0.50,4.501,0.2).tolist() # normW
    #args = [1.3] # normW
    results = pool.map(whatyouwant, args)

    print "Simulation took :",time()-t0,"s"
