#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore",category=Warning)
from pylab import *
from time import time
from numpy import load as np_load, save as np_save
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

def showshow(W,mi):
    figure()
    abreviations = ['ENT','PARH','TP','FP','FUS','TT','LOCC','SP','IT','IP',\
                    'SMAR','BTST','MT','ST','PSTC','PREC','CMF','POPE','PTRI','RMF',\
                    'PORB','LOF','CAC','RAC','SF','MOF','LING','PCAL','CUN','PARC',\
                    'ISTC','PCUN','PC']
    abre_ind = ['r'+R for R in abreviations]
    abreviations.reverse()
    for i in range(33): abre_ind.append('l'+ abreviations[i])
    yticks(arange(66), abre_ind[:66], fontsize=10)
    imshow(W,interpolation='nearest', vmin=mi)
    colorbar()
    show()

def whatyouwant(paramin):
    global dir_pri2

    dens_moy = arange(0.02,0.981,0.03) # (33)
    W=[]
    for d in range(len(dens_moy)):
        W.append(loadConnectivity2(dir_pri2, paramin, dens_moy[d]))

    return W

if __name__ == '__main__':
    dir_pri2 = './MatricesG30TVarXAda3'
    if not os.path.exists(dir_pri2):
        print 'Pas de dossier ' + dir_pri2
        quit()

    t0 = time()
    pool = multiprocessing.Pool()
    args = arange(0.50,4.501,0.2).tolist() # normW

    results = pool.map(whatyouwant, args)
    results= array(results)
    np_save(dir_pri2 + '/results.npy', results)

    print "Simulation took :",time()-t0,"s"

    W = []
    for j in range(results.shape[1]):
        figure()
        for i in range(results.shape[0]):
            subplot(3,7,i+1)
            title(i)
            imshow(results[i,j],interpolation='nearest')
        show()
        mats = input()
        for i in range(results.shape[0]):
            if i in mats:
                W.append(results[i,j])
    W = array(W)

    showshow(W.sum(0)/len(W),0)

    np_save(dir_pri2 + '/W.npy', array(W))


