#!/usr/bin/env python
#-*- coding:Utf-8 -*-

############################################################################################'
# **Libraries**
import sys, os.path
sys.path.append( os.path.expanduser('/home/mathieuG') )

from evaCure import main
from Tools.ext import Pdir, array2data, mkdir, data2array, adressExists, paramExplo
from Tools.functions import similarity, sortBy, preClustering, similarity_Euclidean, fPearsonCorrelation
from pylab import zeros, arange, rand
import os


############################################################################################'
# **Parameters**
while len(sys.argv) > 1:
    option = sys.argv[1];                                                del sys.argv[1]
    if   option == '-G'  :    G = float(sys.argv[1].replace(',','.'));   del sys.argv[1]
    elif option == '-revert': revert = int(sys.argv[1]);                 del sys.argv[1]
    elif option == '-dir':    dir_= str(sys.argv[1]);                    del sys.argv[1]
    elif option == '-con':    dco = str(sys.argv[1]);                    del sys.argv[1]
    else:
        print 'Options invalides :',option,'->',sys.argv[0]

mdir = '/home/mathieuG/SearchForPatt/%s' %dir_
dir_prix = mdir + '/Potentiels'
dir_priA = mdir + '/Activity'
dir_priT = mdir + '/Tendances'
dir_priC = mdir + '/Criticality'
mmax = 3300
sim_coef = 0.9
multi_dens = arange(0.02,0.981,0.03)


############################################################################################'
# **Function**
def funcNbofAttr(P, G):

    endDir = 'G_%.3f_P_%.5f.npy' %(G,P) 
    try:
        pattsx = data2array(dir_prix + '/patterns_' + endDir, mmap_mode="r+")
        pattsA = data2array(dir_priA + '/patterns_' + endDir, mmap_mode="r+")
        return len(pattsx)

    except:
        try:
            pattsx = data2array(dir_prix + '/allPatt_' + endDir)
            pattsA = data2array(dir_priA + '/allPatt_' + endDir)
                    
        except:
            N = data2array(dco).shape[0]
            pattsx = zeros((mmax,N))
            pattsA = zeros((mmax,N))

            conn = {'connAd': dco,
                    'normType': '1'}

            noise = {'colors': None}

            model = {'model': 'HopfieldBasedStatic',
                    'threshold': 'local',
                    'tauT': 0,
                    'P': P,
                    'G': G} 

            out = []

            other = {'init': 'rand',
                    'dens': rand(),} #p2,!!!!!!!!!!!!!!  RAND

            for d in range(mmax):
                eva = main.evaCure(evaCon=conn, evaNoi=noise, evaMod=model, out=out, **other)
                eva.toEquilibrium()
                pattsx[d] = eva.evaMod.x.copy()
                pattsA[d] = eva.evaMod.A.copy()
                 
            array2data(pattsx, dir_prix + '/allPatt_' + endDir)
            array2data(pattsA, dir_priA + '/allPatt_' + endDir)
        
        patts = pattsx  # !!!!!!!!!
        S = sortBy(patts.mean(1) - patts.mean(), inverse=1)[0]
        C1, freq = preClustering(patts[S],                sim_coef=sim_coef, sim_func=similarity_Euclidean)
        C2, freq = preClustering(patts[S][C1], freq=freq, sim_coef=sim_coef, sim_func=fPearsonCorrelation)
        SC, freq = sortBy(freq, inverse=1)

        array2data(pattsx[S][C1][C2][SC], dir_prix + '/patterns_' + endDir)
        array2data(pattsA[S][C1][C2][SC], dir_priA + '/patterns_' + endDir)
        array2data(freq, dir_priT + '/tendances_' + endDir)
        os.system('rm ' + dir_prix + '/allPatt_' + endDir)
        os.system('rm ' + dir_priA + '/allPatt_' + endDir)
    
        return len(pattsx[S][C1][C2][SC])


############################################################################################'
# **Main**
dir_sub = '/dP_G_%.3f_r_%i.npy' %(G, revert)

if not adressExists(dir_priC + dir_sub):
    mkdir(dir_prix)
    mkdir(dir_priA)
    mkdir(dir_priT)
    mkdir(dir_priC)
    
    if revert: x = [0.95, 2.00]
    else:      x = [0.00, 1.05]
    dx = paramExplo(funcNbofAttr, nb=[7,1], ax='x', x=x, y=[G], revert=revert)
    array2data(dx, dir_priC + dir_sub)
    
    
    