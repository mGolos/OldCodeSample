#!/usr/bin/env python
#-*- coding:Utf-8 -*-

############################################################################################'
# **Libraries**
import sys, os.path
sys.path.append( os.path.expanduser('/home/mathieuG') )

from Tools.ext import array2data, mkdir, findParameters, data2array, adressExists
from Tools.display import tic, tac
from Tools.functions import sortBy, preClustering, similarity_Euclidean, fPearsonCorrelation
import os


############################################################################################'
# **Parameters**
while len(sys.argv) > 1:
    option = sys.argv[1];                                             del sys.argv[1]
    if   option == '-p1':  p1 = float(sys.argv[1].replace(',','.'));  del sys.argv[1]
    elif option == '-p2':  p2 = float(sys.argv[1].replace(',','.'));  del sys.argv[1]
    elif option == '-p3':  p3 = str(sys.argv[1]);                     del sys.argv[1]
    else:
        print 'Options invalides :',option,'->',sys.argv[0]

main = '/home/mathieuG/SearchForPatt/%s' %p3
typ2 = '_pC'
mkdir(main + typ2)
dir1S = main + '/Cand_S/'
dir1X = main + '/Cand_X/'
dir2S = main + typ2 + '/Cand_S/'
dir2X = main + typ2 + '/Cand_X/'
mkdir(dir2S)
mkdir(dir2X)
sim_coef = 0.9
#endFName = '_P_%.2f_f0_%.2f.npy' %(p1,p2)
endFName = '_P_%.2f_G_%.2f.npy' %(p1,p2)
fileName = 'patterns' + endFName


############################################################################################'
# **Main**
if not adressExists(dir2X + fileName):
    tic()
    pattsX = data2array(dir1X + fileName, mmap_mode="r+")#.reshape((3300,998))
    pattsS = data2array(dir1S + fileName, mmap_mode="r+")#.reshape((3300,998))
    patts = pattsX  #!!!

    S = sortBy(patts.mean(1) - patts.mean(), inverse=1)[0]
    C1, freq = preClustering(patts[S],                sim_coef=sim_coef, sim_func=similarity_Euclidean)
    C2, freq = preClustering(patts[S][C1], freq=freq, sim_coef=sim_coef, sim_func=fPearsonCorrelation)
    SC, freq = sortBy(freq, inverse=1)

    array2data(pattsX[S][C1][C2][SC], dir2X + fileName)
    array2data(pattsS[S][C1][C2][SC], dir2S + fileName)
    array2data(freq, dir2X + '/tendances' + endFName)
    tac()
