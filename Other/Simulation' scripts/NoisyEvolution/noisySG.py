#!/usr/bin/env python
# -*- coding: utf-8 -*-


############################################################################################'
# **Libraries**
import sys, os.path
sys.path.append( os.path.expanduser('/home/mathieuG') )

from evaCure import main
from evaCure import parameters
from Tools.ext import data2array, Pdir, array2data, adressExists, findParameters
from Tools.display import tic, tac, mapMatrices
from pylab import randint, array, concatenate, zeros


############################################################################################'
# **Parameters**
while len(sys.argv) > 1:
    option = sys.argv[1];                                           del sys.argv[1]
    if option == '-p1':    p1 = int(sys.argv[1].replace(',','.'));  del sys.argv[1]
    else:
        print 'Options invalides :',option,'->',sys.argv[0]

if not adressExists('./SG/x_%i.npy'%p1):

    conn = {'connAd': Pdir('Connectomes/SC_FB_D_998_0.npy'),
            'normType': '1'}
    
    noise = {'stdD_x': 0.0,
             'stdD_T': 0.0,
             'colors': ['white','white']}

    model = {'model': 'HopfieldBasedStatic',
            'threshold': 'global',
            'tauT': 80,
            'P': 0.86,
            'G': 900.}

    out = ['x', 'A']

    other = {'init': 'ext',
            'rperiod': 100}

    dir_ = '../SearchForPatt/SG_900_pC/Cand_S'
    pars = findParameters(dir_)
    ksix = zeros((0,998))
    ksiA = zeros((0,998))
    for f0 in pars['f0']:
        pattAd = dir_ + '/patterns_P_0.86_f0_%.2f.npy' %f0
        kx = data2array(pattAd.replace('Cand_S','Cand_X'))
        kA = data2array(pattAd.replace('Cand_X','Cand_S'))
        ksix = concatenate((ksix, kx))
        ksiA = concatenate((ksiA, kA))
    T0 = 20000 # 0.1 ms
    TF = 80000 # 0.1 ms
    sx = 0.2
    sT = 0.2
    x, A = [], []


    ############################################################################################'
    # **Loop**
    tic()
    for iii in range(100):
        # Init
        m = randint(ksix.shape[0])
        other.update({'x':ksix[m], 
                      'A':ksiA[m]})
        eva = main.evaCure(evaCon=conn, evaNoi=noise, evaMod=model, out=out, **other)

        #Run
        for i in range(T0):
            eva.update()
        eva.evaNoi.updateNoise(stdD_x=sx, stdD_T=sT)
        for i in range(TF):
            eva.update()
            
        #Save
        x.append( array(eva.out['x']) )
        A.append( array(eva.out['A']).mean(1) )
        
        del eva
        print iii
    tac()


    ############################################################################################'
    # **Save**
    array2data(array(x), './SG/x_%i.npy'%p1); del x
    array2data(array(A), './SG/A_%i.npy'%p1); del A


else:
    print 'File exists.'
