#!/usr/bin/env python
# -*- coding: utf-8 -*-


############################################################################################'
# **Libraries**
import sys, os.path
sys.path.append( os.path.expanduser('/home/mathieuG') )

from evaCure import main
from evaCure import parameters
from Tools.ext import tic, tac, loadPatterns, Pdir, array2data, adressExists
from Tools.display import mapMatrices
from pylab import randint, array, rand
import sys


############################################################################################'
# **Parameters**
while len(sys.argv) > 1:
    option = sys.argv[1];                                           del sys.argv[1]
    if option == '-p1':    p1 = int(sys.argv[1].replace(',','.'));  del sys.argv[1]
    else:
        print 'Options invalides :',option,'->',sys.argv[0]

if not adressExists('./DL/x_%i.npy'%p1):

    conn = {'connAd': Pdir('Connectomes/SC_FB_D_998_0.npy'),
            'normType': '1'}
    
    noise = {'stdD_x': 0.0,
             'stdD_T': 0.0,
             'colors': ['white','white']}

    model = {'model': 'HopfieldBasedDynamic',
            'threshold': 'local',
            'tauT': 80,
            'P': 10.,
            'G': 900.}

    out = ['x', 'A']

    other = {'init': 'rand',
            'dens': 0.5,
            'rperiod': 100,}


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
        other.update({'dens':rand()})
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
    array2data(array(x), './DL/x_%i.npy'%p1); del x
    array2data(array(A), './DL/A_%i.npy'%p1); del A
    
    
else:
    print 'File exists.'
