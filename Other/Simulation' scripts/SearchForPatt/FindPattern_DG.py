#!/usr/bin/env python
#-*- coding:Utf-8 -*-

############################################################################################'
# **Libraries**
import sys, os.path
sys.path.append( os.path.expanduser('/home/mathieuG') )

from evaCure import main
from Tools.ext import Pdir, array2data, mkdir, data2array, adressExists
from Tools.functions import similarity
from pylab import zeros, arange, rand
import os


############################################################################################'
# **Parameters**
while len(sys.argv) > 1:
    option = sys.argv[1];                                              del sys.argv[1]
    if   option == '-p1':   p1 = float(sys.argv[1].replace(',','.'));  del sys.argv[1]
    elif option == '-p2':   p2 = float(sys.argv[1].replace(',','.'));  del sys.argv[1]
    elif option == '-dir':  dir_= str(sys.argv[1]);                    del sys.argv[1]
    elif option == '-con':  dco = str(sys.argv[1]);                    del sys.argv[1]
    else:
        print 'Options invalides :',option,'->',sys.argv[0]

mdir = '/home/mathieuG/SearchForPatt/%s' %dir_
dir_priA = mdir + '/Cand_S'
dir_prix = mdir + '/Cand_X'
#dir_sub = '/patterns_P_%.2f_f0_%.2f.npy' %(p1,p2)
dir_sub = '/patterns_P_%.2f_G_%.2f.npy' %(p1,p2)  
mmax = 3300

############################################################################################'
# **Main**
if not adressExists(dir_prix + dir_sub):
    mkdir(dir_priA)
    mkdir(dir_prix)

    multi_dens = arange(0.02,0.981,0.03) # (33)
    N = data2array(dco).shape[0]
    all_pattx = zeros((mmax,N))
    all_pattA = zeros((mmax,N))

    conn = {'connAd': dco,
            'normType': '1'}

    noise = {'colors': None}

    model = {'model': 'HopfieldBasedDynamic',
             'threshold': 'global',
             'tauT': 80,
             'P': p1,
             'G': p2} #900}

    out = []

    other = {'init': 'rand',
             'dens': 0.5, #rand(), #p2,
            'rperiod': 100}

    for d in range(mmax):
        eva = main.evaCure(evaCon=conn, evaNoi=noise, evaMod=model, out=out, **other)
        eva.toEquilibrium()
        all_pattx[d] = eva.evaMod.x.copy()
        all_pattA[d] = eva.evaMod.A.copy()
            
    array2data(all_pattx, dir_prix + dir_sub)
    array2data(all_pattA, dir_priA + dir_sub)
