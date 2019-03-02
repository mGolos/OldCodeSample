#!/usr/bin/env python
#-*- coding:Utf-8 -*-

############################################################################################'
# **Libraries**
import sys, os.path
sys.path.append( os.path.expanduser('/home/mathieuG') )

from evaCure import main
from Tools.ext import Pdir, array2data, mkdir
from Tools.functions import similarity
from pylab import zeros, arange
import os


############################################################################################'
# **Parameters**
while len(sys.argv) > 1:
    option = sys.argv[1];                                             del sys.argv[1]
    if   option == '-p1':  p1 = float(sys.argv[1].replace(',','.'));  del sys.argv[1]
    elif option == '-p2':  p2 = float(sys.argv[1].replace(',','.'));  del sys.argv[1]
    else:
        print 'Options invalides :',option,'->',sys.argv[0]

dir_priA = './DL/Cand_S'
dir_prix = './DL/Cand_X'
mkdir(dir_priA)
mkdir(dir_prix)
dir_sub = '/patterns_P_%.2f_f0_%.2f.npy' %(p1,p2)


############################################################################################'
# **Main**
multi_dens = arange(0.02,0.981,0.03) # (33)
all_pattx = zeros((33,100,998))
all_pattA = zeros((33,100,998))
  
conn = {'connAd': Pdir('Connectomes/SC_FB_D_998_0.npy'),
        'normType': '1'}

noise = {'colors': None}

model = {'model': 'HopfieldBasedDynamic',
            'threshold': 'local',
            'tauT': 80,
            'P': p1,
            'G': 900}

out = []

other = {'init': 'rand',
            'dens': p2,
            'rperiod': 100}

for d in range(33):
    for i in range(100):
        eva = main.evaCure(evaCon=conn, evaNoi=noise, evaMod=model, out=out, **other)
        eva.toEquilibrium()
        all_pattx[d,i] = eva.evaMod.x.copy()
        all_pattA[d,i] = eva.evaMod.A.copy()
        
array2data(all_pattx, dir_prix + dir_sub)
array2data(all_pattA, dir_priA + dir_sub)