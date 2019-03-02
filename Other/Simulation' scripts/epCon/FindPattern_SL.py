#!/usr/bin/env python
#-*- coding:Utf-8 -*-

############################################################################################'
# **Libraries**
import sys, os.path
sys.path.append( os.path.expanduser('/home/mathieuG') )

from evaCure import main
from Tools.ext import Pdir, array2data, mkdir
from Tools.functions import similarity
from pylab import zeros, arange, rand
import os


############################################################################################'
# **Parameters**
while len(sys.argv) > 1:
    option = sys.argv[1];                                             del sys.argv[1]
    if   option == '-p1':  p1 = float(sys.argv[1].replace(',','.'));  del sys.argv[1]
    elif option == '-p2':  p2 = float(sys.argv[1].replace(',','.'));  del sys.argv[1]
    else:
        print 'Options invalides :',option,'->',sys.argv[0]

dir_priA = './SL/Cand_S'
dir_prix = './SL/Cand_X'
mkdir(dir_priA)
mkdir(dir_prix)
dir_sub = '/patterns_P_%.2f_f0_%.2f.npy' %(p1,p2)
#dir_sub = '/patterns_P_%.2f_G_%.2f.npy' %(p1,p2)


############################################################################################'
# **Main**
multi_dens = arange(0.02,0.981,0.03) # (33)
all_pattx = zeros((3300,512)) #998))
all_pattA = zeros((3300,512)) #998))

conn = {'connAd': Pdir('Connectomes/ab-p21-epi_1.txt'),
        'normType': '1'}

noise = {'colors': None}

model = {'model': 'HopfieldBasedStatic',
         'threshold': 'local',
         'tauT': 0,
         'P': p1,
         'G': 900}

out = []

other = {'init': 'rand',
         'dens': p2,
         'rperiod': 100}

for d in range(3300):
    eva = main.evaCure(evaCon=conn, evaNoi=noise, evaMod=model, out=out, **other)
    eva.toEquilibrium()
    all_pattx[d] = eva.evaMod.x.copy()
    all_pattA[d] = eva.evaMod.A.copy()
        
array2data(all_pattx, dir_prix + dir_sub)
array2data(all_pattA, dir_priA + dir_sub)
