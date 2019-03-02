#!/usr/bin/env python
#-*- coding:Utf-8 -*-

############################################################################################'
# **Libraries**
import sys, os.path
sys.path.append( os.path.expanduser('/home/mathieuG') )

from evaCure import main
from Tools.ext import array2data, data2array, mkdir, adressExists
from Tools.display import tic, tac
from Tools.functions import fPearsonCorrelation
from pylab import array


############################################################################################'
# **Parameters**
while len(sys.argv) > 1:
    option = sys.argv[1];                                            del sys.argv[1]
    if   option == '-P':   P = float(sys.argv[1].replace(',','.'));  del sys.argv[1]
    elif option == '-G':   G = float(sys.argv[1].replace(',','.'));  del sys.argv[1]
    elif option == '-sx':  sx = float(sys.argv[1].replace(',','.')); del sys.argv[1]
    elif option == '-sT':  sT = float(sys.argv[1].replace(',','.')); del sys.argv[1]
    elif option == '-dir': dir_= str(sys.argv[1]);                   del sys.argv[1]
    else:
        print 'Options invalides :',option,'->',sys.argv[0]

mdir = '/home/mathieuG/TCS/%s' %dir_
name = 'G_%.2f_sx%.3f_sT%.3f.npy' %(G, sx, sT)  

dir_FC = mdir + '/FC_998'
dir_TC = mdir + '/TC_998'
try:
    mkdir(dir_FC)
    mkdir(dir_TC)
except:
    pass


############################################################################################'
# **Main**
if not adressExists(dir_TC + '/TC_998_' + name):

    conn = {'connAd': "/home/mathieuG/Connectomes/Hagmann/SC_D_998_0.npy",
            'normC': 1.}

    model = {'model': 'HopfieldBasedStatic',
            'threshold': 'local',
            'tauT': 80,
            'P': 1.,
            'G': G}

    noise = {'stdD_x': sx,
             'stdD_T': sT,
             'colors': ['white','white']}

    out = ['A']

    other = {'init': 'rand',
            'dens': 0.5,
            'rperiod': 100,
            'dur':'20m'}

    eva = main.evaCure(evaCon=conn, evaNoi=noise, evaMod=model, out=out, **other)
    eva.updateTillEnd()
    
    # !!!!!!!!!!!!!!!! BINARISATION
    TC = eva.out['A'] > 0.5
    # !!!!!!!!!!!!!!!! 
    array2data(TC[len(TC)/20:], dir_TC + '/TC_998_' + name)
    
if not adressExists(dir_FC + '/FC_998_' + name):
    try:
        TC = array(TC)
    except:
        TC = data2array(dir_TC + '/TC_998_' + name)
        
    tic()
    FC = fPearsonCorrelation(TC)
    tac('h')
    array2data(FC, dir_FC + '/FC_998_' + name)
    
    
