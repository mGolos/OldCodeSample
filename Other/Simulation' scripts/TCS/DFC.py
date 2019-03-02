#!/usr/bin/env python
#-*- coding:Utf-8 -*-

############################################################################################'
# **Libraries**
import sys, os.path
sys.path.append( os.path.expanduser('/home/mathieuG') )

from evaCure import main
from Tools.ext import array2data, data2array, mkdir, adressExists
from Tools.display import tic, tac
from Tools.functions import fPearsonCorrelation, windowedFCs, windowedCorrelations
from pylab import array


############################################################################################'
# **Parameters**
while len(sys.argv) > 1:
    option = sys.argv[1];                                            del sys.argv[1]
    if option == '-G':   G = float(sys.argv[1].replace(',','.'));  del sys.argv[1]
    elif option == '-sx':  sx = float(sys.argv[1].replace(',','.')); del sys.argv[1]
    elif option == '-sT':  sT = float(sys.argv[1].replace(',','.')); del sys.argv[1]
    elif option == '-dir': dir_= str(sys.argv[1]);                   del sys.argv[1]
    else:
        print 'Options invalides :',option,'->',sys.argv[0]

mdir = '/home/mathieuG/TCS/%s' %dir_
name = 'G_%.2f_sx_%.3f_sT_%.3f.npy' %(G, sx, sT)  

dir_DFC = mdir + '/DFC_998/'
dir_TC = mdir + '/TC_998/'
try:
    mkdir(dir_DFC)
except:
    pass


############################################################################################'
# **Main**
if not adressExists(dir_DFC + name):

    TC = data2array(dir_TC + name)
    FCs = windowedFCs(TC, window=6000, jump=100)
    DFC = windowedCorrelations(FCs)
    array2data(DFC, dir_DFC + name)
    
    