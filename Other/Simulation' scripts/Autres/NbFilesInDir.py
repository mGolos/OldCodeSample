#!/usr/bin/env python


import sys
from pylab import *
import os

p1 = arange(0.5,30,1)
p2 = arange(0.1,20,0.7)
dir_pri = './Simu2/Patt_S'
repartition = zeros((len(p1), len(p2)))

for i in arange(len(p1)):
    for j in arange(len(p2)):
        dir_ = dir_pri + '/Norm_%.2f_ThetaRef_%.2f' %(p1[i], p2[j])
        #nb_patt = len(os.listdir(dir_pri + '/Norm_%.2f_DensMoy_%.2f' %(i,j)))
        nb_patt = len(os.listdir(dir_)) - 1
        #nb_patt = len(os.listdir(dir_pri + '/Norm_%.2f_DensMoy_%.2f' %(i,j))) / 2

        repartition[i,j] = nb_patt

print repartition
