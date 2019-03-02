#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore",category=Warning)
from pylab import linspace, array, zeros, arange, corrcoef, isfinite, where, norm
from time import time
from numpy import load as np_load, save as np_save
import multiprocessing
import commands
import os
import sys
sys.path.append('../evaCure/.')
from FromToFiles import *


def whatyouwant(paramin):
    global dir_pri1, dir_pri2, dir_pri3
    normW2 = arange(1.4,2.601,0.05).tolist()
    dens_moy = arange(0.02,0.981,0.03)

    for d1 in range(len(dens_moy)):
        dir_sub1 = dir_pri1 + '/Norm_%.2f_DensMoy_%.2f' %(paramin,dens_moy[d1])
        nbpatterns1 = len(os.listdir(dir_sub1))

        for p1 in range(nbpatterns1):
            d_patty1 = dir_sub1 + '/pattern_%.3i.npy' %p1
            patty1 = np_load(d_patty1)

            for w2 in range(len(normW2)):
                dir_sub2 = dir_pri2 + '/Norm_%.2f' %(normW2[w2])
                nbpatterns2 = len(os.listdir(dir_sub2)) / 2

                for p2 in range(nbpatterns2):
                    d_patty2 = dir_sub2 + '/pattern_%.3i.npy' %p2
                    patty2 = np_load(d_patty2)
                    d_patty0 = dir_sub2 + '/states_%.3i.jpeg' %p2

                    if  (corrcoef(patty1, patty2)[0,1] > 0.99) \
                    and (abs(norm(patty1) - norm(patty2)) / norm(patty2) < 0.0005):
                        direction = dir_pri3 + '/C_W%.2fD%.2fP%.3i_D_W%.2fP%.3i' %(paramin,dens_moy[d1],p1,normW2[w2],p2)
                        cmd = commands.getoutput('cp ' + d_patty2 + " " + direction + '.npy')
                        cmd = commands.getoutput('cp ' + d_patty0 + " " + direction + '.jpeg')
                        print direction, abs(norm(patty1) - norm(patty2)) / norm(patty2)
    return 1

if __name__ == '__main__':
    dir_pri1 = './Patt_Gamma_30'
    dir_pri2 = './Premieres/VarNorm_Gamma_30_ss_vartheta'
    dir_pri3 = './Similaire_gamma30_DiscEtCont'

    if not os.path.exists(dir_pri1):
        print 'Pas de dossier ' + dir_pri1
        quit()
    if not os.path.exists(dir_pri2):
        print 'Pas de dossier ' + dir_pri2
        quit()
    if not os.path.exists(dir_pri3):
        os.mkdir(dir_pri3)

    t0 = time()
    pool = multiprocessing.Pool()
    args = arange(0.50,4.501,0.2).tolist()
    results = pool.map(whatyouwant, args)

    print "Simulation took :",time()-t0,"s"
