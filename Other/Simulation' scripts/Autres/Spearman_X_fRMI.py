#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore",category=Warning)
from pylab import *
from time import time
from scipy.misc import toimage
from numpy import load as np_load, save as np_save
import multiprocessing
import commands
import os
import sys
sys.path.append('../evaCure/.')
from FromToFiles import *
from scipy.stats import spearmanr


def whatyouwant(paramin):
    global dir_pri1, dir_pri2, dir_pri3, patient, NormW

    ofi1 = open(dir_pri2 + '/TC%ia.dat' %patient,'r')
    ofi2 = open(dir_pri2 + '/TC%ib.dat' %patient,'r')
    A = array([line.split() for line in ofi1]).astype(float)
    B = array([line.split() for line in ofi2]).astype(float)
    C = concatenate([A,B], axis=1)
    ofi1.close()
    ofi2.close()

    dir_sub1   = dir_pri1 + '/Norm_%.2f_DensMoy_%.2f' %(NormW, paramin)
    nbpatterns = len(os.listdir(dir_sub1)) - 1
    nbtimeshap = size(C, axis=1)
    nruter     = zeros((nbpatterns, nbtimeshap))

    for p1 in range(nbpatterns):
        d_patty = dir_sub1 + '/pattern_%.3i.npy' %p1
        patty = np_load(d_patty)

        for p2 in range(nbtimeshap):
            nruter[p1,p2] = spearmanr(patty,C[:,p2])[0]

    return nruter

if __name__ == '__main__':
    patient  = 5
    NormW    = 4
    dir_pri1 = './AT_tT1t_X_Patt'
    dir_pri2 = '../FilesIN/Time_Courses/only_positives'
    dir_pri3 = './fRMI_only_positives_Spearman'

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
    args = arange(0.02,0.981,0.03).tolist()
    resultslist = pool.map(whatyouwant, args)
    resultsarra = resultslist[0]
    for ll in range(1, len(resultslist)):
        resultsarra = append(resultsarra, resultslist[ll], axis=0)

    np_save(dir_pri3 + '/correlations_W%.2f_p%i.npy' %(NormW, patient), resultsarra)

    print "Simulation took :",time()-t0,"s"
    print "min= %.2f ; max= %.2f" %(resultsarra.min(), resultsarra.max())
    imshow(abs(resultsarra), cmap=get_cmap('Blues'), interpolation='nearest')
    show()


