#!/usr/bin/env python
#-*- coding:Utf-8 -*-
''' Calcule les corrélations entre les patterns des différents jeux de paramètres avec les time courses
    BOLD fRMI data issues des 5 patients (pas de 2 secondes) et enregistre le tableau dans un fichier
    correlations_*.npy dans le dossier choisi.
    @author: mathieu.golos@gmail.com
'''

import warnings
warnings.filterwarnings("ignore",category=Warning)
from pylab import *
from time import time
from scipy.misc import toimage
from numpy import load as np_load, save as np_save
from os import path as os_path, mkdir as os_mkdir, listdir as os_listdir
from sys import path as sys_path, argv as sys_argv
import multiprocessing
import commands
sys_path.append('../evaCure/.')
from FromToFiles import *
from Wcorrelation import *

def mkdir(adress):
    ''' Crée un dossier s'il n'existe pas
    '''
    if not os_path.exists(adress):
        os_mkdir(adress)

def isdir(adress, f=None):
    if not os_path.exists(adress):
        print 'Pas de dossier ' + adress
        if f != None: f()

def numberOfPatterns(adress):
    ''' Renvoie le nombre de patterns dans le dossier "adress"
    '''
    nb_patt = 0
    list_tri = os_listdir(adress + '/')
    for i in xrange(len(list_tri)):
        if 'patt' in list_tri[i]: nb_patt = nb_patt + 1
    return nb_patt

def whatyouwant(paramin):
    global dir_pri1, dir_pri2, dir_pri3, patient, nbrNoeu, noise
    Omega = arange(1,30.001,1)

    t1 = time()
    ofi1 = open(dir_pri2 + '/TC_%i_%ia.dat' %(nbrNoeu, patient),'r')
    ofi2 = open(dir_pri2 + '/TC_%i_%ib.dat' %(nbrNoeu, patient),'r')
    A = array([line.split() for line in ofi1]).astype(float)
    B = array([line.split() for line in ofi2]).astype(float)
    TC = array(concatenate([A,B], axis=1)).T
    ofi1.close()
    ofi2.close()

    nruter = zeros(len(Omega))
    for d1 in xrange(len(Omega)):
        dir_sub1   = dir_pri1 + '/Norm_%.2f_Omega_%.2f' %(paramin, Omega[d1])
        nbpatterns = numberOfPatterns(dir_sub1)
        if nbpatterns == 0:
            nruter[d1] = 0
        else:
            nbtimeshap = size(TC, axis=1)
            ksi = zeros((nbpatterns, nbrNoeu))
            for p1 in xrange(nbpatterns):
                d_patty = dir_sub1 + '/pattern_%.3i.npy' %p1
                ksi[p1] = np_load(d_patty)

            W = Wcorrelation(TC, ksi)
            RTC = dot(W, ksi)
            CCC = [corrcoef(TC[t], RTC[t])[0,1] for t in xrange(len(TC))]
            nruter[d1] = sum(CCC) / len(CCC)
    print 'Noise:%i, Patient:%i, Norm:%.2f, Time:%.2f' %(noise, patient, paramin, time()-t1)

    return nruter

if __name__ == '__main__':
    if len(sys_argv) == 4:
        patient = int(sys_argv[1])
        noise   = int(sys_argv[2])
        nbrNoeu = int(sys_argv[3])

        dir_pri1 = './66fbdenRS/At%.2i_S' %noise
        dir_pri2 = '../FilesIN/Time_Courses'
        dir_pri3 = './fRMI_66RS/At%.2i' %noise

        isdir(dir_pri1, quit)
        isdir(dir_pri2, quit)
        mkdir(dir_pri3)

        t0 = time()
        pool = multiprocessing.Pool()
        args = arange(1,30.001,1).tolist()
        results = pool.map(whatyouwant, args)
        results = array(results)
        results = where(isnan(results), 0, results)

        np_save(dir_pri3 + '/correlations_S_p%i.npy' %patient, results)

        ind_a, ind_b = unravel_index(results.argmax(), results.shape)

        print "Patient %i, noise %.2i, min %.2f, max %.2f (%.2i;%.2i), time %.2fs" \
            %(patient, noise, results.min(), results.max(), args[ind_a], args[ind_b], time()-t0)
        #imshow(results.T, cmap=get_cmap('RdYlGn'), interpolation='nearest', origin="lower", vmin=0,vmax=1)
        #xlabel(r'$||W||$')
        #ylabel(r'$\Omega$')
        #colorbar()
        #show()
    else:
        print '(Arguments needed) 1:Patient, 2: Noise, 3: N'