#!/usr/bin/env python
#-*- coding:Utf-8 -*-
''' Downsampling spatial for Time Courses from 998 nodes to 66 nodes
'''

'#################################################  LIBRAIRIES  ######################################################'
import sys
sys.path.append('../evaCure/.')
from FromToFiles import *
from pylab import zeros, append, linspace, sqrt, exp, sin, convolve
from numpy import load as np_load, save as np_save
from os import path as os_path, mkdir as os_mkdir, listdir as os_listdir
from time import time, sleep
from scipy.stats import ss

'##################################################  FUNCTIONS  #######################################################'
def mkdir(adress):
    ''' Crée un dossier s'il n'existe pas
    '''
    if not os_path.exists(adress):
        os_mkdir(adress)

def isdir(adress, f=None):
    if not os_path.exists(adress):
        print 'Pas de dossier ' + adress
        if f != None: f()

def PearsonCorrelation(TC):
    #TODO change with coravriance, see Plot.py
    ''' Return the Pearson Correlation. The calculus is done for HALF of the matrix and duplicate for symetry.
    TC need to have a shape like [time, nodes]
    Attention for TC which is modified by the function
    '''
    from scipy.stats import ss
    from pylab import sqrt as np_sqrt, dot as np_dot
    TC -= TC.mean(axis=0)[(slice(None,None,None),None)].T
    TCss= np_sqrt(ss(TC, axis=0))
    N   = TC.shape[1]
    corr= zeros((N, N))
    for i in xrange(N):
        corr[i,i:]  = np_dot(TC[:,i:].T, TC[:,i])
        corr[i,i:] /= TCss[i:] * TCss[i]
        corr[i:,i]  = corr[i,i:]
    return corr

def fastCovariance(TC):
    TCm = (TC - TC.mean(axis=0)).T
    N   = TC.shape[1]
    cov= zeros((N, N))
    for i in xrange(N):
        cov[i,i:]  = TCm[i:].dot(TCm[i].T)
        cov[i:,i]  = cov[i,i:]
        #print i,'/',N
    return cov, TCm

def fastPearsonCorrelation(TC):
    N        = TC.shape[1]
    corr,TCm = fastCovariance(TC)
    TCss     = sqrt(ss(TCm, axis=1))
    for i in xrange(N):
        corr[i,i:] /= TCss[i:] * TCss[i]
        corr[i:,i]  = corr[i,i:]
    return where(isfinite(corr), corr, 0.)

def numberOfTimeCourses(adress):
    ''' Retourne le nombre de time courses dans le dossier "adress"
    '''
    nb_TC = 0
    list_tri = os_listdir(adress)
    for i in xrange(len(list_tri)):
        if 'TC' in list_tri[i]:
            nb_TC += 1
    return nb_TC

def HRF(ms = 1.):
    ''' The Heamodynamic response function
    '''
    T     = 10
    tau_s = 0.8
    tau_f = 0.4
    scale = 1000. / ms
    v_time= linspace(0., T, scale * T)
    sqrt_tfts = sqrt(1./ tau_f - 1./ (4.* tau_s ** 2))
    exp_ts    = exp(-0.5 * (v_time / tau_s))
    h         = exp_ts * sin(v_time * sqrt_tfts) / sqrt_tfts
    return h


'##################################################  PARAMETERS  ######################################################'
while len(sys.argv) > 1:
    option = sys.argv[1];                                             del sys.argv[1]
    if option == '-p1':    p1 = float(sys.argv[1].replace(',','.'));  del sys.argv[1]
    elif option == '-p2':  p2 = float(sys.argv[1].replace(',','.'));  del sys.argv[1]
    elif option == '-w':   wait = int(sys.argv[1].replace(',','.'));  del sys.argv[1]
    else:
        print 'Options invalides :',option,'->',sys.argv[0]

sleep(wait)

D998_D66 = array([ \
     0,  0,  1,  1,  1,  1,  1,  1,  2,  2,  2,  3,  3,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,\
     4,  5,  5,  5,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,\
     7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,\
     8,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, 10, 10, 10, 10, 10,\
    10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,\
    12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 14, 14,\
    14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15,\
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 16, 16, 16,\
    16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 19,\
    19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,\
    21, 21, 21, 21, 21, 21, 21, 22, 22, 22, 22, 23, 23, 23, 23, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24,\
    24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 25, 25, 25, 25, 25, 25, 25,\
    25, 25, 25, 25, 25, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 28, 28,\
    28, 28, 28, 28, 28, 28, 28, 28, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 30, 30, 30, 30, 30, 30, 30, 30, 31, 31, 31, 31, 31, 31,\
    31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 32, 32, 32, 32, 32, 32, 32, 33, 33, 33, 33, 33, 33, 33, 34, 34, 34,\
    34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 35, 35, 35, 35, 35, 35, 35, 35, 36, 36, 36, 36, 36, 36,\
    36, 36, 36, 36, 36, 37, 37, 37, 37, 37, 37, 37, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39,\
    39, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41,\
    41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 42, 42,\
    42, 42, 43, 43, 43, 43, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 45, 45, 45, 45, 45, 45, 46, 46,\
    46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 47, 47, 47, 47, 47, 47, 47, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48,\
    48, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,\
    50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51,\
    51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52,\
    52, 52, 52, 52, 52, 52, 52, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 54, 54, 54, 54, 54, 55, 55, 55,\
    55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56,\
    56, 56, 56, 56, 56, 56, 56, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58,\
    58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59,\
    59, 59, 59, 59, 59, 60, 60, 60, 60, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 62, 62, 63,\
    63, 63, 63, 64, 64, 64, 64, 64, 64, 65, 65, 65])
AvgN_D66 = array([ \
    2,6,3,2,22,3,19,27,19,28,16,7,20,28,31,36,13,10,8,22,6,19,4,4,46,12,17,10,10,12,8,23,7,7,\
    23,8,11,8,9,16,12,50,4,4,20,6,19,7,11,13,36,30,29,19,5,19,25,17,27,22,4,22,2,4,6,3])
N1, N2   = 998, 66
step     = 10000
avgOn    = 2000
dir_TC1  = './TC_%i/' %N1;   isdir(dir_TC1)
dir_TC2  = './TC_%i/' %N2;   mkdir(dir_TC2)
dir_BO1  = './BO_%i/' %N1;   mkdir(dir_BO1)
dir_BO2  = './BO_%i/' %N2;   mkdir(dir_BO2)
dir_FC1  = './FC_%i/' %N1;   mkdir(dir_FC1)
dir_FC2  = './FC_%i/' %N2;   mkdir(dir_FC2)
#name     = 'Psi_%.2f_sigmax_%.2f' %(p1,p2)
#name     = 'tauT_%.2f_sigmaT_%.2f' %(p1,p2)
#name     = 'Omega_%.2f_sigmaT_%.2f' %(p1,p2)
name     = 'sigma2x_%.2f_sigma2T_%.2f' %(p1,p2)
dir_sub  = name + '/'

nTC  = numberOfTimeCourses(dir_TC1 + dir_sub)
tmax = nTC * step
TC1  = zeros((tmax, N1))
TC2  = zeros((tmax, N2))
BO1  = zeros((tmax / avgOn, N1))
BO2  = zeros((tmax / avgOn, N2))
Hrv  = HRF()


'#####################################################  MAIN  ###########################################################'

'####################  TC_998'
t0 = time()
for i in range(nTC):
    TC1[i*step : (i+1)*step, :] = np_load(dir_TC1 + dir_sub + 'TC_%i.npy'%i)
print 'TC 998 loaded      : %.2fs'%(time()-t0)

'####################  BO_998'
t0 = time()
for k in range(N1):
    Bk = convolve(TC1[:,k], Hrv)[:tmax]
    for o in range(avgOn):
        BO1[:,k] += Bk[o::avgOn]
BO1 /= avgOn
np_save(dir_BO1 + 'BO_%i_'%N1 + name + '.npy', BO1)
print 'BO 998 generated   : %.2fs'%(time()-t0)

'####################  BO_66, FC_B66, FC_B998'
t0 = time()
for n in range(N1):
    BO2[:, D998_D66[n]] += BO1[:,n]
for t in range(BO1.shape[0]):
    BO2[t,:] /= AvgN_D66
np_save(dir_BO2 + 'BO_%i_'%N2 + name + '.npy', BO2)
FC2 = fastPearsonCorrelation(BO2);                          del BO2
np_save(dir_FC2 + 'FC_B%i_'%N2 + name + '.npy', FC2);   del FC2
FC1 = fastPearsonCorrelation(BO1);                          del BO1
np_save(dir_FC1 + 'FC_B%i_'%N1 + name + '.npy', FC1);   del FC1
print 'BO66, FCB66, FCB998: %.2fs'%(time()-t0)

'####################  TC_66, FC_A66, FC_A998'
t0 = time()
for n in range(N1):
    TC2[:, D998_D66[n]] += TC1[:,n]
for t in range(TC1.shape[0]):
    TC2[t,:] /= AvgN_D66
np_save(dir_TC2 + 'TC_%i_'%N2 + name + '.npy', TC2)
FC2 = fastPearsonCorrelation(TC2);                          del TC2
np_save(dir_FC2 + 'FC_A%i_'%N2 + name + '.npy', FC2);   del FC2
FC1 = fastPearsonCorrelation(TC1);                          del TC1
np_save(dir_FC1 + 'FC_A%i_'%N1 + name + '.npy', FC1);   del FC1
print 'BO66, FCA66, FCA998: %.2fs'%(time()-t0)
