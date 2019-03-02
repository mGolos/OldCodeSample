#!/usr/bin/env python

import warnings
warnings.filterwarnings("ignore",category=Warning)
from pylab import linspace, array, zeros, arange, corrcoef, isfinite, where
from time import time
from numpy import load as np_load
import multiprocessing
import commands
import os
import sys
sys.path.append('../evaCure/.')
from FromToFiles import *



def whatyouwant(paramin):
    global dir_pri, dir_pri2

    dens_moy = arange(0.02,0.981,0.03) # (33)
    #dens_moy = [0.74] # (33)
    for d in range(len(dens_moy)):
        dir_sub = '/Norm_%.2f_DensMoy_%.2f' %(paramin,dens_moy[d])
        if not os.path.exists(dir_pri + dir_sub):
            print 'Pas de sous-dossier ' + dir_sub
            return 0
        if not os.path.exists(dir_pri2 + dir_sub):
            os.mkdir(dir_pri2 + dir_sub)
        if not os.path.exists(dir_pri2 + dir_sub + '/ok'):
            nb_patterns = 0
            err = -1
            nb_de_candidats = len(os.listdir(dir_pri + dir_sub)) / 2

            if nb_de_candidats < 250:
                for i in range(nb_de_candidats):
                    candidat = dir_pri + dir_sub + '/pattern_%.3i.npy' %i
                    #destination = '../FilesIN/pattern_000.npy'
                    ##destination = '../FilesIN/initialization.npy'
                    #cmd = commands.getoutput('cp ' + candidat + " " + destination)

                    if candidat.sum() == 0:
                        break

                    #print i,nb_de_candidats,paramin,dens_moy[d]
                    prog = "../evaCure/Main.py"
                    sys.argv = [prog , "-dens_moy",str(dens_moy[d]), \
                                    '-normW',str(paramin), \
                                    '-fileIN', str(candidat), \
                                '-N','66', \
                                '-m','1', \
                                '-dt','0.01', \
                                '-dur','100', \
                                '-a','0.', \
                                '-b','1.', \
                                '-theta','0.5', \
                            '-nper','3000', \
                            '-err','0.00001', \
                                '-GUI','0', \
                                '-save','0', \
                            '-test','2', \

                                '-netw','0', \
                                '-patt','1', \
                            '-init','3', \
                                '-who','0', \
                            #'-fileIN','../FilesIN/initialization.npy', \
                                '-dens','0.5', \
                            '-noise','0.', \
                                '-conn','5', \
                                '-upd','0', \
                                '-Dx','0.1', \
                                '-posit','0', \
                                '-ap0e0','0', \
                                '-mT','0', \
                                '-p_var','0.', \
                            #'-normW','2.', \
                                '-tau_theta','1', \
                                '-theta_ref','0.5', \
                            #'-dens_moy','0.3', \

                                '-priT','0', \
                                '-resh','1', \
                                '-nbcol','1', \

                                '-thre','0', \
                            '-gamma','30', \
                                '-derr','0.2', \
                                '-tau','1', \
                                '-intFDx','0']
                    execfile(prog)

                    #Network.closers()
                    #coooo = corrcoef(np_load(candidat), Network.states)[0,1]
                    #err = 1. - where(isfinite(coooo), coooo, 0)
                    #print i,nb_de_candidats,paramin,dens_moy[d] , coooo, err
                    err = 1 - corrcoef(np_load(candidat), Network.states)[0,1]
                    if err < 1e-3:
                    #if Network.dpos < 1e-3:
                        destination2 = dir_pri2 + dir_sub + '/pattern_%.3i.npy' %nb_patterns
                        nb_patterns += 1

                        cmd = commands.getoutput('cp ' + candidat + " " + destination2)
            #print 'Norm %.2f, DensMoy %.2f, candidats %i, retenu %i, dcorr %.2e' %(paramin,dens_moy[d],nb_de_candidats,nb_patterns,Network.dpos)
            cmd = commands.getoutput('touch ' + dir_pri2 + dir_sub + '/ok')
            print 'Norm %.2f, DensMoy %.2f, candidats %i, retenu %i, dcorr %.2e' %(paramin,dens_moy[d],nb_de_candidats,nb_patterns,err)
    return 1

if __name__ == '__main__':
    dir_pri  = './Cand_Disc_Dyn'
    dir_pri2 = './Patt_Disc_Dyn'

    if not os.path.exists(dir_pri):
        print 'Pas de dossier ' + dir_pri
        quit()
    if not os.path.exists(dir_pri2):
        os.mkdir(dir_pri2)

    t0 = time()
    pool = multiprocessing.Pool()
    args = arange(0.50,4.501,0.2).tolist() # normW
    #args = [1.3] # normW
    results = pool.map(whatyouwant, args)

    print "Simulation took :",time()-t0,"s"
