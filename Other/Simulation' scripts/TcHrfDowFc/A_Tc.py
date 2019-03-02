#!/usr/bin/env python
#-*- coding:Utf-8 -*-
''' Lance une simulation de 1000s et génère 100 fichiers TC_i.npy représentant la Time Courses totale moyennée sur
10 pas de temps
'''

'#################################################  LIBRAIRIES  ######################################################'
import sys
sys.path.append('../evaCure/.')
from FromToFiles import *
from os import path as os_path, mkdir as os_mkdir, listdir as os_listdir
from time import time


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


'##################################################  PARAMETERS  ######################################################'
while len(sys.argv) > 1:
    option = sys.argv[1];                                             del sys.argv[1]
    if option == '-p1':    p1 = float(sys.argv[1].replace(',','.'));  del sys.argv[1]
    elif option == '-p2':  p2 = float(sys.argv[1].replace(',','.'));  del sys.argv[1]
    else:
        print 'Options invalides :',option,'->',sys.argv[0]

dir_pri = './TC_998'
#dir_sub = dir_pri + '/tauT_%.2f_sigmaT_%.2f/' %(p1,p2)
#dir_sub = dir_pri + '/Psi_%.2f_sigmax_%.2f/' %(p1,p2)
#dir_sub = dir_pri + '/Omega_%.2f_sigmaT_%.2f/' %(p1,p2)
dir_sub = dir_pri + '/sigma2x_%.2f_sigma2T_%.2f/' %(p1,p2)
isdir(dir_pri)
mkdir(dir_sub)


'#####################################################  MAIN  ###########################################################'
t0 = time()
prog = "../evaCure/Main.py"
sys.argv = [prog ,  '-sigma2x',str(p1), \
                    '-sigma2T',str(p2), \
                    '-TCOUT',str(dir_sub), \

                '-GUI','0', \
                '-netw','1', \
                '-init','0', \
                '-patt','0', \
                '-local','0', \
                '-N','998', \
                '-m','0', \

                '-dt','0.01', \
            '-dur','100010', \
            '-nper','0', \
                '-err','1e-15', \

                '-taux','1.0', \
                '-tauT','5.0', \
                '-tauW','0.0', \
                '-P','1.0', \
                '-G','60.0', \
                #'-sigma2x','0.0', \
                #'-sigma2T','0.0', \
                '-sigma2A','0.0', \
                '-sigma20','0.0', \
                '-dens','0.0', \
                '-theta','0.5', \

                '-who','0', \
                '-fileIN','../FilesIN/initialization.npy', \
                '-direIN','../FilesIN/', \
                '-Cnnctm','../Connectomes/SC_FB_D_998_0.npy', \
                '-priT','0', \
            '-DSTC','1', \
                '-TtS','100000']#, \
                #'-TCOUT','../FilesIN/Time_Courses/',]

execfile(prog)

print 'sigma2x:%.2f ; sigma2T:%.2f ; dur:%.2f'%(p1,p2,time()-t0)
#print 'Omega:%.2f ; sigmaT:%.2f ; dur:%.2f'%(p1,p2,time()-t0)
#print 'tauT:%.2f ; sigmaT:%.2f ; dur:%.2f'%(p1,p2,time()-t0)

