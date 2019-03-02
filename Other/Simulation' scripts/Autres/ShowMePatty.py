#!/usr/bin/env python
#-*- coding:Utf-8 -*-
''' Affiche les patterns en fonction de gamma, normW et la densité
@author: gmoaltos@hotmail.com
'''

from pylab import *
from numpy import load as LOAD
import sys

nb_pattern = len(sys.argv)-1
patterns = []
for i in range(nb_pattern):
    patterns.append( LOAD(str(sys.argv[i+1])) )
    #print patterns[i-1]


abreviations = ['ENT','PARH','TP','FP','FUS','TT','LOCC','SP','IT','IP',\
                'SMAR','BTST','MT','ST','PSTC','PREC','CMF','POPE','PTRI','RMF',\
                'PORB','LOF','CAC','RAC','SF','MOF','LING','PCAL','CUN','PARC',\
                'ISTC','PCUN','PC']
abre_ind = []
for i in range(len(abreviations)): abre_ind.append('r'+ abreviations[i])
abreviations.reverse()
for i in range(len(abreviations)): abre_ind.append('l'+ abreviations[i])


fig = figure()
for i in range(nb_pattern):
    a = fig.add_subplot(1, nb_pattern, i+1)
    a.imshow(array([patterns[i]]).T, interpolation="nearest", cmap=get_cmap('gist_yarg'), vmin=0, vmax=1)
    a.set_yticks(range(66))
    a.set_yticklabels(abre_ind, fontsize=9)
    a.xaxis.set_visible(False)

show()
