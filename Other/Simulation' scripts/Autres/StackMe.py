#!/usr/bin/env python
from pylab import *
from numpy import load, save

p1 = arange(0.02,0.981,0.03) #dens
#p2 = arange(0.001,0.3011,0.01) # Dx
#p2 = arange(0.001,0.6011,0.02) # Dx
p2 = arange(0.001,1.2011,0.04) # Dx
#p1 = arange(0.05,0.95,0.05) # dens (18)
#p2 = arange(0.001,0.301,0.02) # Dx (15)
n= 33
results = zeros((len(p1), len(p2), n))

for i in range(len(p1)):
    for j in range(len(p2)):
        results[i,j] = load('result_part_p1_%.5f_p2_%.5f.npy' %(p1[i],p2[j]))

save('results.npy', results)