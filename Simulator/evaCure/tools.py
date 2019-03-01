#!/usr/bin/env python
#-*- coding:Utf-8 -*-

''' @author: mathieu.golos@gmail.com
'''

from numpy import array, save as np_save, load as np_load, where, linspace, array
from numpy.random import randn
from pylab import disp
from time import time
import os


def loadPatterns(adress, m=-1):
    '''m define the maximum number of patterns loaded.'''
    try:
        S = data2array(adress + '/patterns.npy')
        
    except:
        S = []
        if os.path.exists(adress):
            mx = numberOfPatterns(adress)
            if m == -1 or m > mx:
                m = mx
                
            for i in xrange(m):
                S.append( data2array(adress + '/pattern_%.3i.npy' %i) )
            S = array(S)
        
        else:
            print "%s doesn't exists" %adress
            
    return S


def numberOfPatterns(adress):
    ''' Retourne le nombre de patterns dans le dossier "adress"
    '''
    try:
        nb_patt = len(data2array(adress))
    except:
        nb_patt = 0
        list_tri = os.listdir(adress + '/')
        for i in xrange(len(list_tri)):
            if 'patt' in list_tri[i]: nb_patt = nb_patt + 1
    return nb_patt


def data2array(adress, mmap_mode=None, dic=False, other=None, finite=False, struc=None, typ=float):
    ''' Return a data file ("*.dat/txt" or "*.npy") into a 2D numpy array.'''
    try:
        if adress.endswith('.npy'):
            fl = np_load(adress, mmap_mode=mmap_mode)
            if dic:
                return dict(fl.all())
            else:
                return fl
        
        elif adress.endswith('.txt')\
        or   adress.endswith('.asc')\
        or   adress.endswith('.dat'):
            ofi = open(adress, 'r')
            arr = array([line.split() for line in ofi])
            ofi.close()
            if finite:
                arr = where(arr == 'NA', 'nan', arr)
            try: 
                return arr.astype(typ)
            except:
                return arr
        
        elif adress.endswith('.mat'):
            if struc == 'np':
                def groupLoop(f):
                    nruter = {}
                    for _ in f.keys():
                        k = str(_)
                        if 'Dataset' in str(type(f[k])):
                            nruter[k] = f[k].value
                        elif 'Group' in str(type(f[k])):
                            nruter[k] = groupLoop(f[k])
                    return nruter
                
                import h5py
                f = h5py.File(adress)
                return groupLoop(f)
                
            elif struc == 'hdf5':
                import h5py
                f = h5py.File(adress)
                return f
                
            else:
                from scipy.io import loadmat
                if other is None:
                    return loadmat(adress)
                else:
                    return loadmat(adress)[other]
        
    except ImportError:
        print '%s does not exists or the format is wrong (*.npy, *.txt, *.dat), \
               the file should be written with tabulation between values.' %adress
        return 0


def array2data(arr, adress, other='C', delimiter='\t'):
    '''Save an array (max 2D for dat/txt) as "*.dat/txt" or "*.npy".'''
    if adress.endswith('.npy'):
        np_save(adress, arr)
        
    elif adress.endswith('.txt')\
    or   adress.endswith('.asc')\
    or   adress.endswith('.dat'):
        ofi = open(adress, 'w')
        try:
            for i in xrange(arr.shape[0]):
                for j in xrange(arr.shape[1]):
                    ofi.write(str(arr[i,j]) + delimiter)
                ofi.write('\n')
        except:
            for i in xrange(arr.shape[0]):
                ofi.write(str(arr[i]) + delimiter)
        ofi.close()
        
    elif adress.endswith('.mat'):
        from scipy.io import savemat
        if type(arr) == dict:
            savemat(adress, arr)
        else:
            savemat(adress, {other:arr})


def estimAndPercent(i, Tmax, avg=100, t0=None):
    global time_t0
    
    if t0 != None: tUsed = t0
    else:          tUsed = time_t0
    
    if i == 0:
        return tic()
        
    t1 = tac(ret=True, t0=tUsed)
        
    if not((i+1)%(Tmax/10)):
        print '(%.2i%%), time elapsed: %.2fs' %(100*(i+1)/Tmax, t1)
       
    elif i == avg:
        Te = Tmax * t1 / avg
        disp('Time estimated: %ih%.2imin%.2is' %((Te/3600)%24, (Te/60)%60, Te%60))
    

def tic(ret=False):
    global time_t0
    time_t0 = time()
    
    if ret == True:
        return time_t0


def tac(ret=False, t0=None, st=""):
    global time_t0
    
    if t0 != None: tUsed = t0
    else:          tUsed = time_t0
        
    if ret == True:
        return time() - tUsed
    
    elif ret == 'h':
        Te = time() - tUsed
        print 'Time elapsed: %ih%.2imin%.2is' %((Te/3600)%24, (Te/60)%60, Te%60) + st
        
    else:
        print "Time elapsed: %.2fs" %(time() - tUsed) + st


def whiteNoise(dim1=1000, dim2=None, stdD=1):
    '''Return a 1D or 2D array of white noise with a variance var.
    '''
    if dim2 is None:
        data = randn(dim1)
    else:
        data = randn(dim1, dim2)
    return data * stdD


def noiseClass(ntype='white', dim1=1000, stdD=1., **kwa):
    if ntype is None:
        return noNoise()
        
    elif ntype == 'white':
        a = whiteNoiseClass(dim1=dim1, stdD=stdD, **kwa)
        return a
        
    elif ntype == 'pink':
        return pinkNoiseClass(dim2=dim1, stdD=stdD, **kwa)
        
    elif ntype == 'brown':
        return brownNoiseClass(dim2=dim1, stdD=stdD, **kwa)


class noNoise():
    def update(self):
        return 0.
        
        
class whiteNoiseClass():
    def __init__(self, dim1=1000, stdD=1., **kwa):
        self.dim1 = dim1
        self.kwa = kwa
        self.stdD = stdD
        
    def update(self):
        return whiteNoise(self.dim1, stdD=self.stdD, **self.kwa)
     
    
def pinkNoiseClass(**kwa):
    '''Return brownNoiseClass with AonF=0.'''
    kwa.update({'AonF':0})
    return brownNoiseClass(**kwa)


class brownNoiseClass():
    def __init__(self, dim2=None, stdD=1, AonF=0):
        sTree = int(5 + AonF)
        self.magic_number = (sTree+1) **.5
        self.updateTable = equalTree(sTree)
        self.lenT = len(self.updateTable) 
        self.stdD = stdD #/ magic_number 
        self.dim2 = dim2
        self.tableInd = 0
        
        if dim2:
            self.whiteValues = randn(sTree, dim2)
            self.update = self.update2D
        else:
            self.whiteValues = randn(sTree)
            self.update = self.update1D
    
    def update1D(self):
        self.whiteValues[self.updateTable[(self.tableInd+1) % self.lenT]] = randn()
        self.tableInd += 1
        return (self.whiteValues.sum() + randn()) * (self.stdD / self.magic_number)
    
    def update2D(self):
        self.whiteValues[self.updateTable[(self.tableInd+1) % self.lenT]] = randn(self.dim2)
        self.tableInd += 1
        return (self.whiteValues.sum(0) + randn(self.dim2)) * (self.stdD / self.magic_number)


def equalTree(n=5):
    '''Return a mirror list composed like:
    0,1,0,2,0,1,0,3,0,1,0,2,0,1 (for n=4)
    '''
    tree = []
    for i in range(n):
        tree.append(i)
        tree.extend(tree[:-1])
    del tree[-1]
    return tree
    

