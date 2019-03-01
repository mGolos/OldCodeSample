#!/usr/bin/env python
#-*- coding:Utf-8 -*-

from scipy.stats import ranksums, mannwhitneyu
from scipy.fftpack import dct
from scipy.stats.stats import pearsonr
from scipy.signal import hilbert
from numpy import tanh, arctanh, sqrt, isfinite, where, zeros, array, var as variance, diag, empty, fill_diagonal
from numpy.random import randn
from pylab import corrcoef, ones, unravel_index, convolve, concatenate, newaxis, c_
from pylab import permutation, argmax, mean, arange, rand, find, floor, triu, dot, tril
from Tools.lists import sortBy
from Tools.stats import equalTree, triSupPearsonCor
from Tools.matrices import triSup, triInf
from itertools import combinations
from matplotlib.mlab import normpdf
from sklearn.decomposition import NMF


gaussian = normpdf

dct2D = lambda x: dct(dct(x.T, norm='ortho').T, norm='ortho')

envelope = hilbert


def ss(x, axis=None):
    return sum(x**2, axis=axis)


def significance(l1, l2, fun=ranksums):
    pval = fun(l1, l2)[1]
    if pval > 0.05:
        return 'n.s'
    elif pval > 0.01 and pval <= 0.05:
        return '*'
    elif pval > 0.001 and pval <= 0.01:
        return '**'
    else:
        return '***'


def timeTest(fs, args, N):
    #from timeit import timeit
    from time import time
    for f, i in zip(fs, range(len(fs))):
        #timeit('"-".join(f(args[i]))', number=10000)
        t0 = time()
        for n in range(N):
            f(*args[i])
        print 'Time elapsed: %es (%s)' %(time()-t0, f.func_name)
    

def sigmoid(x, H=.5, Q=1, G=60, P=1, T=0.5):
    '''Sigmoidal thresholding function.'''
    return H * (Q + tanh(G * (P * x - T)))


def sigmoidM1(A, H=.5, Q=1, G=60, P=1, T=0.5):
    '''Inverse of sigmoidal thresholding function.'''
    #TODO Add epsilon error, Q and H
    return (arctanh(1.99 * (A - 0.5)) / G + T) / P


def similarity_Manhattan(arr, mat=None):
    '''Return the similarity between one vector and a vector/matrix using 
    the Manhattan distance as : 1 / (1 + d) with adjusting parameters for 
    the minimum/maximum similarity of 0/1.'''
    if mat == None:
        L = len(arr)
        nruter = ones((L,L))
        for p in range(L):
            nruter[p,p+1:] = nruter[p+1:,p] = similarity_Manhattan(arr[p], arr[p+1:])
        return array(nruter)
    else:
        return 1. / (1. + abs(arr-mat).sum(axis=-1))


def similarity_Euclidean(arr, mat=None):
    '''Return the similarity between one vector and a vector/matrix using 
    the Euclidean distance as : 1 / (1 + d) with adjusting parameters for 
    the minimum/maximum similarity of 0/1.'''
    if mat == None:
        L = len(arr)
        nruter = ones((L,L))
        for p in range(L):
            nruter[p,p+1:] = nruter[p+1:,p] = similarity_Euclidean(arr[p], arr[p+1:])
        return array(nruter)
    else:
        return 1. / (1. + ((arr-mat)**2).sum(axis=-1)**.5)


def distance_Hamming(arr, mat=None):
    if mat == None:
        L = len(arr)
        nruter = ones((L,L))
        for p in range(L):
            nruter[p,p+1:] = nruter[p+1:,p] = distance_Hamming(arr[p], arr[p+1:])
        return array(nruter)
    else:
        return abs(arr-mat).sum(axis=-1)


def distance_Euclidean(arr, mat=None):
    if mat == None:
        L = len(arr)
        nruter = ones((L,L))
        for p in range(L):
            nruter[p,p+1:] = nruter[p+1:,p] = distance_Euclidean(arr[p], arr[p+1:])
        return array(nruter)
    else:
        return ((arr-mat)**2).sum(axis=-1)**.5


similarity = similarity_Euclidean


def otherDistance(arr, mat):
    '''TODO Define.'''
    return mat.dot(arr) / (arr**2).sum(axis=-1)**(0.5) / (mat**2).sum(axis=-1)**(0.5)


def similAndPear(A, B=None, mmm='max'):
    '''A have to be 1 dimensional.'''
    sP = fPearsonCorrelation(A.T, B)
    sE = similarity_Euclidean(A, B)
    
    if A.ndim == 2:
        if mmm == 'max':
            return concatenate((sP[:,newaxis], sE[:,newaxis]), 1).max(1)
        elif mmm == 'mean':
            return concatenate((sP[:,newaxis], sE[:,newaxis]), 1).mean(1)
        elif mmm == 'min':
            return concatenate((sP[:,newaxis], sE[:,newaxis]), 1).min(1)
        elif mmm == None:
            return sP[:,newaxis], sE[:,newaxis]
    else:
        if mmm == 'max':
            return max(sP, sE)
        elif mmm == 'mean':
            return mean((sP, sE))
        elif mmm == 'min':
            return min(sP, sE)
        elif mmm == None:
            return sP, sE
    

def fPearsonCorrelation(TC1, finite=False):
    '''Return the temporal correlation...
    finite=True change NaN values to 0.
    '''
    corr,TCm1 = fCovariance(TC1)
    Tavg1     = sqrt(ss(TCm1, axis=1))
    for i in xrange(TC1.shape[1]):
        corr[i,i:] /= Tavg1[i:] * Tavg1[i]
        corr[i:,i]  = corr[i,i:]
    if finite:
        return where(isfinite(corr), corr, 0.)
    else:
        return corr
    

def corrTCs(TC1, TC2):
    corr = empty((len(TC1), len(TC2)))
    for i in range(len(TC1)):
        corr[i] = fCovariance(TC1[i], TC2)
    return corr
       

def fPearsonCorrelationTriangle(TC, finite=False, ind=None):
    '''Return the temporal correlation...
    finite=True change NaN values to 0.
    '''
    if ind == None:
        ind = triSup(empty((TC.shape[1], TC.shape[1])), ind=True)
    
    corr, TCm = fCovarianceTriangle(TC)
    Tavg      = sqrt(ss(TCm, axis=1))
    for i in xrange(TC.shape[1]):
        corr[i,i:] /= Tavg[i:] * Tavg[i]
        
    corr = corr.take(ind)
    
    if finite:
        return where(isfinite(corr), corr, 0.)
    else:
        return corr
        

def fCovarianceTriangle(TC):
    '''
    '''
    N   = TC.shape[1]
    TCm = (TC - TC.mean(axis=0)).T
    cov = empty((N, N))
    for i in xrange(N):
        cov[i,i:]  = TCm[i:].dot(TCm[i].T)
    return cov, TCm
       

def fCovariance(TC1, TC2=None):
    '''
    '''
    if TC2 == None:
        if TC1.ndim == 2:
            N    = TC1.shape[1]
            TCm1 = (TC1 - TC1.mean(axis=0)).T
            cov  = zeros((N, N))
            for i in xrange(N):
                cov[i,i:]  = TCm1[i:].dot(TCm1[i].T)
                cov[i:,i]  = cov[i,i:]
            return cov, TCm1
        
    else:
        if TC1.ndim == 2:
            TCm1 = (TC1 - TC1.mean(axis=0)).T
            TCm2 = (TC2 - TC2.mean(axis=0)).T
            N1, N2 = TC1.shape[1], TC2.shape[1]
            cov  = zeros((N1, N2))
            
            if N1 == N2:
                for i in xrange(N1):
                    #cov[i,i:] = TCm1[i:].dot(TCm2[i].T)
                    cov[i,:] = TCm1[:].dot(TCm2[i].T)
                    #cov[i:,i] = cov[i,i:]
            else:
                for i in xrange(N1):
                    cov[i,:] = TCm1[[i]].dot(TCm2[:].T)
                    
            return cov, TCm1, TCm2
        
        elif TC1.ndim > 2:
            TCm1 = (TC1 - TC1.mean(axis=0)).T
            TCm2 = (TC2 - TC2.mean(axis=0)).T
            N1, N2 = TC1.shape[1], TC2.shape[1]
            cov  = ones((N1, N2))
            
            if N1 == N2:
                for i in xrange(N1):
                    cov[i,i:] = TCm1[i:].dot(TCm2[i].T)
                    cov[i:,i] = cov[i,i:]
            else:
                for i in xrange(N1):
                    cov[i,:] = TCm1[[i]].dot(TCm2[:].T)
                    
            return cov, TCm1, TCm2

        elif TC1.ndim == 1:
            if TC2.ndim == 2:
                V   = TC2 - TC2.mean(axis=1)[(slice(None,None,None),None)]
                V0  = TC1 - TC1.mean()
                Vss = sqrt(ss(V, axis=1))
                cov = V.dot(V0)
                cov/= Vss * sqrt(ss(V0))
                return cov
            else:
                return pearsonr(TC1,TC2)[0]


def fPCA(TC):
    from scipy.linalg import eigh
    CO = fCovariance(TC)[0]
    V,PC = eigh(CO)
    iPC, V = sortBy(V, inverse=True)
    return V, PC.T[iPC]


def matricesCorrelation(M1, M2, corrFunc=pearsonr, diag=True, finite= False, posit=False, avg=True, **kwa):
    ''' Return the correlation between the two matrices, defined by the mean across
    correlations for each columns. If posit=True, the correlation is done for positive matrices.
    '''
    if M1.shape != M2.shape:
        print 'Functional Connectomes shape are not the same'
        return 0
    
    N = M1.shape[0]
    corr = zeros(N)
    
    if diag:
        if corrFunc == corrcoef:
            if posit:
                for c in xrange(N):
                    corr[c] = corrcoef(where(M1[c]<0, 0, M1[c]), where(M2[c]<0, 0, M2[c]))[0,1]
            else:
                for c in xrange(N):
                    corr[c] = corrcoef(M1[c], M2[c])[0,1]
        elif corrFunc == pearsonr:
            if posit:
                for c in xrange(N):
                    corr[c] = pearsonr(where(M1[c]<0, 0, M1[c]), where(M2[c]<0, 0, M2[c]))[0]
            else:
                for c in xrange(N):
                    corr[c] = pearsonr(M1[c], M2[c])[0]
        else:
            if posit:
                for c in xrange(N):
                    corr[c] = corrFunc(where(M1[c]<0, 0, M1[c]), where(M2[c]<0, 0, M2[c]), **kwa)
            else:
                for c in xrange(N):
                    corr[c] = corrFunc(M1[c], M2[c], **kwa)
    else:
        ind = arange(N)
        if corrFunc == corrcoef:
            if posit:
                for c in xrange(N):
                    corr[c] = corrcoef(where(M1[c][ind!=c]<0, 0, M1[c][ind!=c]), 
                                       where(M2[c][ind!=c]<0, 0, M2[c][ind!=c]))[0,1]
            else:
                for c in xrange(N):
                    corr[c] = corrcoef(M1[c][ind!=c], M2[c][ind!=c])[0,1]
        elif corrFunc == pearsonr:
            if posit:
                for c in xrange(N):
                    corr[c] = pearsonr(where(M1[c][ind!=c]<0, 0, M1[c][ind!=c]),
                                       where(M2[c][ind!=c]<0, 0, M2[c][ind!=c]))[0]
            else:
                for c in xrange(N):
                    corr[c] = pearsonr(M1[c][ind!=c], M2[c][ind!=c])[0]
        else:
            if posit:
                for c in xrange(N):
                    corr[c] = corrFunc(where(M1[c][ind!=c]<0, 0, M1[c][ind!=c]),
                                       where(M2[c][ind!=c]<0, 0, M2[c][ind!=c]), **kwa)
            else:
                for c in xrange(N):
                    corr[c] = corrFunc(M1[c][ind!=c], M2[c][ind!=c], **kwa)
            
    if finite:
        corr = where(isfinite(corr), corr, 0.)

    if avg:
        return corr.mean()
    else:
        return corr


def windowedFCs(TC, window=0.1, jump=1, nodes=None, sym=False, **kwa):
    ''' Returns the functionnal connectivity matrices based on a certain
    temporal window specified in percentage(<0) or number of step(>0).
    The time have to be the first dimension of TC array.
    '''
    T, N = TC.shape
    
    if window < 1:
        lWin = int(T * window)
    else:
        lWin = window
        
    if sym:
        mask = triSup(zeros((N,N)), ind=True)
        FCs = zeros(((T - lWin) / jump, (N*(N-1))/2))
        for i in xrange(0, FCs.shape[0]):
            FCs[i] = fPearsonCorrelation(TC[i*jump : i*jump + lWin], **kwa).take(mask)
        
    else:
        if nodes == None:
            FCs = zeros(((T - lWin) / jump, N, N))
            for i in xrange(0, FCs.shape[0]):
                FCs[i] = fPearsonCorrelation(TC[i*jump : i*jump + lWin], **kwa)
        
        else:
            try:
                n = len(nodes)
                FCs = zeros(((T - lWin) / jump, n, N))
                for i in xrange(0, FCs.shape[0]):
                    FCs[i] = fPearsonCorrelation(TC[i*jump : i*jump + lWin, nodes], 
                                                TC[i*jump : i*jump + lWin], **kwa)
                
            except:
                FCs = zeros(((T - lWin) / jump, N))
                for i in xrange(0, FCs.shape[0]):
                    FCs[i] = fPearsonCorrelation(TC[i*jump : i*jump + lWin, nodes], 
                                                TC[i*jump : i*jump + lWin].T, **kwa)        
        
    return FCs


def windowedFCsSym(TC, window=0.1, jump=1,  **kwa):
    ''' Returns the functionnal connectivity matrices based on a certain
    temporal window specified in percentage(<0) or number of step(>0).
    The time have to be the first dimension of TC array.
    30% faster than windowedFCs function.
    '''
    T, N = TC.shape
    ind = triSup(empty((N,N)), ind=True)
    
    if window < 1:
        lWin = int(T * window)
    else:
        lWin = window
        
    FCs = empty(((T - lWin) / jump, (N*(N-1))/2))
    for i in xrange(0, FCs.shape[0]):
        FCs[i] = fPearsonCorrelationTriangle(TC[i*jump : i*jump + lWin], ind=ind, **kwa)
        
    return FCs

    
def windowedCorrelations(FCs1, FCs2=None, mask=None, **kwa):
    ''' Returns a correlation matrix between time dependant functional connectivities.
    '''
    L1,n,N = FCs1.shape
    
    if FCs2 != None:
        L2 = len(FCs2)
        comp = True
    else:
        L2 = L1
        comp = False
    
    if mask == None and n == N: # or n == N:
        #if n != N:
            #mask = xrange(n * N)
            
        nruter = zeros((L1,L2))
        if not comp:
            for i in xrange(L1):
                for j in xrange(i,L1):
                    nruter[i,j] = triSupPearsonCor(FCs1[i,:,:], FCs1[j,:,:], mask=mask, **kwa)
            for i in xrange(L1):
                nruter[i:,i] = nruter[i,i:]
        
        else:
            for i in xrange(L1):
                for j in xrange(L2):
                    nruter[i,j] = triSupPearsonCor(FCs1[i,:,:], FCs2[j,:,:], mask=mask, **kwa)
                    
    elif n != N:
        nruter = zeros((L1,L1))
        for i in xrange(L1):
            for j in xrange(i,L1):
                nruter[i,j] = triSupPearsonCor(FCs1[i,:,:], FCs1[j,:,:], mask=mask, **kwa)
        for i in xrange(L1):
            nruter[i:,i] = nruter[i,i:]
        
    
    elif len(mask) == N:
        nruter = zeros((N,L1,L1))
        for n in xrange(N):
            nruter[n,:,:] = fPearsonCorrelation(FCs1[:,n,:].T, **kwa)
        
    else:
        nruter = zeros((L1,L1))
        for i in xrange(L1):
            for j in xrange(i,L1):
                nruter[i,j] = triSupPearsonCor(FCs1[i,:,:], FCs1[j,:,:], mask=mask, **kwa)
        for i in xrange(L1):
            nruter[i:,i] = nruter[i,i:]
    
    return nruter


def whiteNoise(dim1=1000, dim2=None, stdD=1):
    '''Return a 1D or 2D array of white noise with a variance var.
    '''
    if dim2 == None:
        data = randn(dim1)
    else:
        data = randn(dim1, dim2)
    return data * stdD
    
    
def pinkNoise(dim1=1000, **kwa):
    '''Return brownNoise with AonF=0.'''
    kwa.update({'AonF':0})
    return brownNoise(dim1=dim1, **kwa)
  
       
def brownNoise(dim1=1000, dim2=None, stdD=1, AonF=0):
    sTree = int(5 + AonF)
    magic_number = (sTree+1) **.5
    updateTable = equalTree(sTree)
    lenT = len(updateTable)
    
    if dim2:
        whiteValues = randn(sTree, dim2)
        data = zeros((dim1, dim2))
        for i in xrange(dim1):
            whiteValues[updateTable[(i+1) % lenT]] = randn(dim2)
            data[i] = (whiteValues.sum(0) + randn(dim2)) / magic_number
    else:
        whiteValues = randn(sTree)
        data = zeros(dim1)
        for i in xrange(dim1):
            whiteValues[updateTable[(i+1) % lenT]] = randn()
            data[i] = (whiteValues.sum() + randn()) / magic_number
            
    return data * stdD 


def noiseClass(ntype='white', dim1=1000, stdD=1., **kwa):
    if ntype == None:
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
    
    
def fft(TS, tmax=None, sr=1, nf=1):
    '''Return the absolute of the real part of the fast fourrier transformation 
    from a time serie 'TS' with a sampling rate 'sr' and a fft resolution 'nf'.
    '''
    from numpy import linspace
    from numpy.fft import fft
    
    N = len(TS)
    if not tmax:
        tmax = N
    #frequencies = sr * linspace(0, tmax/2., N/2., endpoint=True) / nf
    frequencies = linspace(0, tmax/2., N/2., endpoint=True) / nf
    amplitudes = abs(fft(TS, n=nf*N)[:N/2] / N) 
    return frequencies, amplitudes    
    
 
def patternChecking(candidats, patterns, out='patterns', test='equal', areIn=False, opt=0.9):
    ''' Return the candidats that are not in patterns or the extended patterns.
    out - "patterns": pattens extended with the different candidats
        - "candidats": the different candidats from patterns
    test - "equal": equality
         - "corrcoef": Pearson correlation.
    '''

    if not len(patterns):
        try:
            candidats[0][0]
            return candidats
        except:
            return array([candidats])
    
    else:
        if   out == 'patterns':  nruter = list(patterns)
        elif out == 'candidats': nruter = []

        if   test == 'equal':    simil = lambda x,y: x in y
        elif test == 'corrcoef': simil = lambda x,y: corrcoef(x,y)[0,1:].max()
        elif test == 'fPearson': simil = lambda x,y: fPearsonCorrelation(x,y).max()
        
        if areIn:
            testFunc = lambda x,y: simil(x,y) >= opt
        else:
            testFunc = lambda x,y: simil(x,y) <= opt
        
        try:
            candidats[0][0]
            for i in range(len(candidats)):
                if testFunc(candidats[i], patterns):
                    nruter.append(candidats[i])
        except:
            if testFunc(candidats, patterns):
                nruter.append(candidats)

        return nruter


def newList(l0):
    '''Return a list of list and not adresses. (so far for 2D).'''
    l1 = []
    try:
        for i in range(len(l0)):
            l1.append(l0[i][:])
        return l1
    
    except :
        return list(l0)
    

def clustering(patts=None, corr=None, setPatts=None, 
               sim_coef=0.9, sim_func=fPearsonCorrelation, miea='max'):
    """From patts with or without correlations / sets, or from correlation w/wout sets, without patterns.
    NB: Similarity can be calculated from correlations or distances."""
    
    if corr == None:
        corr = sim_func(patts.T)    
        corr-= diag(diag(corr))
    
    if setPatts == None:
        try:
            setP = [[k] for k in range(len(patts))]
        except:
            setP = [[k] for k in range(len(corr))]
    else:
        setP = newList(setPatts)
        
    if   miea == 'min':  miea = type(array(int())).min
    elif miea == 'mean': miea = type(array(int())).mean
    elif miea == 'max':  miea = type(array(int())).max

    sim_final = 1
    while len(setP) > 1 and sim_final > sim_coef:
        sim = sim_coef#-1
        for i in xrange(len(setP)):
            for j in xrange(i+1, len(setP)):
                simEF = miea(corr[setP[i]][:,setP[j]])
                #simEF = corr[setP[i]][:,setP[j]].max()
                if simEF > sim:
                    indE = i
                    indF = j
                    sim = simEF
        sim_final = sim
        
        if sim_final > sim_coef:
            try:
                for i in xrange(len(setP[indF])):
                    setP[indE].append(setP[indF].pop(0))
                setP.remove([])
            except:
                pass

    return setP

  
def clustering2(patts=None, corr=None, setPatts=None, 
                sim_coef=0.9, sim_func=fPearsonCorrelation, acc=True, L=150, miea='max', one=False): 
    """From patts with or without correlations / sets, or from correlation w/wout sets, without patterns.
    NB: Similarity can be calculated from correlations or distances."""
    
    if one:
        f_combi = lambda L: list(combinations(xrange(L), 2))[slice(L - 1)]
        acc = False
    else:
        f_combi = lambda L: list(combinations(xrange(L), 2))
    
    if corr == None:
        corr = sim_func(patts.T)    
        corr-= diag(diag(corr))
    
    if setPatts == None:
        try:
            setP = [[k] for k in range(len(patts))]
        except:
            setP = [[k] for k in range(len(corr))]
    else:
        setP = newList(setPatts)
        
    if len(setP) > L and acc:
        kwa = {'patts': patts,
               'corr': corr,
               'sim_coef': sim_coef,
               'sim_func': sim_func,
               'L': L,
               'miea': miea,
               'one': one}
        set1 = clustering2(setPatts=setPatts[:L], **kwa)
        set2 = clustering2(setPatts=setPatts[L:], **kwa)
        return clustering2(setPatts=set1+set2, acc=False, **kwa)
    
    elif len(setP) == 1:
        return setP
    
    if   miea == 'min':  miea = type(array(int())).min
    elif miea == 'mean': miea = type(array(int())).mean
    elif miea == 'max':  miea = type(array(int())).max
    
    sim_final = 1
    while len(setP) > 1 and sim_final > sim_coef:
        
        combi = f_combi(len(setP))
        sims = []
        for i,j in combi:
            sims.append( miea(corr[setP[i]][:,setP[j]]) )
        sim_final = max(sims)
        
        if sim_final > sim_coef:
            i,j = combi[argmax(sims)]
            setP[i] += setP[j]
            del setP[j]
            
    return setP


def preClustering(S, freq=None, sim_coef=0.90, sim_func=fPearsonCorrelation):
    '''Pre clusterize the set of patterns S, sorting it by its nodes mean
    and removing the elements with a similarity higher than 0.99.
    Return the indice list and the predominances.
    '''
    iS = sortBy(S.mean(1), inverse=1)[0]
    SS = S[iS]
    to_test = range(len(SS))
    if freq == None:
        freq = [1 for f in to_test]
    
    for i_p in range(len(SS)):
        if len(to_test) > i_p:
            corrs = sim_func(SS[to_test[i_p]], SS[to_test[i_p+1:]])
            to_delete = (corrs >= sim_coef).nonzero()[0] + 1
            
            for i_d in to_delete[::-1]:
                freq[i_p] += freq[i_p + i_d]
                del freq[i_p + i_d]
                del to_test[i_p + i_d]
                
    return to_test, freq


def indicesFrom2DList(AllClust):
    '''Return a list of indices.'''
    iClust = []
    for ic in range(len(AllClust)):
        iClust.append(range(len(AllClust[ic])))
    return iClust


def winClus2StackHist0(AllClust, AllFreqs, AllPatts, sim_coef=0.9, sim_func=fPearsonCorrelation, rnorm=1.):
    '''Windowed Clustering to create Stacked Histogram: 
    Returns two lists of lenght equal to the number of attractors for all P values,
    using a simple clustering (Pearson correlation). The two lists contains the predominances
    and their associate P value.
    '''
    freqAssoc = []
    PvalAssoc = []
    NP = len(AllClust)
    
    #Create a list or indices that should become a list of empty lists at the end
    iClust = indicesFrom2DList(AllClust)
    
    #Loop while there is patterns to be compared
    while iClust != NP * [[]]:
        
        #Create lists for the global patterns and initialize correlated patterns indices
        gloPatF = []
        gloPatP = []
        gloPatN = []
        ibs = None
        
        #Look for bloodlink between parameter P(x)=P0 and P(x+1)=P1
        for iP0 in range(NP-1):
            
            #Look if there is patterns for P0
            if iClust[iP0] == []:
                continue
            
            #Reach the last P1 val for P0 if ibs != None
            if ibs != None and iP0 < iP1:
                continue
            
            #Load indices of all patterns for P0 or the previous pattern found for P1
            if ibs == None:
                iPatts0 = newList(iClust[iP0])
            else:
                iPatts0 = [iClust[iP0][ibs[1]]]
            
            #Try to find indices of all patterns for a P1
            for iP1 in range(iP0+1, NP):
                iPatts1 = newList(iClust[iP1])
                if iPatts1 != []:
                    #P1 found with patterns
                    break
            
            #If no more patterns after P0
            if iPatts1 == []:
                
                #Update global pattern informations for good similarity previously found
                if ibs != None:
                    gloPatF.append( AllFreqs[iP0][iPatts0][ibs[1]] / rnorm )
                    gloPatP.append( iP0 )
                    gloPatN.append( iPatts0[ibs[1]] )
                    
                break
            
            #There is patterns for P0 and P1
            else:
                
                #Calcul of similarities between sets of P0 and P1
                sims = zeros((len(iPatts0), len(iPatts1)))
                for ic in range(len(iPatts0)):
                    #sims[ic,:] = sim_func(AllPatts[iP0][iPatts0][ic], AllPatts[iP1][iPatts1])
                    sP = fPearsonCorrelation(AllPatts[iP0][iPatts0][ic], AllPatts[iP1][iPatts1])
                    sE = similarity_Euclidean(AllPatts[iP0][iPatts0][ic], AllPatts[iP1][iPatts1])
                    sims[ic,:] = concatenate((sP[:,newaxis],sE[:,newaxis]), 1).max(1)
                #print sims
                    
                #Tests the bigest similarity
                if sims.max() > sim_coef:
                    #Indices of the two attractors with the biggest similarity
                    ibs = unravel_index(sims.argmax(), (len(iPatts0), len(iPatts1)))
                    
                    #Update global pattern informations
                    gloPatF.append( AllFreqs[iP0][iPatts0][ibs[0]] / rnorm )
                    gloPatP.append( iP0 )
                    gloPatN.append( iPatts0[ibs[0]] )
                    
                #Add uncorrelated patterns as global patterns
                elif ibs == None:
                    for i in iClust[iP0][::-1]:
                        freqAssoc.append( [AllFreqs[iP0][i] / rnorm] )
                        PvalAssoc.append( [iP0] )
                        iClust[iP0].remove( i )
                    continue
                    
                #Failed test: change P0 to the next one, even if possible next concordances
                else:
                    ibs = None
                    continue
                    
        #Put the pattern informations into output
        if gloPatF != []:
            freqAssoc.append( gloPatF )
            PvalAssoc.append( gloPatP )
            for i in range(len(gloPatF))[::-1]:
                iClust[gloPatP[i]].remove( gloPatN[i] )
        
        #For patterns of the last P when no more patterns before
        if iClust[:-1] == (NP-1) * [[]]:
            for i in iClust[-1][::-1]:
                freqAssoc.append( [AllFreqs[-1][i] / rnorm] )
                PvalAssoc.append( [NP-1] )
                iClust[-1].remove( i )
                
    return freqAssoc, PvalAssoc


def winClus2StackHist42(AllPatts, AllTends, sim_coef=0.9, sim_func=similAndPear, rnorm=1., miea='max'):
    '''Windowed Clustering to create Stacked Histogram: 
    Returns two lists of lenght equal to the number of attractors for all P values,
    using a simple clustering (Pearson correlation). The two lists contains the predominances
    and their associate P value.
    '''
    freqAssoc = []
    PvalAssoc = []
    NP = len(AllPatts)
    
    #Create a list or indices that should become a list of empty lists at the end
    nBoxPatt = indicesFrom2DList(AllPatts)
    
    #Loop while there is patterns to be compared
    while nBoxPatt.count([]) != NP:
        
        #Create lists for the global patterns and initialize correlated patterns indices
        gloPatF = []  # Frequencies
        gloPatP = []  # P value
        gloPatN = []  # Pattern indice
        target = None # [P value, indices of last set of patts]
            
        #Look for bloodlink between parameter P(x)=P0 and P(x+1)=P1
        for iPval0 in range(NP):
            
            #Look if there is patterns for P0
            if nBoxPatt[iPval0] == []:
                continue
            
            #Load indices of all patterns for P0
            if target == None:
                
                #More than one pattern in nBoxPatt[iPval0]
                if len(nBoxPatt[iPval0]) > 1:
                    
                    #Names of patterns  in the box for P0
                    nPatts0 = nBoxPatt[iPval0]
                    
                    #Clusterize patterns for P0
                    corr = sim_func(AllPatts[iPval0][nPatts0])
                    corr-= diag(diag(corr))
                    iClustPatt = clustering(corr=corr, sim_coef=sim_coef, miea=miea)
                
                    #Look for the biggest cluster
                    iPatts0 = iClustPatt[0]
                    for clust in iClustPatt:
                        if len(clust) > len(iPatts0):
                            iPatts0 = list(clust)
                    #patts0 = array(nBoxPatt[iPval0])[iPatts0].tolist()
                    
                #Only one pattern in nBoxPatt[iPval0]
                else:
                    #patts0 = nBoxPatt[iPval0]
                    iPatts0 = [0]
                                    
                #Updated names of patterns in the box for P0
                nPatts0 = array(nBoxPatt[iPval0])[iPatts0].tolist()
                
                #Update global pattern informations
                #iFreqs = nBoxPatt[iPval0][iPatts0]
                gloPatF.append( AllTends[iPval0][nPatts0].sum() / rnorm )
                gloPatP.append( iPval0 )
                gloPatN.append( nPatts0 )
                                    
                #Define target with P1 and iCorrelated Patts
                target = iPval0, nPatts0
                
            #Load the previous indices found for P1
            else:
                if iPval0 != target[0]:
                    continue
                else:
                    #patts0 = target[1]
                    nPatts0 = target[1]
            
            
            #Try to find indices of all patterns for a P1
            for iPval1 in range(iPval0 + 1, NP):
                
                ##Stop after 3 P higher thant P0
                #if iPval1 > iPval0 + 2:
                    #target = None
                    #break
                
                #Look if there is patterns for P1
                if nBoxPatt[iPval1] == []:
                    continue
                
                #Names of patterns  in the box for P1
                nPatts1 = nBoxPatt[iPval1]
                       
                ##Load indices of all patterns for P1
                #patts1 = list(nBoxPatt[iPval1])
              
                #SET OF SET for clustering renaming indices of patts1
                setOfset = [ range(len(nPatts0)) ]
                set1pls0 = [ [k + len(nPatts0)] for k in range(len(nPatts1)) ]
                setOfset.extend(set1pls0)
               
                #CORRELATION
                corr = sim_func(concatenate((AllPatts[iPval0][nPatts0], AllPatts[iPval1][nPatts1]))) #[patts1])))
                corr-= diag(diag(corr))
               
                #CLUSTERING
                iClustPatt = clustering(corr=corr, sim_coef=sim_coef, setPatts=setOfset, miea=miea)
            
                      
                #Indices of patts0 in iClustPatt
                #whereiP0 = [patts0[0] in j for j in iClustPatt]
                whereiP0 = [0 in j for j in iClustPatt]
                jiP0Clust = whereiP0.index(True)
                iiP01 = iClustPatt[jiP0Clust]
                iiP0 = range(len(nPatts0))
                    
                ##No previous correlations
                #if target == None:
                                    
                    ##Define target with P1 and iCorrelated Patts
                    #target = iPval0, iPatts0
                    
                    ##Update global pattern informations
                    #iFreqs = [ nBoxPatt[iPval0].index(ip) for ip in patts0 ]
                    #gloPatF.append( AllTends[iPval0][iFreqs].sum() / rnorm )
                    #gloPatP.append( target[0] )
                    #gloPatN.append( target[1] )
                    
                #print "iPval0, iPval1:", iPval0, iPval1
                #If correlation between iPatts0 and set1
                #if patts0 != iiP01:
                if iiP01 != iiP0:
                    
                    #Remove elements of patts0
                    #nPatts1 = iiP01
                    #for iP0 in patts0:
                        #nPatts1.remove(iP0)
                    #nPatts1 = (array(nPatts1) - len(AllPatts[iPval0])).tolist()
                    for iP0 in iiP0:
                        iiP01.remove(iP0)
                    iiP01 = (array(iiP01) - len(iiP0)).tolist()
                    #nPatts1 = iiP01 #array(nBoxPatt[iPval1])[iiP01].tolist()
                    nPatts1 = array(nBoxPatt[iPval1])[iiP01].tolist()
                    
                    #Update global pattern informations
                    #iFreqs = [ nBoxPatt[iPval1].index(ip) for ip in nPatts1 ]
                    gloPatF.append( AllTends[iPval1][nPatts1].sum() / rnorm )
                    gloPatP.append( iPval1 )
                    gloPatN.append( nPatts1 )
                    
                    #Define target with P1 and iCorrelated Patts
                    target = iPval1, nPatts1
                    
                    #Go to an other value of P0
                    break
                               
            #If no P1 available
            if nBoxPatt[iPval1] == []:
                
                #Update global pattern informations
                gloPatF.append( AllTends[iPval0][nPatts0].sum() / rnorm )
                gloPatP.append( iPval0 )
                gloPatN.append( nPatts0 )
                                    
                #Define target with P1 and iCorrelated Patts
                target = iPval0, nPatts0
                
            
            if target == None:
                break
               
                    
        #Put the pattern informations into output
        if gloPatF != []:
            freqAssoc.append( gloPatF )
            PvalAssoc.append( gloPatP )
            for i in range(len(gloPatF))[::-1]:
                for j in range(len(gloPatN[i]))[::-1]:
                    nBoxPatt[gloPatP[i]].remove( gloPatN[i][j] )
                        
        print "len(freqAssoc):", len(freqAssoc), " (%i set remaining)"%(NP-nBoxPatt.count([]))
    return freqAssoc, PvalAssoc


def winClus2StackHist(AllPatts, AllTends, 
                      sim_coef=0.9, sim_func=similAndPear, rnorm=1., miea='max', one=False):
    '''Windowed Clustering to create Stacked Histogram: 
    Returns two lists of lenght equal to the number of attractors for all P values,
    using a simple clustering (Pearson correlation). The two lists contains the predominances
    and their associate P value.
    '''
    freqAssoc = []
    PvalAssoc = []
    IvalAssoc = []
    NP = len(AllPatts)
    
    #Create a list or indices that should become a list of empty lists at the end
    nBoxPatt = indicesFrom2DList(AllPatts)
    
    #Loop while there is patterns to be compared
    while nBoxPatt.count([]) != NP:
        
        #Create lists for the global patterns and initialize correlated patterns indices
        gloPatF = []  # Frequencies
        gloPatP = []  # P value
        gloPatN = []  # Pattern indice
        target = None # [P value, indices of last set of patts]
            
        #Look for bloodlink between parameter P(x)=P0 and P(x+1)=P1
        for iPval0 in range(NP):
            
            #Look if there is patterns for P0
            if nBoxPatt[iPval0] == []:
                continue
            
            #Load indices of all patterns for P0
            if target == None:
                
                #More than one pattern in nBoxPatt[iPval0]
                if len(nBoxPatt[iPval0]) > 1:
                    
                    #Names of patterns  in the box for P0
                    nPatts0 = nBoxPatt[iPval0]
                    
                    #Clusterize patterns for P0
                    corr = sim_func(AllPatts[iPval0][nPatts0])
                    corr-= diag(diag(corr))
                    iClustPatt = clustering2(corr=corr, sim_coef=sim_coef, miea=miea, one=one)
                
                    #Look for the biggest cluster
                    iPatts0 = iClustPatt[0]
                    for clust in iClustPatt:
                        if len(clust) > len(iPatts0):
                            iPatts0 = list(clust)
                    
                #Only one pattern in nBoxPatt[iPval0]
                else:
                    iPatts0 = [0]
                                    
                #Updated names of patterns in the box for P0
                nPatts0 = array(nBoxPatt[iPval0])[iPatts0].tolist()
                
                #Update global pattern informations
                gloPatF.append( AllTends[iPval0][nPatts0].sum() / rnorm )
                gloPatP.append( iPval0 )
                gloPatN.append( nPatts0 )
                                    
                #Define target with P1 and iCorrelated Patts
                target = iPval0, nPatts0
                
            #Load the previous indices found for P1
            else:
                if iPval0 != target[0]:
                    continue
                else:
                    #patts0 = target[1]
                    nPatts0 = target[1]
            
            
            #Try to find indices of all patterns for a P1
            for iPval1 in range(iPval0 + 1, NP):
                
                #Look if there is patterns for P1
                if nBoxPatt[iPval1] == []:
                    continue
                
                #Names of patterns in the box for P1
                nPatts1 = nBoxPatt[iPval1]
              
              
                #SET OF SET for clustering renaming indices of nPatts1
                setOfset = [ range(len(nPatts0)) ]
                set1pls0 = [ [k + len(nPatts0)] for k in range(len(nPatts1)) ]
                setOfset.extend(set1pls0)
               
                #CORRELATION
                corr = sim_func(concatenate((AllPatts[iPval0][nPatts0], AllPatts[iPval1][nPatts1])))
                corr-= diag(diag(corr))
               
                #CLUSTERING
                iClustPatt = clustering2(corr=corr, sim_coef=sim_coef, setPatts=setOfset, miea=miea, one=one)
            
            
                #Indices of nPatts0 in iClustPatt
                whereiP0 = [0 in j for j in iClustPatt]
                jiP0Clust = whereiP0.index(True)
                iiP01 = iClustPatt[jiP0Clust]
                iiP0 = range(len(nPatts0))
                                        
                #If correlation between nPatts0 and nPatts1
                if iiP01 != iiP0:
                    
                    #Remove elements of iiP01
                    for iP0 in iiP0:
                        iiP01.remove(iP0)
                    iiP01 = (array(iiP01) - len(iiP0)).tolist()
                    nPatts1 = array(nBoxPatt[iPval1])[iiP01].tolist()
                    
                    #Update global pattern informations
                    gloPatF.append( AllTends[iPval1][nPatts1].sum() / rnorm )
                    gloPatP.append( iPval1 )
                    gloPatN.append( nPatts1 )
                    
                    #Define target with P1 and iCorrelated Patts
                    target = iPval1, nPatts1
                    
                    #Go to an other value of P0
                    break
                               
            #If no P1 available
            if nBoxPatt[iPval1] == []:
                
                #Update global pattern informations
                gloPatF.append( AllTends[iPval0][nPatts0].sum() / rnorm )
                gloPatP.append( iPval0 )
                gloPatN.append( nPatts0 )
                                    
                #Define target with P1 and iCorrelated Patts
                target = iPval0, nPatts0
               
            if target == None:
                break
              
              
        #Put the pattern informations into output
        if gloPatF != []:
            freqAssoc.append( gloPatF )
            PvalAssoc.append( gloPatP )
            IvalAssoc.append( gloPatN )
            for i in range(len(gloPatF))[::-1]:
                for j in range(len(gloPatN[i]))[::-1]:
                    try:nBoxPatt[gloPatP[i]].remove( gloPatN[i][j] )
                    except:pass
                        
        print "len(freqAssoc):", len(freqAssoc), " (%i set remaining)"%(NP-nBoxPatt.count([]))
    return freqAssoc, PvalAssoc, IvalAssoc


def findMaxDim(l):
    L0, L1 = len(l), 0
    for i in range(L0):
        test = len(l[i])
        if test > L1:
            L1 = test
    return L0, L1


def sortByIrregularDimension(l, inverse=False):
    L0, L1 = findMaxDim(l)
    with0 = zeros((L0,L1))
    
    for i in range(L0):
        with0[i, :len(l[i])] = l[i]
        
    return sortBy(with0.sum(1), inverse=inverse)


def mvAvg(data, window=10, same=False):
    '''Return an average moved over a window with a the same shape
    or a downsampled one.
    '''
    if same:
        window = ones(int(window)) / float(window)
        if data.ndim == 2:
            nruter = zeros(data.shape)
            for i in range(data.shape[1]):
                nruter[:,i] = convolve(data[:,i], window, 'same')
            return nruter
        else:
            return convolve(data, window, 'same')

    else:
        s = list(data.shape)
        s[0] /= int(window)
        nruter = zeros(s)
        for i in range(window):
            nruter += data[i::window]
            
        return nruter / window


def switch(l, i1, i2):
    '''Switch two position of a list.'''
    tmp = l[i1]
    l[i1] = l[i2]
    l[i2] = tmp
    

def HRe_from_LDown(L):
    '''Return the array of indices to reorder a matrix of high dimension H
    using the order of a matrix of lower dimension L.
    '''
    r = array(L)
    s = 0
    for k in xrange(max(L)+1):
        n = (L==k).sum()
        r[L==k] = s + arange(n)
        s += n
    nruter = arange(len(L))
    nruter[r] = arange(len(L))
    return nruter


def AALtoHemi(N, sym=True, H=None):
    '''Return the array of indices to reorder from AAL order 
    to an order starting from the Left hemisphere and finishing
    from the Right one.
    - By default the right hemisphere is inversely reorder to have
    a symmetry from the center. For a matrix the antidiagonal
    corresponds then to the same region but from different hemispheres.
    - TODO: If the dimension is higher than AAL nomenclature, it needs 
    the downsampling law.
    '''
    if H == None:
        H = arange(N)
    nruter = arange(N)
    
    if (H == arange(N)).all():
        if sym:
            inv = where(H%2==0, H/2, H.max() - H/2)
        else:
            inv = where(H%2==0, H/2, H/2 + H.max()/2)
        nruter[inv] = H
        return nruter
            
    else:
        lu, ld = [], []
        su, sd = 0, 0
        avg = array([(H == i).sum() for i in range(max(H)+1)])
        
        for i in range(len(avg)):
            if i%2 == 0:
                lu.append(su + arange(avg[i]))
            else:
                ld.append(sd + arange(avg[i]))
            su += avg[i]
            sd += avg[i]
        if sym:
            nruter = concatenate(lu).tolist() + concatenate(ld).tolist()[::-1]
        else:
            nruter = concatenate(lu).tolist() + concatenate(ld).tolist()
        return array(nruter)


def invOrder(o1):
    '''Return the opposite law of reordering.'''
    o2 = arange(len(o1))
    o2[o1] = arange(len(o1))
    return o2


def nullDirectedNetwork(mat, ntry=None):
    """The adjacency matrix of a randomized network with the same set of in- and out-degrees
    as the original one.
    - mat: the adjacency matrix of an directed network  
    - ntry: (optional) the number of rewiring steps. If none is given ntry=4*(# of edges in the network)
    """
    
    nrew = 0
    null = array(mat)
    i_null, j_null = unravel_index(find(null), mat.shape) # indices of non-null elements
    i_null = array(i_null)
    j_null = array(j_null)
    Ne = len(i_null)
    if ntry == None:
        ntry = 4 * Ne

    for i in range(ntry):
        e1 = floor(Ne * rand())
        e2 = floor(Ne * rand())
        v1 = i_null[e1]
        v2 = j_null[e1]
        v3 = i_null[e2]
        v4 = j_null[e2]
            
        if (v1!=v3) & (v1!=v4) & (v2!=v4) & (v2!=v3):
            if (null[v1,v4] == 0) & (null[v3,v2] == 0):
                null[v1,v4] = null[v1,v2]
                null[v3,v2] = null[v3,v4]
                null[v1,v2] = 0
                null[v3,v4] = 0
                nrew = nrew + 1
                i_null[e1] = v1
                j_null[e1] = v4
                i_null[e2] = v3
                j_null[e2] = v2
    return null


def nullUndirectedNetwork(mat, ntry=None):
    """The adjacency matrix of a randomized network with the same set of in- and out-degrees
    as the original one.
    - mat:  the adjacency matrix of an undirected network  
    - ntry: (optional) the number of rewiring steps. If none is given ntry=4*(# of edges in the network)
    """

    nrew = 0
    null = array(mat)
    i_null, j_null = unravel_index(find(null), mat.shape) # indices of non-null elements
    aux = find(i_null > j_null) #select only values of one triangle
    i_null = i_null[aux]
    j_null = j_null[aux]
    Ne = len(i_null)
    if ntry == None:
        ntry = 4 * Ne

    for i in range(ntry):
        e1 = floor(Ne * rand())
        e2 = floor(Ne * rand())
        v1 = i_null[e1]
        v2 = j_null[e1]
        v3 = i_null[e2]
        v4 = j_null[e2]
            
        if (v1!=v3) & (v1!=v4) & (v2!=v4) & (v2!=v3):
            if rand() > 0.5:
                if not ((null[v1,v3] == 0) & (null[v2,v4] == 0)):
                    v5 = v3
                    v3 = v4
                    v4 = v5
                    del v5
                    
                if (null[v1,v3] == 0) & (null[v2,v4] == 0):
                    null[v1,v3] = null[v1,v2]
                    null[v2,v4] = null[v4,v3]
                    null[v3,v1] = null[v2,v1]
                    null[v4,v2] = null[v3,v4]  
                    null[v1,v2] = 0
                    null[v4,v3] = 0
                    null[v2,v1] = 0
                    null[v3,v4] = 0     
                    nrew = nrew + 1
                    i_null[e1] = v1
                    j_null[e1] = v3
                    i_null[e2] = v2
                    j_null[e2] = v4
    return null


def nullUndirectedHemiNetwork(A):
    """Return undirected null-model matrix using two hemispheres as conditions.
    The intrahemispheric links are swapped using 'nullUndirectedNetwork'
    where the interhemispheric links using 'nullDirectedNetwork'.
    """
    
    N = len(A)
    nruter = zeros((N,N))
    nruter[:N/2,:N/2] = nullUndirectedNetwork(A[:N/2,:N/2])
    nruter[N/2:,N/2:] = nullUndirectedNetwork(A[N/2:,N/2:])
    nruter[N/2:,:N/2] = nullDirectedNetwork(A[N/2:,:N/2])
    nruter[:N/2,N/2:] = nruter[N/2:,:N/2].T

    return nruter


def switchSymetricalValues(A):
    """Return a matrix with the same structure as the input (symetrical matrix)
    but with shuffled values.
    """
    
    nruter = zeros(A.shape)
    indNZ = array(A.nonzero()).T
    indNZ = indNZ[indNZ[:,0] > indNZ[:,1]]
    swiNZ = permutation(indNZ)
    for i,j,k,l in c_[indNZ, swiNZ]:
        nruter[i,j] = nruter[j,i] = A[k,l]
        
    return nruter


def randomMatrixPermutations(mat):
    """Return a matrix with elements of a triangle permuted without saving the degrees.
    """
    
    s = mat.shape
    nruter = array(mat)
    indices = find( triu( ones(s), k=1))    
    aleaInd = permutation(indices)
    nruter[unravel_index(indices, s)] = mat[unravel_index(aleaInd, s)]
    for i in range(s[0]):
        nruter[i:,i] = nruter[i,i:]
        
    return nruter


def random2HemiPermutations(mat):
    """Return a matrix with elements of a triangle permuted without saving the degrees.
    """
    s = mat.shape
    N  = s[0]
    N2 = N/2
    L  = N*(N-2)/8
    nruter = array(mat)
    
    Sind = triSup(arange(N2*N2).reshape((N2,N2)) + arange(0, N2*N2, N2).reshape((N2,1)))
    Mind = (arange(N2*N2).reshape((N2,N2)) + arange(N2, N2*N2+1, N2).reshape((N2,1))).reshape(N2*N2)
    Iind = triSup(ones((N,N)), ind=1)[-L:]
    Salea = permutation(Sind)
    Malea = permutation(Mind)
    Ialea = permutation(Iind)
    
    nruter[unravel_index(Sind, s)] = mat[unravel_index(Salea, s)]
    nruter[unravel_index(Mind, s)] = mat[unravel_index(Malea, s)]
    nruter[unravel_index(Iind, s)] = mat[unravel_index(Ialea, s)]
    for i in range(s[0]):
        nruter[i:,i] = nruter[i,i:]
        
    return nruter
    
   
def binned_srand_internal(sa, bedg1, bedg2):
    sa2 = sa - diag(diag(sa))
    nb1 = len(bedg1)
    nb2 = len(bedg2)
    k_out_a = sa2.T.sum(0)
    k_in_a = sa2.sum(0)

    k_1 = k_out_a
    k_2 = k_in_a

    i1, j1 = unravel_index( find(sa2), sa2.shape)
    n_1_2 = zeros((nb1, nb2))

    for i in range(len(i1)):
        kc1 = k_1[i1[i]]
        kc2 = k_2[j1[i]]
        if (kc1 * kc2) > 0:
            b1 = min(find((bedg1 - kc1) > 0)) - 1
            b2 = min(find((bedg2 - kc2) > 0)) - 1
            n_1_2[b1,b2] += 1
            
    return n_1_2

   
def sym_correlation_profile(s1, Nstat=3, edge_bins=None, disp=False):
    ''' 
    INPUT:
    - s1: the adjacency matrix of an undirected network  
    - srand: nullUndirectedNetwork(s1)
    - Nstat: (optional) the number of randomized networks in the ensemble. Default: 3
    - edge_bins: (otional) the array to bin degrees. Default: [1,3,10,30,100...]
    
    OUTPUT:
    - n_12: number of edges connecting different bins to each other
    - nr_12: same averaged over Nstat realizations of a randomized network
    - nsr_12: sqaure of nr_12 averaged over Nstat realizations of a randomized network
    - R_12: correlation profile ratio: R_12=n_12./nr_12;
    - Z_12: correlation profile Z-score: Z_12=(n_12-nr_12)./sqrt(nsr_12-nr_12.^2);'''

    srand = 1 * (0 < abs( s1 - diag(diag(s1))))
    srand = 1 * (0 < (srand + srand.T))

    k2 = sum(srand)
    k2_max = k2.max()

    if edge_bins == None:
        edge_bins = [1,3]
        m1 = 1
        while edge_bins[m1] <= k2_max:
            edge_bins.append( 10 * edge_bins[m1-1])
            m1 += 1
        
    bedg1 = list(edge_bins)

    if k2_max > bedg1[-1]:
        bedg1.append(k2_max)

    n_1_2_sym_orig = binned_srand_internal(srand, bedg1, bedg1)

    print 'randomized network #',1
    srand = nullUndirectedNetwork(srand)
    n_1_2 = binned_srand_internal(srand,bedg1,bedg1)

    aver_n_1_2_sym = n_1_2
    aver_sq_n_1_2_sym = n_1_2 ** 2
    
    for k in range(1, Nstat):
        print 'randomized network #',k+1
        srand = nullUndirectedNetwork(srand)
        n_1_2 = binned_srand_internal(srand,bedg1,bedg1)
        aver_n_1_2_sym = aver_n_1_2_sym + n_1_2
        aver_sq_n_1_2_sym = aver_sq_n_1_2_sym + n_1_2 ** 2
        
    aver_n_1_2_sym = aver_n_1_2_sym / Nstat
    aver_sq_n_1_2_sym = aver_sq_n_1_2_sym / Nstat
    err_n_1_2_sym = sqrt(aver_sq_n_1_2_sym - aver_n_1_2_sym ** 2)

    sym_ratio_1_2_sym = n_1_2_sym_orig / (aver_n_1_2_sym + 0.0001 * (aver_n_1_2_sym == 0))
    dev_n_1_2_sym_orig = (n_1_2_sym_orig - aver_n_1_2_sym) / (err_n_1_2_sym + 0.0001 * (aver_n_1_2_sym == 0))

    sym_ratio_1_2_sym = sym_ratio_1_2_sym[:-1,:-1]
    dev_n_1_2_sym_orig = dev_n_1_2_sym_orig[:-1,:-1]


    R_12 = sym_ratio_1_2_sym
    Z_12 = dev_n_1_2_sym_orig
    n_12 = n_1_2_sym_orig[:-1,:-1]
    nr_12 = aver_n_1_2_sym[:-1,:-1]
    nsr_12 = aver_sq_n_1_2_sym[:-1,:-1]

    sym_ratio_1_2_sym = concatenate((sym_ratio_1_2_sym, sym_ratio_1_2_sym[-1,:][newaxis]), axis=0)
    sym_ratio_1_2_sym = concatenate((sym_ratio_1_2_sym, sym_ratio_1_2_sym[:,-1][:,newaxis]), axis=1)
    dev_n_1_2_sym_orig = concatenate((dev_n_1_2_sym_orig, dev_n_1_2_sym_orig[-1,:][newaxis]), axis=0)
    dev_n_1_2_sym_orig = concatenate((dev_n_1_2_sym_orig, dev_n_1_2_sym_orig[:,-1][:,newaxis]), axis=1)

    if disp:
        from pylab import figure, meshgrid, pcolor, title, xlabel, ylabel, xscale, yscale, colorbar
        x, y = meshgrid(bedg1, bedg1)

        figure(); title('R(K_1,K_2) for network w/ %s nodes, %s links'%(srand.shape[1], srand.sum()/2))
        pcolor(x, y, sym_ratio_1_2_sym); colorbar()
        xlabel('K_1'); xscale('log')
        ylabel('K_2'); yscale('log')

        figure(); title('Z(K_1,K_2) for network w/ %s nodes, %s links'%(srand.shape[1], srand.sum()/2))
        pcolor(x, y, dev_n_1_2_sym_orig); colorbar()
        xlabel('K_1'); xscale('log')
        ylabel('K_2'); yscale('log')

    return  R_12, Z_12, n_12, nr_12, nsr_12


def lm(l1, obj):
    l2 = list(l1)
    l2.remove(obj)
    return l2
