#!/usr/bin/env python
#-*- coding:Utf-8 -*-

from scipy.stats import ranksums
from scipy.stats.stats import pearsonr
from numpy import sqrt, isfinite, where, zeros, array, diag, empty, sum
from numpy.random import randn
from pylab import corrcoef, ones, unravel_index, convolve, concatenate, newaxis
from pylab import r_, c_, permutation, argmax, mean, arange, find, triu, dot
from itertools import combinations

from Tools.matrices import triSup
from Tools.lists import equalTree, sortBy


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
    
    
def FWER(pval, s=0.1, prin=False):
    '''Family-Wise Error Rate with Bonferroni correction.
    '''
    N = len(pval)
    nruter = array(pval) * N
    for n in range(1,N+1)[::-1]:
        if ((nruter / n) < s).sum() >= n:
            if prin:
                print(n)
            return nruter / n
    if prin:
        print(n)
    return nruter


def slidingAverage(x, N):
    tmp = convolve(x, ones((N,))/N, mode='valid', )
    return r_[tmp[0]*ones(N/2), tmp, tmp[-1]*ones(N/2-1+N%2)]

    
def significativity(pvals,sec=[0.1,0.05,0.01],ind=["n.s","*","**","***"]):
    '''Return a numpy array of strings completed by n.s, *, ** or ***
    depending on the significativity or pvals.
    '''
    tmp = array(pvals)
    return where(tmp > sec[0], ind[0], where(tmp > sec[1], ind[1], where(tmp > sec[2], ind[2], ind[3])))


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
    Tavg      = sqrt(sum(TCm**2, axis=1))
    for i in range(TC.shape[1]):
        corr[i,i:] /= Tavg[i:] * Tavg[i]
        
    corr = corr.take(ind)
    
    if finite:
        return where(isfinite(corr), corr, 0.)
    else:
        return corr
        
   
def fPearsonCorrelation(TC1, TC2=None, finite=False):
    '''Return the temporal correlation...
    finite=True change NaN values to 0.
    '''
    if TC1.ndim == 2:
        if TC2 != None:
            corr,TCm1,TCm2 = fCovariance(TC1, TC2)
            Tavg1 = sqrt(sum(TCm1**2, axis=1))
            Tavg2 = sqrt(sum(TCm2**2, axis=1))
            if len(Tavg1) == len(Tavg2):
                for i in range(TC1.shape[1]):
                    corr[i,i:] /= Tavg1[i:] * Tavg2[i]
                    corr[i:,i]  = corr[i,i:]
            else:
                for i in range(TC1.shape[1]):
                    corr[i,:] /= Tavg1[i] * Tavg2[:]
                    
            if finite:
                return where(isfinite(corr), corr, 0.)
            else:
                return corr
            
        else:
            corr,TCm1 = fCovariance(TC1, TC2)
            Tavg1     = sqrt(sum(TCm1**2, axis=1))
            for i in range(TC1.shape[1]):
                corr[i,i:] /= Tavg1[i:] * Tavg1[i]
                corr[i:,i]  = corr[i,i:]
            if finite:
                return where(isfinite(corr), corr, 0.)
            else:
                return corr
    else:
        return fCovariance(TC1, TC2)
        

def fCovarianceTriangle(TC):
    '''
    '''
    N   = TC.shape[1]
    TCm = (TC - TC.mean(axis=0)).T
    cov = empty((N, N))
    for i in range(N):
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
            for i in range(N):
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
                for i in range(N1):
                    #cov[i,i:] = TCm1[i:].dot(TCm2[i].T)
                    cov[i,:] = TCm1[:].dot(TCm2[i].T)
                    #cov[i:,i] = cov[i,i:]
            else:
                for i in range(N1):
                    cov[i,:] = TCm1[[i]].dot(TCm2[:].T)
                    
            return cov, TCm1, TCm2
        
        elif TC1.ndim > 2:
            TCm1 = (TC1 - TC1.mean(axis=0)).T
            TCm2 = (TC2 - TC2.mean(axis=0)).T
            N1, N2 = TC1.shape[1], TC2.shape[1]
            cov  = ones((N1, N2))
            
            if N1 == N2:
                for i in range(N1):
                    cov[i,i:] = TCm1[i:].dot(TCm2[i].T)
                    cov[i:,i] = cov[i,i:]
            else:
                for i in range(N1):
                    cov[i,:] = TCm1[[i]].dot(TCm2[:].T)
                    
            return cov, TCm1, TCm2

        elif TC1.ndim == 1:
            if TC2.ndim == 2:
                V   = TC2 - TC2.mean(axis=1)[(slice(None,None,None),None)]
                V0  = TC1 - TC1.mean()
                Vss = sqrt(sum(V**2, axis=1))
                cov = V.dot(V0)
                cov/= Vss * sqrt(sum(V0**2))
                return cov
            else:
                return pearsonr(TC1,TC2)[0]


def fPCA(TC):
    from scipy.linalg import eigh
    CO = fCovariance(TC)[0]
    V,PC = eigh(CO)
    iPC, V = sortBy(V, inverse=True)
    return V, PC.T[iPC]


def matricesCorrelation(M1, M2, corrFunc=pearsonr, finite= False, posit=False, avg=True, **kwa):
    ''' Return the correlation between the two matrices, defined by the mean across
    correlations for each columns. If posit=True, the correlation is done for positive matrices.
    '''
    if M1.shape != M2.shape:
        print('Functional Connectomes shape are not the same')
        return 0
    
    N = M1.shape[0]
    corr = zeros(N)
    
    if corrFunc == corrcoef:
        if posit:
            for c in range(N):
                corr[c] = corrcoef(where(M1[c]<0, 0, M1[c]), where(M2[c]<0, 0, M2[c]))[0,1]
        else:
            for c in range(N):
                corr[c] = corrcoef(M1[c], M2[c])[0,1]
    elif corrFunc == pearsonr:
        if posit:
            for c in range(N):
                corr[c] = pearsonr(where(M1[c]<0, 0, M1[c]), where(M2[c]<0, 0, M2[c]))[0]
        else:
            for c in range(N):
                corr[c] = pearsonr(M1[c], M2[c])[0]
    else:
        if posit:
            for c in range(N):
                corr[c] = corrFunc(where(M1[c]<0, 0, M1[c]), where(M2[c]<0, 0, M2[c]), **kwa)
        else:
            for c in range(N):
                corr[c] = corrFunc(M1[c], M2[c], **kwa)
            
    if finite:
        corr = where(isfinite(corr), corr, 0.)

    if avg:
        return corr.mean()
    else:
        return corr


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
        for i in range(dim1):
            whiteValues[updateTable[(i+1) % lenT]] = randn(dim2)
            data[i] = (whiteValues.sum(0) + randn(dim2)) / magic_number
    else:
        whiteValues = randn(sTree)
        data = zeros(dim1)
        for i in range(dim1):
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
        for i in range(len(setP)):
            for j in range(i+1, len(setP)):
                simEF = miea(corr[setP[i]][:,setP[j]])
                #simEF = corr[setP[i]][:,setP[j]].max()
                if simEF > sim:
                    indE = i
                    indF = j
                    sim = simEF
        sim_final = sim
        
        if sim_final > sim_coef:
            try:
                for i in range(len(setP[indF])):
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
        f_combi = lambda L: list(combinations(list(range(L)), 2))[slice(L - 1)]
        acc = False
    else:
        f_combi = lambda L: list(combinations(list(range(L)), 2))
    
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
    to_test = list(range(len(SS)))
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


def randomMatrixPermutations(mat):
    """Return a matrix with elements of a triangle permuted without saving the degrees.
    """
    
    s = mat.shape
    nruter = array(mat)
    indices = find( triu( ones(s), k=1))    
    aleaInd = permutation(indices)
    nruter[unravel_index(indices, s)] = nruter[unravel_index(aleaInd, s)]
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
    
    nruter[unravel_index(Sind, s)] = nruter[unravel_index(Salea, s)]
    nruter[unravel_index(Mind, s)] = nruter[unravel_index(Malea, s)]
    nruter[unravel_index(Iind, s)] = nruter[unravel_index(Ialea, s)]
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

    print('randomized network #',1)
    srand = nullUndirectedNetwork(srand)
    n_1_2 = binned_srand_internal(srand,bedg1,bedg1)

    aver_n_1_2_sym = n_1_2
    aver_sq_n_1_2_sym = n_1_2 ** 2
    
    for k in range(1, Nstat):
        print('randomized network #',k+1)
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


def triSupPearsonCor(mat1, mat2, finite= False, norm=False, mask=None):
    if mask == None:
        mask = triSup(mat1, ind=True)
    
    corr = corrcoef(mat1.take(mask), mat2.take(mask))[0,1]
    
    if finite:
        return where(isfinite(corr), corr, 0.)
    else:
        return corr
    
    

def Wcorrelation(TC, ksi, posit=False):
    ''' Find the correlation between a Time Course 'TC' and different spatial patterns 'ksi'
    using the adjoint vectors.
    Returns the weights W and the reconstructed time course RTC.
    '''
    adj = adjointVectors(ksi)
    W = weights(TC, adj)
    W[isfinite(W) == 0] = 0
    if posit:
        W[W<0] = 0
    RTC = dot(W, ksi)
    return W, RTC

    
def goodnessOfFit(TC, RTC):
    '''Return the Goodness Of Fit between the Time Course and the Reconstructed Time Course
    reconstructed using the weighs of the adjoint vectors.
    '''
    difTC = TC-RTC
    GOF = 1.- (difTC*difTC).sum(0) / (TC*TC).sum(0)
    GOF[isfinite(GOF) == 0] = 0
    return GOF


def Kmeans(FC, n=10):
    '''kmeans clustering.
    n: number of clusters.
    '''
    from scipy.cluster.vq import kmeans
    from Tools.functions import sortBy
    Sk = kmeans(FC, n)[0]
    return Sk[sortBy(abs(Sk).mean(1), inverse=1)[0]]


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
        for i in range(0, FCs.shape[0]):
            FCs[i] = fPearsonCorrelation(TC[i*jump : i*jump + lWin], **kwa).take(mask)
        
    else:
        if nodes == None:
            FCs = zeros(((T - lWin) / jump, N, N))
            for i in range(0, FCs.shape[0]):
                FCs[i] = fPearsonCorrelation(TC[i*jump : i*jump + lWin], **kwa)
        
        else:
            try:
                n = len(nodes)
                FCs = zeros(((T - lWin) / jump, n, N))
                for i in range(0, FCs.shape[0]):
                    FCs[i] = fPearsonCorrelation(TC[i*jump : i*jump + lWin, nodes], 
                                                TC[i*jump : i*jump + lWin], **kwa)
                
            except:
                FCs = zeros(((T - lWin) / jump, N))
                for i in range(0, FCs.shape[0]):
                    FCs[i] = fPearsonCorrelation(TC[i*jump : i*jump + lWin, nodes], 
                                                TC[i*jump : i*jump + lWin].T, **kwa)        
        
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
            #mask = range(n * N)
            
        nruter = zeros((L1,L2))
        if not comp:
            for i in range(L1):
                for j in range(i,L1):
                    nruter[i,j] = triSupPearsonCor(FCs1[i,:,:], FCs1[j,:,:], mask=mask, **kwa)
            for i in range(L1):
                nruter[i:,i] = nruter[i,i:]
        
        else:
            for i in range(L1):
                for j in range(L2):
                    nruter[i,j] = triSupPearsonCor(FCs1[i,:,:], FCs2[j,:,:], mask=mask, **kwa)
                    
    elif n != N:
        nruter = zeros((L1,L1))
        for i in range(L1):
            for j in range(i,L1):
                nruter[i,j] = triSupPearsonCor(FCs1[i,:,:], FCs1[j,:,:], mask=mask, **kwa)
        for i in range(L1):
            nruter[i:,i] = nruter[i,i:]
        
    
    elif len(mask) == N:
        nruter = zeros((N,L1,L1))
        for n in range(N):
            nruter[n,:,:] = fPearsonCorrelation(FCs1[:,n,:].T, **kwa)
        
    else:
        nruter = zeros((L1,L1))
        for i in range(L1):
            for j in range(i,L1):
                nruter[i,j] = triSupPearsonCor(FCs1[i,:,:], FCs1[j,:,:], mask=mask, **kwa)
        for i in range(L1):
            nruter[i:,i] = nruter[i,i:]
    
    return nruter


