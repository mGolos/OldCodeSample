#!/usr/bin/env python
#-*- coding:Utf-8 -*-

from pylab import array, diag, sort, where, dot, zeros, isfinite
from pylab import arange, unravel_index, r_, sqrt, empty, NaN
from pylab import find, ones, concatenate, newaxis, c_

from Tools.lists import dimension, dimension, sortBy


def HRe_from_LDown(L):
    '''Return the array of indices to reorder a matrix of high dimension H
    using the order of a matrix of lower dimension L.
    '''
    r = array(L)
    s = 0
    for k in range(max(L)+1):
        n = (L==k).sum()
        r[L==k] = s + arange(n)
        s += n
    nruter = arange(len(L))
    nruter[r] = arange(len(L))
    return nruter


def invOrder(o1):
    '''Return the opposite law of reordering.'''
    o2 = arange(len(o1))
    o2[o1] = arange(len(o1))
    return o2

    
def seuil(m, s, val=0, over=False):
    n = array(m)
    if over:
        n[n > s] = val
    else:
        n[n < s] = val
       
    return n


def seuilLiens(m, s):
    L = m.shape[0] * m.shape[1]
    h = sort(m.reshape(L))
    v = h[int(s * (L-1))]
    return where(m < v, 0., m)


def fProximity(A, B=None, zeroDiag=True):
    ''' Return the proximity (similarity x correlation) as :
    - 2D nparray scalar between 2D nparray vectors (filled with zeros for diagonal and symetrix terms)
    - 1D nparray scalar between 1D nparray vector and 2D nparray vectors
    - 2D nparray scalar between 2D nparray vectors
    '''
    sA = A.shape
    sB = B.shape
    
    if B == None:
        corr = zeros((sA[0], sA[0]))
        for i in range(sA[0]):
            corr[i, i+1:] = fProximity(A[i], A[i+1:])
        return corr
    
    elif A.ndim == 1:
        dif = 1.- abs(A-B).sum(axis=-1) / (1.*sA[0])
        sim = B.dot(A) / (A**2).sum(axis=-1)**(0.5) / (B**2).sum(axis=-1)**(0.5)
        return where(isfinite(sim), dif * sim, dif)
    
    elif B.ndim == 1:
        dif = 1.- abs(A-B).sum(axis=-1) / (1.*sB[0])
        sim = A.dot(B) / (A**2).sum(axis=-1)**(0.5) / (B**2).sum(axis=-1)**(0.5)
        return where(isfinite(sim), dif * sim, dif)
    
    else:
        corr = zeros((sA[0], sB[0]))
        for i in range(sA[0]):
            corr[i] = fProximity(A[i], B)
        return corr - zeroDiag * diag(diag(corr))


def distance(a,b):
    ''''''
    return abs(a-b).sum(axis=-1) / float(len(a))


def distanceBasedSim(a,b):
    ''''''
    return 1./ (1.+ distance(a,b))


def triHemi(mat, h=0, ind=False, anti=False, sup=True, h2flat=True):
    try:
        N = len(mat)
    except:
        N = mat
    if sup:
        fSI = triSup
    else:
        fSI = triInf
        
    if h == 0:
        return fSI(mat, ind=ind, anti=anti)
    
    elif h == 1:
        part = arange(N * N).reshape((N, N))
        indices = r_[fSI(part[:N/2, :N/2], anti=anti), 
                         fSI(part[N/2:, N/2:], anti=anti)]
        if ind:
            return indices
        else:
            return mat[unravel_index(indices, mat.shape)]
    
    elif h == 2:
        if sup:
            sli = getSlices(N, h=2)[2:4]
        else:
            sli = getSlices(N, h=2)[0:2]
        if ind:
            if h2flat:
                return arange(N * N).reshape((N, N))[sli].flatten()
            else:
                return arange(N * N).reshape((N, N))[sli]
        else:
            return mat[sli]

        
def meanHemisph(mat, h=0):
    '''Return the matrix average across intra or inter-hemispheres.
    h=0: Both hemispheres, complete slices,
    h=1: Intra-hemisphere slices,
    h=2: inter-hemisphere slices.
    '''
    sli = getSlices(len(mat), h)
    if sli == None:
        return mat
    else:
        return (mat[sli[0:2]] + mat[sli[2:4]]) *.5
    
    
def getSlices(N, h=0):
    '''Return slices to average across intra or inter hemispheres.
    N: number of nodes.
    '''
    if h== 0:
        return None
    elif h == 1:
        return slice(0,N/2), slice(0,N/2), slice(N-1,N/2-1,-1), slice(N-1,N/2-1,-1)
    elif h == 2:
        return slice(N/2,N), slice(0,N/2), slice(-N/2-1,-N-1,-1), slice(N-1,N/2-1,-1)


def TwoTri(m1, m2):
    if m1.ndim == 1:
        L = len(m1) * 2
        N = int(sqrt(L)) + 1
    else:
        N = len(m1)
    
    nruter = empty((N,N)) * NaN
    indS = triSup(nruter, ind=True)
    indI = triInf(nruter, ind=True)
    if m1.ndim == 1:        
        v1, v2 = m1, m2
    else:
        v1, v2 = m1.take(indS), m2.take(indS)
        
    nruter[unravel_index(indS, nruter.shape)] = v1
    nruter[unravel_index(indI, nruter.shape)] = v2
    return nruter


def triSup(mat, ind=False, anti=False):
    try:
        N = len(mat)
    except:
        N = mat
    IND = zeros(N*(N-1)/2, dtype=int)
    s = 0
    if anti:
        for i in range(N):
            ds = N -i -1
            IND[s:s+ds] = arange(i, i+ds*N, N)
            s += ds
    else:
        for i in range(N):
            ds = N -i -1
            IND[s:s+ds] = arange(1+i*(N+1), (i+1)*N)
            s += ds
    if ind:
        return IND
    else:
        return mat.take(IND)


def triSupWithDiag(mat, ind=False):
    try:
        N = len(mat)
    except:
        N = mat
    IND = zeros(N*(N+1)/2, dtype=int)
    s = 0
    
    for i in range(N):
        ds = N - i
        IND[s:s+ds] = arange(i*(N+1), (i+1)*N)
        s += ds
        
    if ind:
        return IND
    else:
        return mat.take(IND)


def triInf(mat, ind=False, anti=False):
    try:
        N = len(mat)
    except:
        N = mat
    IND = []
    s = 0
    if anti:
        for i in range(N):
            IND[s:s+i] = arange(N+i*(N-1), (i+1)*N)
            s += i
    else:
        for i in range(N):
            ds = N -i -1
            IND[s:s+ds] = arange(i+(i+1)*N, i+N*N, N)
            s += ds
    if ind:
        return IND
    else:
        return mat.take(IND)
    
    
def triIndices(size, nodes=[], nodes2=None):
    '''Return the mask using the selected nodes (list of integer) to obtain the corresponding lines
    to apply on a vectorized triangle matrix.
    Apply the mask like vect[mask] and not vect.take(mask).
    '''
    i,j = triToMat(size)
    mask = zeros(len(i)).astype(bool)
    if nodes2 == None:
        for n in nodes:
            mask = mask | (i==n) | (j==n)
    else:
        for n in nodes:
            for n2 in nodes2:
                mask = mask | ((i==n) & (j==n2))
                mask = mask | ((i==n2) & (j==n))
    return mask
        

def triToMat(N, tri=True):
    '''Return the indices i,j from the upper right triangle array (triSup) as a 2*L np.array.
    '''
    if tri:
        #L = N * (N-1) / 2
        #nruter = zeros((2,L), dtype=int)
        #k = 0
        #for i in range(N-1):
            #for j in range(i+1, N):
                #nruter[:, k] = [i,j]
                #k += 1
        nruter = array(unravel_index(triSup(N, ind=1), (N,N)))
    else:
        #nruter = zeros((2,N**2), dtype=int)
        #k = 0
        #for i in range(N):
            #for j in range(N):
                #nruter[:, k] = [i,j]
                #k += 1
        nruter = array(unravel_index(arange(N**2), (N,N)))
            
    return nruter


def adjointVectors(ksi):
    ''' Return the adjoint vectors of ksi where ksi is composed of m vector of size > 1
    '''
    from pylab import matrix, eye
    A = matrix( dot(ksi, ksi.T) )
    A = A + 1e-1 * eye(A.shape[0])
    try:
        #A = matrix( dot(ksi, ksi.T) ).I
        A = A.I
    except:
        A = zeros(dot(ksi,ksi.T).shape)
        print('Non inversible dot(ksi, ksi.T) building adjoint vectors where ksi are the patterns')
    return dot(A, ksi)


def weights(TC, adj):
    ''' Return the weights of each m vectors for each t time : W[t][m]
    '''
    T = TC.shape[0]
    W = zeros((T, adj.shape[0]))
    for t in range(T):
        W[t] = dot(TC[t], adj.T)# / norm(TC[t])
    return W


def reorderByMaxs(corrAB, corrAC=None, iA=None):
    L = len(corrAB)
    
    if corrAC is None:
        if iA is None:
            iA = sortBy(corrAB.max(0), inverse=True)[0]
            retA = True
        else:
            retA = False
            
        iB = zeros(L, dtype=int)
        tmp = list(range(L))
        order = sortBy(corrAB[iA].max(1), inverse=1)[0]
        for i in range(L):
            iB[order[i]] = tmp[corrAB[iA][order[i]][tmp].argmax()]
            tmp.remove(iB[order[i]])
            
        if retA: return iA, iB
        else:    return iB
        
    else:
        if iA is None:
            ma = array([corrAB.max(1), corrAC.max(1)])
            iA = sortBy(where(ma.argmax(0), ma[1], ma[0]), inverse=True)[0]
            retA = True
        else:
            retA = False
            
        iB = corrAB[iA].argmax(1)
        iC = corrAC[iA].argmax(1)
        if not len(set(iB)) == len(iB):
            print("To implement")
        if not len(set(iC)) == len(iC):
            print("To implement")
            
        if retA: return iA, iB, iC
        else:    return iB, iC        

