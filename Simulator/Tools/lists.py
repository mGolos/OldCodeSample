#!/usr/bin/env python
#-*- coding:Utf-8 -*-

from pylab import zeros


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


def switch(l, i1, i2):
    '''Switch two position of a list.'''
    tmp = l[i1]
    l[i1] = l[i2]
    l[i2] = tmp
    
    
def sortBy(root, inverse=False):
    '''Return the indice list to rearange 'root' and the rearange array.
    '''
    rootb = array(root)
    si = rootb.argsort()
    rootb.sort()
    if inverse:
        return si[::-1], rootb[::-1]
    else:
        return si, rootb
    
    
def sortByIrregularDimension(l, inverse=False):
    L0, L1 = findMaxDim(l)
    with0 = zeros((L0,L1))
    
    for i in range(L0):
        with0[i, :len(l[i])] = l[i]
        
    return sortBy(with0.sum(1), inverse=inverse)


def findMaxDim(l):
    L0, L1 = len(l), 0
    for i in range(L0):
        test = len(l[i])
        if test > L1:
            L1 = test
    return L0, L1


def indicesFrom2DList(AllClust):
    '''Return a list of indices.'''
    iClust = []
    for ic in range(len(AllClust)):
        iClust.append(list(range(len(AllClust[ic]))))
    return iClust


def newList(l0):
    '''Return a list of list and not adresses. (so far for 2D).'''
    l1 = []
    try:
        for i in range(len(l0)):
            l1.append(l0[i][:])
        return l1
    
    except :
        return list(l0)
    

def dimension(a):
    '''Return the dimensions of a list/array of lists/arrays.'''
    pos = a[0]
    dim = 1
    while True:
        try:
            pos = pos[0]
            dim += 1
        except:
            break
    return dim


def dimensions(a):
    '''Return the dimensions of a list/array of lists/arrays.'''
    pos = a
    dim = [len(pos)]
    while True:
        try:
            pos = pos[0]
            dim.append(len(pos))
        except:
            break
    return dim


def profondeur(dic):
    dim = 0
    try:
        pos = dict(dic)
    except:
        return dim
    while True:
        try:
            pos = pos[pos.keys()[0]]
            dim += 1
        except:
            break
    return dim


def flattenDict(dic, ndim=None):
    if ndim == None:
        ndim = profondeur(dic)
    if ndim == 0:
        return [[dic]]
    else:
        l = []
        for k in sorted(dic.keys()):
            for j in flattenDict(dic[k], ndim=None):  
                l.append([k])
                l[-1].extend(j)
        return l    

