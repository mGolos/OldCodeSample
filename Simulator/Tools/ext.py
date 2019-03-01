#!/usr/bin/env python
#-*- coding:Utf-8 -*-

from numpy import array, save as np_save, load as np_load, where, linspace, array
from ast import literal_eval
import os



def humanSize(nbytes):
    suffixes = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']
    if nbytes == 0: return '0 B'
    i = 0
    while nbytes >= 1024 and i < len(suffixes)-1:
        nbytes /= 1024.
        i += 1
    f = ('%.2f' % nbytes).rstrip('0').rstrip('.')
    return '%s %s' % (f, suffixes[i])


def dir2array(adress, mmap_mode=None, dic=False, **kwa):
    ll = os.listdir(adress)
    ll.sort()
    if '.' in ll:
        ll.remove('.')
    
    if dic == False:
        nruter = []
        for l in ll:
            nruter.append(data2array(adress+'/'+l, mmap_mode=mmap_mode, **kwa))
        
    else:
        nruter = {}
        for l in ll:
            nruter[l.rsplit('.')[-2]] = data2array(adress+'/'+l, mmap_mode=mmap_mode, **kwa)
    
    return nruter


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
                if other == None:
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


def loadTendance(Folder):
    ''' Return tendances.txt of patterns in a ndarray.'''
    filetxt = Folder + '/tendance.txt'
    try:
        ofi = open(filetxt, 'r')
        tendance = array([line.split() for line in ofi]).astype(float)
        ofi.close()
        return tendance
    except:
        print '%s does not exists' %filetxt
        return array(0)


#def findParameters(adress, p1=1, p2=3, dOf="file", dic=True):
    #''''''
    #l0 = os.listdir(adress)
    #for i in range(len(l0))[::-1]:
        #if (dOf == "file" and (not ('patterns' in l0[i] or 'BO' in l0[i])))\
        #or (dOf == "dir"  and isfile(adress +'/'+ l0[i]))\
        #or l0[i].startswith('.'):
            #l0.remove(l0[i])

    #if l0 == []:
        #return l0

    #lp1,lp2 = [],[]
    #for i in range(len(l0)):
        #splited = l0[i].replace('.npy','').rsplit('_')
        #if splited[p1+1] not in lp1:
            #lp1.append( splited[p1+1] )
        #if splited[p2+1] not in lp2:
            #lp2.append( splited[p2+1] )
    #lp1 = array(lp1, dtype=float);  lp1.sort()
    #lp2 = array(lp2, dtype=float);  lp2.sort()
    
    #if dic:
        #return {splited[p1]:lp1, splited[p2]:lp2}
    #else:
        #return lp1, lp2
    
    
def findParameters(adress, spl='_', stwith='', edwith='', atype=str, dic=True):
    '''Find the parameters of files inside a directory using a specific character (default:"_").
    An starting (stwith) or ending (edwith) string can be defined for specification.
    The fina type can be define by normal types or something like "lambda x: '%.3f'%float(x)".
    '''
    keys = []
    nruter = []
    l0 = os.listdir(adress)
    
    for f in l0:
        if not f.startswith(stwith) or not f.endswith(edwith):
            continue
            
        fs = f[:-len(f.split('.')[-1]) -1].split(spl)
        for ip, p in enumerate(fs):
            try:
                float(p)
                try:
                    nruter[keys.index(fs[ip-1])].add(atype(p))
                except:
                    keys.append(fs[ip-1])
                    nruter.append(set([atype(p)]))
            except:
                pass
                
    for i in range(len(nruter)):
        nruter[i] = sorted(nruter[i])
        #nruter[i] = natSort(nruter[i])
    
    if not dic:
        return keys, nruter
    else:
        dic = {}
        for i in range(len(keys)):
            dic[keys[i]] = nruter[i]
        return dic


def paramExplo(f, ax='x', x=[0.,1.], y=[0.,1.], nb=[10,10], val0=1., revert=False, loop=False):
    '''Explore the parameter in a dicotomic way.
    '''
    
    Ai, Af = [x,y][ax!='x']
    Bs = [y,x][ax!='x']
    if len(Bs) == 2:
        Bs = linspace(Bs[0], Bs[1], nb[1])
    else:
        Bs = array(Bs)
        
    dots = []
    for B in Bs:
        A0, A1, val = Ai, Af, val0
        d = []
        Am = [A1, A0][revert]

        while len(d) < nb[0]:
            #if A0 == A1 and not loop:
                #return paramExplo(f, ax=ax, x=x, y=y, nb=nb, val0=val0, revert=not revert, loop=1)
            
            args = [Am, B][::(-1)**(ax!='x')]
            val = f(*args)
            d.append(args + [val])

            if (val == val0) != (revert):
                A0 = Am
            elif (val == val0) != (not revert):
                A1 = Am
                
            Am = (A0 + A1) / 2.

        dots.extend(d)
    return array(dots)


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
  

def Pdir(sdir='', host=None):
    '''Return the 'Programmation' directory depending on the computer.'''
    if host == None:
        import socket
        host = socket.gethostname()     
    
    if host == 'HP': #Kubuntu
        if 'Simulation' in sdir:
            prog = '/home/golos/Main/'
        elif 'Main' in sdir:
            prog = '/home/golos/'
        else:
            prog = '/home/golos/Main/Programs/'
        return prog + sdir.replace('\\',"/")
        
    elif host == 'Arkas': #W8
        if 'Simulation' in sdir:
            prog = 'D:\\'
        else:
            prog = 'C:\\Users\\Mathieu\\Documents\\'
        return prog + sdir.replace('/',"\\")
        
    elif host == 'Sleipnir': #W8
        prog = 'D:\\Documents\\Main\\Programs\\'
        return prog + sdir.replace('/',"\\")
        
    elif host == 'loginuser.cluster': #Linux
        prog = '/home/mathieuG/'
        return prog + sdir.replace('\\',"/")
        
    elif host == 'Windows':
        print 'Computer not defined.'
        prog = '..\\'
        return prog + sdir.replace('/',"\\")
        
    else:
        print 'Computer not defined.'
        prog = '../'
        return prog + sdir.replace('\\',"/")

            
def mkdir(adress):
    ''' Crée un dossier s'il n'existe pas.'''
    if not os.path.exists(adress):
        os.mkdir(adress)

def isdir(adress, f=None, disp=False):
    ''' Teste si un dossier existe.'''
    nruter = os.path.exists(adress)
    if not nruter and disp:
        print 'Pas de dossier ' + adress
        if f != None:
            f()
    return nruter

isfile = os.path.isfile
            
adressExists = isdir

          
def forFreeSurfer(vert, tria, l_adr, r_adr=None):
    if r_adr == None:
        ofi = open(l_adr, 'w')
        ofi.write('%i  %i\n' %(len(vert),len(tria)))
        for i in xrange(len(vert)):
            ofi.write('%f  %f  %f  0\n' %(vert[i][0],vert[i][1],vert[i][2]))
        for i in xrange(len(tria)):
            ofi.write('%i  %i  %i  0\n' %(tria[i][0],tria[i][1],tria[i][2]))
        ofi.close()
        
    else:
        Lv = len(vert)/2
        Lt = len(tria)/2
        lvert = vert[:Lv]
        rvert = vert[Lv:]
        ltria = tria[:Lt]
        rtria = tria[Lt:] - Lv
        for adr, V, T in [l_adr, lvert, ltria], [r_adr, rvert, rtria]:
            ofi = open(adr, 'w')
            ofi.write('%i  %i\n' %(len(V),len(T)))
            for i in xrange(len(V)):
                ofi.write('%f  %f  %f  0\n' %(V[i][0],V[i][1],V[i][2]))
            for i in xrange(len(T)):
                ofi.write('%i  %i  %i  0\n' %(T[i][0],T[i][1],T[i][2]))
            ofi.close()
    
    
def fromFreeSurfer(lhS_adr, rhS_adr, excent=125.):
    from numpy import append
    lhS = data2array(lhS_adr).tolist()
    rhS = data2array(rhS_adr).tolist()
    nl_vert, nl_tria = array(lhS[1]).astype(int)
    nr_vert, nr_tria = array(rhS[1]).astype(int)
    del lhS[:2], rhS[:2]
    
    l_vert = array(lhS[0:nl_vert]).astype(float)
    r_vert = array(rhS[0:nr_vert]).astype(float)
    vert = array(append(l_vert, r_vert, axis=0).tolist()).astype(float)[:,:3]
    vert[:nl_vert, 0] -= excent
    vert[nl_vert:, 0] += excent
    
    l_tria = array(lhS[nl_vert:nl_vert+nl_tria]).astype(int)
    r_tria = array(rhS[nr_vert:nr_vert+nr_tria]).astype(int) + nl_vert
    tria = array(append(l_tria, r_tria, axis=0).tolist()).astype(int)[:,:3]
    
    return vert, tria


def natSort(lst):
    #return lst.sort(key=lambda x:map(int, str(float(x)).split(".")))
    return sorted(lst, key = lambda x: (float(x.partition(' ')[0]) if x[0].isdigit() else float('inf'), x))


