#!/usr/bin/env python
#-*- coding:Utf-8 -*-
''' Affiche les patterns, leur densité, leur répartition ou leur tendance en fonction des jeux de paramètres
    @author: mathieu.golos@gmail.com
'''

'#################################################  LIBRARIES  ##################################################'
import os
import copy
import Tools.ext as Te
import Tools.matrices as Tm
import Tools.functions as Tf
from pylab import figure, array, zeros, arange, get_cmap, corrcoef, where, average, size, linspace, pi, cm
from pylab import log, imread, subplot, rand, unravel_index, ones, sqrt, norm, newaxis, isfinite, diag
from numpy import load as np_load, save as np_save
from pickle import load as pk_load, dump as pk_dump
from matplotlib import colors, rcParams, gridspec
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.ticker import FormatStrFormatter
from scipy.stats import ss
from Tkinter import *
from time import time
rcParams.update({'font.size': 9})

CNTRF = 420


'##############################################  BASIC FUNCTIONS  ###############################################'
def textentry(root, variable, label, width):
    global row_counter
    l= Label(root, text=label)
    l.grid(column=0, row=row_counter, sticky='w')
    widget = Entry(root, textvariable=variable, width=width)
    widget.grid(column=1, row=row_counter)
    row_counter += 1
    return widget

def Quit():
    root.quit()
    W_frame.destroy()

def listaDef(e):
    p['p1'].set(float(lista.get(lista.curselection())))
    Show()

def listbDef(e):
    p['p2'].set(float(listb.get(listb.curselection())))
    Show()

def listcDef(e):
    p['cmap'].set(str(listc.get(listc.curselection())))
    Show()

def dirPriPlusSub():
    ''' Retourne l'adresse où les différentes jeux de patterns sont présents
    '''
    if 'Attr' in p['dir_sub'].get():
        if 'X' in p['dir_sub'].get():
            dir_ = p['dir_pri'].get() + '/At' + p['attr'].get() + '_X/'
        else:
            dir_ = p['dir_pri'].get() + '/At' + p['attr'].get() + '_S/'
    else:
        dir_ = p['dir_pri'].get() + p['dir_sub'].get() + '/'
    return dir_

def dirTri(p1, p2):
    ''' Retourne une liste d'adresses complètes menants aux patterns des paramètres choisis
    '''
    global p1str, p2str
    dir_ = dirPriPlusSub()
    test = os.listdir(dir_)[0]
    p1str= 'Norm' * ("Norm" in test) \
         + 'Psi' * ("Psi" in test) \
         + 'P' * ("P_" in test)
    p2str= 'Gamma' * ("Gamma" in test) \
         + 'gamma' * ("gamma" in test) \
         + 'DensMoy' * ("DensMoy" in test) \
         + 'ThetaRef' * ("ThetaRef" in test) \
         + 'G' * ("G_" in test) \
         + 'f0' * ("f0" in test) \
         + 'Omega' * ("Omega" in test)
     
    params = []
    for i in xrange(len(p1)):
        params.append([])
        for j in xrange(len(p2)):
            if len(str(lp2[1]-lp2[0]))-2 > 2:
                params[i].append( p1str + '_%.2f_'%p1[i] + p2str + '_%.3f'%p2[j] )
            else:
                params[i].append( p1str + '_%.2f_'%p1[i] + p2str + '_%.2f'%p2[j] )
    return params


def paramLists():
    ''' Retourne les adrresses des fichiers "patterns.npy" et "tendances.txt".'''
    global p1str, p2str
    dir_ = dirPriPlusSub()
    parameters = Te.findParameters(dir_, atype=float)

    pK = parameters.keys()
    p1str= 'Norm' * ("Norm" in pK) \
         + 'Psi' * ("Psi" in pK) \
         + 'P' * ("P" in pK)
    p2str= 'Gamma' * ("Gamma" in pK) \
         + 'gamma' * ("gamma" in pK) \
         + 'DensMoy' * ("DensMoy" in pK) \
         + 'ThetaRef' * ("ThetaRef" in pK) \
         + 'G' * ("G_" in pK) \
         + 'f0' * ("f0" in pK) \
         + 'Omega' * ("Omega" in pK)
     
    p1lst = array(parameters[p1str])
    p2lst = array(parameters[p2str])
    return p1lst, p2lst


def adressFromPar(p1, p2, pre='', post=''):
    try:
        order = len(str(p2lst[1] - p2lst[0])) -2 > 2
    except:
        order = 0

    if order:
        return pre + '%s_%.2f_%s_%.3f'% (p1str, p1, p2str, p2) + post
    else:
        return pre + '%s_%.2f_%s_%.2f'% (p1str, p1, p2str, p2) + post


def ShowAtr():
    if not "Attr" in p['dir_sub'].get():
        p['dir_sub'].set("/Attr_S")
    Show()

def mapCanvas(o):
    if not p['C'][o].get_tk_widget().winfo_ismapped():
        for i in range(len(p['C'])):
            p['C'][i].get_tk_widget().pack_forget()
            p['T'][i].pack_forget()
        p['C'][o].get_tk_widget().pack()
        p['T'][o].pack()

def cmap_powerlaw_adjust(cmap, a):
    if a < 0.:
        return cmap
    cdict = copy.copy(cmap._segmentdata)
    fn = lambda x : (x[0]**a, x[1], x[2])
    for key in ('red','green','blue'):
        cdict[key] = map(fn, cdict[key])
        cdict[key].sort()
        assert (cdict[key][0]<0 or cdict[key][-1]>1), \
            "Resulting indices extend out of the [0, 1] segment."
    return colors.LinearSegmentedColormap('colormap', cdict, 1024)

def cmap_cutter(cmap, part):
    cdict = copy.copy(cmap._segmentdata)
    n = len(cdict['red']) / 2 + 1
    ndict = {}
    for key in ('red','green','blue'):
        ndict[key] = zeros((n,3))
        for i in range(n):
            ndict[key][i][0]  = i / (n-1.)
            if part in ['up','U',0]:
                ndict[key][i][1:] = cdict[key][i+n-1][1:]
            elif part in ['down','D',1]:
                ndict[key][i][1:] = cdict[key][i][1:]
    return colors.LinearSegmentedColormap('colormap', ndict, 1024)

def cmap_center_adjust(cmap, vmin, vmax, center=0, relative=0):
    center_ratio = (center - vmin) / (vmax - vmin)
    a = log(abs(center_ratio)) / log(0.5)
    #if center <= vmin:
    if center_ratio <= 0.05:
        cmap.set_under('w', 1.0)
        if relative:
            return cmap_powerlaw_adjust( cmap_cutter(cmap, part='up'), vmax)
        else:
            return cmap_cutter(cmap, part='up')
    #elif vmax <= center:
    elif 0.95 <= center_ratio:
        cmap.set_over('w', 1.0)
        if relative:
            return cmap_powerlaw_adjust( cmap_cutter(cmap, part='down'), -1./ vmin)
        else:
            return cmap_cutter(cmap, part='down')
    else:
        return cmap_powerlaw_adjust(cmap, a)

def parcellation(vec, color, maxs=None):
    N = len(vec)
    if N == 998:
        vec66  = zeros(66)
        for n in range(N):
            vec66[D998_D66[n]] += vec[n]
        vec = vec66 / AvgN_D66
    if maxs == None:
        vmin = vec.min()
        vmax = vec.max()
        if vmin == vmax:
            vmin = vmin * (vmax < 0)
            vmax = vmax * (vmax > 0) + 1.* (vmax==0)
        cRatio= - vmin / (vmax - vmin)
        bord  = where(cRatio < 0.9, vmax, vmin)
    else:
        vmin = maxs[0]
        vmax = maxs[1]
        bord = vmax
    int_R = bord * int_0 + sum([int_N[i]    * (vec[i])    for i in range(33)])
    int_L = bord * int_0 + sum([int_N[-1-i] * (vec[33+i]) for i in range(33)])
    ext_R = bord * ext_0 + sum([ext_N[i]    * (vec[i])    for i in range(33)])
    ext_L = bord * ext_0 + sum([ext_N[-1-i] * (vec[33+i]) for i in range(33)])

    cmap = cmap_center_adjust(color, vmin, vmax, center=0, relative=0)
    grid[0].clear()
    grid[1].clear()
    grid[2].clear()
    grid[3].clear()
    im1 = grid[0].imshow(int_L[:,::-1], cmap=cmap, norm = colors.Normalize(vmin=vmin, vmax=vmax, clip = False))
    im2 = grid[1].imshow(int_R,         cmap=cmap, norm = colors.Normalize(vmin=vmin, vmax=vmax, clip = False))
    im3 = grid[2].imshow(ext_L[:,::-1], cmap=cmap, norm = colors.Normalize(vmin=vmin, vmax=vmax, clip = False))
    im4 = grid[3].imshow(ext_R,         cmap=cmap, norm = colors.Normalize(vmin=vmin, vmax=vmax, clip = False))
    grid.cbar_axes[0].colorbar(im1)
    grid.axes_llc.set_xticks([])
    grid.axes_llc.set_yticks([])
    p['C'][-1].draw()
    #p['F'][-1].savefig('test.png', dpi=300)

def clickParcellation(event, mat, color, maxs=None):
    mapCanvas(-1)
    print 'pattern %i' %(xdata + 0.5)
    parcellation(mat[int(xdata + 0.5)], color, maxs=maxs)
    p['C'][-1].get_tk_widget().bind("<Button-2>", lambda event: Show())

def permitSwitches(mat, color, maxs=None):
    p['C'][o].get_tk_widget().bind("<Button-2>", \
                lambda event: clickParcellation(event, mat, color, maxs=maxs))

def getStatusXY(event):
    global xdata, ydata
    xdata = event.xdata
    ydata = event.ydata

def Abreviations(abrev):
    abre_ind = ['r'+R for R in abrev]
    for i in xrange(33):
        abre_ind.append('l'+ abrev[::-1][i])
    return abre_ind

def sort2by1(one,two,inverse=False):
    oneb = array(one)
    si = oneb.argsort()
    oneb.sort()
    if inverse:
        return oneb[::-1], two[si[::-1]]
    else:
        return oneb, two[si]

def fastPearsonCorrelation1D(V0, V):
    ''' VO is one dimensional and V two dimensional
    '''
    V  -= V.mean(axis=1)[(slice(None,None,None),None)]
    V0 -= V0.mean()
    Vss = sqrt(ss(V, axis=1))
    cov = V.dot(V0)
    cov/= Vss * sqrt(ss(V0))
    return cov

def similAndNotCorr(V0, V):
    sim = zeros(V.shape[0])
    for i in xrange(V.shape[0]):
        sim[i] = V0.dot(V[i]) / norm(V0) / norm(V[i])
    return sim

def fProximity(a,b=None):
    ''' Return the proximity (similarity x correlation) as :
    - 2D nparray scalar between 2D nparray vectors (filled with zeros for diagonal and symetrix terms)
    - 1D nparray scalar between 1D nparray vector and 2D nparray vectors
    - 2D nparray scalar between 2D nparray vectors
    '''
    if b == None:
        corr = zeros((len(a), len(a)))
        for i in range(len(a)):
            corr[i,i+1:] = fProximity(a[i], a[i+1:])
        return corr
    elif a.ndim == 1:
        dif = 1.- abs(a-b).sum(axis=-1) / (1.*len(a))
        sim = b.dot(a) / (a**2).sum(axis=-1)**(0.5) / (b**2).sum(axis=-1)**(0.5)
        return where(isfinite(sim), dif * sim, dif)
    elif b.ndim == 1:
        dif = 1.- abs(a-b).sum(axis=-1) / (1.*len(b))
        sim = a.dot(b) / (a**2).sum(axis=-1)**(0.5) / (b**2).sum(axis=-1)**(0.5)
        return where(isfinite(sim), dif * sim, dif)
    else:
        corr = zeros((len(a), len(b)))
        for i in range(len(a)):
            corr[i] = fProximity(a[i], b)
        return corr

def fCosineSimilarity(a,b=None):
    ''' Return the cosine similarity as :
    - 2D nparray scalar between 2D nparray vectors (filled with zeros for diagonal and symetrix terms)
    - 1D nparray scalar between 1D nparray vector and 2D nparray vectors
    - 2D nparray scalar between 2D nparray vectors
    '''
    if b == None:
        corr = zeros((len(a), len(a)))
        for i in range(len(a)):
            corr[i,i+1:] = fCosineSimilarity(a[i], a[i+1:])
        return corr
    elif a.ndim == 1:
        sim = b.dot(a) / (a**2).sum(axis=-1)**(0.5) / (b**2).sum(axis=-1)**(0.5)
        return where(isfinite(sim), sim, 0)
    elif b.ndim == 1:
        sim = a.dot(b) / (a**2).sum(axis=-1)**(0.5) / (b**2).sum(axis=-1)**(0.5)
        return where(isfinite(sim), sim, 0)
    else:
        corr = zeros((len(a), len(b)))
        for i in range(len(a)):
            corr[i] = fCosineSimilarity(a[i], b)
        return corr

def fOverDistance(a,b=None):
    ''' Return similarity as 1 / (1 + diff) defined on [0.5 ; 1]:
    - 2D nparray scalar between 2D nparray vectors (filled with zeros for diagonal and symetrix terms)
    - 1D nparray scalar between 1D nparray vector and 2D nparray vectors
    - 2D nparray scalar between 2D nparray vectors
    '''
    if b == None:
        corr = zeros((len(a), len(a)))
        for i in range(len(a)):
            corr[i,i+1:] = fOverDistance(a[i], a[i+1:])
        return corr
    elif a.ndim == 1:
        return len(a) / (len(a) + abs(a-b).sum(axis=-1))
    elif b.ndim == 1:
        return len(b) / (len(b) + abs(a-b).sum(axis=-1))
    else:
        corr = zeros((len(a), len(b)))
        for i in range(len(a)):
            corr[i] = fOverDistance(a[i], b)
        return corr


'##############################################  MAIN FUNCTIONS  ################################################'
def refreshLists():
    ''' Actualise les listes de paramètres de façon automatique selon le dossier de simulation donné
    et supprime certaines matrices sauvegarder pour économiser du temps de calcul lors de switch
    entre différents affichages pour un même type de simulation.
    '''
    global lista, listb, G_dens, G_patt, G_tend, G_Ohist, G_Dhist, lp1, lp2
    lp1,lp2 = paramLists()

    lista.destroy()
    listb.destroy()
    lista= Listbox(E2b_frame, width=10, exportselection=0);   lista.pack(side='left', pady=5)
    listb= Listbox(E2b_frame, width=10, exportselection=0);   listb.pack(side='left', pady=5)
    for nW in lp1:
        lista.insert(END, '%.3f' %nW)
    for nW in lp2:
        listb.insert(END, '%.3f' %nW)
    lista.bind('<Double-1>', listaDef)
    listb.bind('<Double-1>', listbDef)

    G_dens  = None
    G_patt  = None
    G_tend  = None
    G_Ohist = None
    G_Dhist = None
    if os.path.exists(p['dir_pri'].get() + '/Cand_S'): CS.config(state = NORMAL)
    else :                                             CS.config(state = DISABLED)
    if os.path.exists(p['dir_pri'].get() + '/Cand_X'): CX.config(state = NORMAL)
    else :                                             CX.config(state = DISABLED)
    if os.path.exists(p['dir_pri'].get() + '/Patt_S'): PS.config(state = NORMAL)
    else :                                             PS.config(state = DISABLED)
    if os.path.exists(p['dir_pri'].get() + '/Patt_X'): PX.config(state = NORMAL)
    else :                                             PX.config(state = DISABLED)

    isAttr = False
    for noise in [1,2,3,4,5,6,7,8,9,10,15,20,25,30,35,40,45,50]:
        if os.path.exists(p['dir_pri'].get() + '/At%.2i_S' %noise):
            globals()['At%.2i' %noise].config(state = NORMAL)
            isAttr = True
        else :
            globals()['At%.2i' %noise].config(state = DISABLED)

    if isAttr:
        AS.config(state = NORMAL)
        AX.config(state = NORMAL)
    else :
        AS.config(state = DISABLED)
        AX.config(state = DISABLED)

    p['p1'].set(lp1[0])
    p['p2'].set(lp2[0])


def paramAndDir(o):                 # (0:11)
    global dir_, p1lst, p2lst, dir3

    dir_ = dirPriPlusSub()
    if o in [0,2,8,9]:
        p1lst = [p['p1'].get()]
        p2lst = [p['p2'].get()]
    elif o in [1,5,6,7,10,11,12,13]:
        p1lst = lp1.copy()
        p2lst = lp2.copy()
    elif o == 3:
        p1lst = lp1.copy()
        p2lst = [p['p2'].get()]
    elif o == 4:
        p1lst = [p['p1'].get()]
        p2lst = lp2.copy()
    dir3 = dirTri(p1lst, p2lst)
    

def loadPatterns(adress):           # (0,1)
    try:
        try:
            nruter = Te.data2array(adress)
        except:
            nruter = Te.data2array(adress + '/patterns.npy')
    except:
        S = []
        if os.path.exists(adress):
            nb_patt = Te.numberOfPatterns(adress)
            for i in xrange(nb_patt):
                S.append( np_load(adress + '/pattern_%.3i.npy' %i) )
        else:
            print "%s doesn't exists" %adress
        nruter = array(S)
    return nruter

def groupSimilarity(E,F):           # (1)
    sim = []
    for i in xrange(len(E)):
        for j in xrange(len(F)):
            sim.append(corrcoef(E[i], F[j])[0,1])
    return sim

def clustering():                   # (1)
    #TODO updgrade and check loadpatterns (only one dim now)
    S = loadPatterns(dir3[0][0])

    sim_final = 1
    while len(S) > 1 and sim_final > p['simcoef'].get():
        sim = -1
        for i in xrange(len(S)):
            for j in xrange(i + 1, len(S)):
                simEF = miea(groupSimilarity(S[i],S[j]))
                if simEF > sim:
                    indE = i
                    indF = j
                    sim = simEF
        sim_final = sim
        for i in xrange(len(S[indF])):
            S[indE].append(S[indF].pop(0))
        S.remove([])
        print len(S)
    return S

def miniClust(S, miea):
    ''' Return arguments for separated and clusterized patterns
    miea argument is min, max or mean
    '''

    correlats = similarity(S,S)
    correlats-= diag(diag(correlats))
    clusters  = [[k] for k in range(len(S))]

    sim_final = 1
    while len(clusters) > 1 and sim_final > p['simcoef'].get():
        sim = -1
        for i in xrange(len(clusters)):
            for j in xrange(i + 1, len(clusters)):
                simEF = miea(correlats[clusters[i]][:,clusters[j]])
                if simEF > sim:
                    indE = i
                    indF = j
                    sim = simEF
        sim_final = sim
        for i in xrange(len(clusters[indF])):
            clusters[indE].append(clusters[indF].pop(0))
        clusters.remove([])

    return clusters

def tendances():                    # (5)
    global G_tend
    G_tend = zeros((len(p1lst), len(p2lst)))
    for i in xrange(len(p1lst)):
        for j in xrange(len(p2lst)):
            params = dir3[i][j]
            if os.path.exists(params):
                tendance = Te.loadTendance(params)
                if len(tendance[0]) == 1:
                    G_tend[i,j] = 1
                else:
                    G_tend[i,j] = array(tendance[:,1:]).sum(0).max() / 3300.
            else:
                print "%s doesn't exists" %params
            print 'Norm',i,'; Omega',j
    G_tend[-1,-1] = 0
    G_tend[0,0] = 1
    print G_tend


def quantities():                   # (6)
    global G_patt
    try:
        G_patt = np_load(dir_ + 'G_Patt.npy')
    except:
        G_patt = zeros((len(p1lst), len(p2lst)))
        for i, p1 in zip(range(len(p1lst)), p1lst):
            for j, p2 in zip(range(len(p2lst)), p2lst):
                
                filename = adressFromPar(p1=p1, p2=p2, pre=dir_+'patterns_', post='.npy')
                
                if Te.isfile(filename):
                    G_patt[i,j] = Te.numberOfPatterns(filename)
                else:
                    print "%s doesn't exists" %filename
                print 'Norm',i,'; Omega',j
        np_save(dir_ + 'G_Patt.npy', G_patt)
        print G_patt.tolist()
    print G_patt.max()


def correlation():                   # (14)
    global G_corr
    try:
        G_corr = np_load(dir_ + 'G_Corr.npy')
    except:
        G_corr = zeros((len(p1lst), len(p2lst)))
        FCE = Te.data2array(Te.Pdir('Connectomes') + '/Hagmann/FC_D_998_0.npy')
        
        for i, p1 in zip(range(len(p1lst)), p1lst):
            for j, p2 in zip(range(len(p2lst)), p2lst):
                
                filename = adressFromPar(p1=p1, p2=p2, pre=dir_+'patterns_', post='.npy')
                if Te.isfile(filename):
                    patts = Te.data2array(filename)
                    FCP = Tf.fPearsonCorrelation(patts, finite=True)
                    Te.array2data(FCP, adressFromPar(p1=p1, p2=p2, pre=dir_+'FC/FC_', post='.npy'))
                    G_corr[i,j] = Tm.triSupPearsonCor(FCP, FCE)
                else:
                    print "%s doesn't exists" %filename
                print 'Norm',i,'; Omega',j
        np_save(dir_ + 'G_Corr.npy', G_corr)
        print G_corr.tolist()
    print G_corr.max()


def clusterizedQuant():             # (12)
    global G_CluPatt
    filename = dir_ + 'G_CluPatt_%.2f_s%i.npy' %(p['simcoef'].get(), p['similType'].get())
    try:
        G_CluPatt = np_load(filename)
    except:
        G_CluPatt = zeros((len(p1lst), len(p2lst)))
        for i in xrange(len(p1lst)):
            for j in xrange(len(p2lst)):
                params = dir3[i][j]
                if os.path.exists(params):
                    patts = loadPatterns(params)
                    clusters = miniClust(patts, miea)
                    G_CluPatt[i,j] = len(clusters)
                else:
                    print "%s doesn't exists" %params
                print 'Norm',i,'; Omega',j
        np_save(filename, G_CluPatt)
        print G_CluPatt.tolist()
    print G_CluPatt.max()


def densities():                    # (7)
    global G_dens
    try:
        G_dens = np_load(dir_ + 'G_dens.npy')        
    except:
        G_dens = zeros((len(p1lst), len(p2lst)))

        for i, p1 in zip(range(len(p1lst)), p1lst):
            for j, p2 in zip(range(len(p2lst)), p2lst):
                
                filename = adressFromPar(p1=p1, p2=p2, pre=dir_+'patterns_', post='.npy')
               
                if Te.isfile(filename):
                    patties = loadPatterns(filename)
                    dens = patties.sum(1) / (1.* patties.shape[1])
    
                    try:
                        G_dens[i,j] = sum(dens) / (1.*len(dens))
                    except:
                        G_dens[i,j] = 0
                else:
                    print "%s doesn't exists"
                print 'Norm',i,'; Omega',j
        np_save(dir_ + 'G_dens.npy', G_dens)
        print G_dens.tolist()

def occurrencesHistograms():        # (10)
    global G_Ohist
    try:
        G_Ohist = np_load(dir_ + 'G_Ohist.npy')
    except:
        G_Ohist = []
        for i in xrange(len(p1lst)):
            G_Ohist.append([])
            for j in xrange(len(p2lst)):
                params = dir3[i][j]
                if os.path.exists(params):
                    tendance = Te.loadTendance(params)
                    if len(tendance[0]) == 1:
                        G_Ohist[i].append(array([1]))
                    else:
                        tendance_resh = []
                        for k in xrange(1, len(tendance[0])):
                            tendance_resh.append(array(tendance[:,k]).sum(0) / 3300.)
                        tendance_resh = array(tendance_resh)
                        tendance_wout = where(tendance_resh < 0.05, -1, tendance_resh)
                        G_Ohist[i].append(tendance_wout)
                else:
                    print "%s doesn't exists" %params
                print 'Norm',i,'; Omega',j
        np_save(dir_ + 'G_Ohist.npy', G_Ohist)

def densitiesHistograms():          # (11)
    global G_Dhist
    try:
        G_Dhist = np_load(dir_ + 'G_Dhist.npy')
    except:
        G_Dhist = []
        for i in xrange(len(p1lst)):
            G_Dhist.append([])
            for j in xrange(len(p2lst)):
                params = dir3[i][j]
                if os.path.exists(params):
                    nb_patt = Te.numberOfPatterns(params)
                    dens = []
                    for k in xrange(nb_patt):
                        patty = np_load(params + '/pattern_%.3i.npy' %k)
                        dens.append(sum(patty) / (1.*len(patty)))
                    if len(dens) == 0:
                        G_Dhist[i].append(array([0]))
                    else:
                        G_Dhist[i].append(array(dens))
                else:
                    print "%s doesn't exists" %params
                print 'Norm',i,'; Omega',j
        np_save(dir_ + 'G_Dhist.npy', G_Dhist)

def showPatterns():                 # 0     Patterns
    filename = adressFromPar(p1=p1lst[0], p2=p2lst[0], pre=dir_+'patterns_', post='.npy')
    S = loadPatterns(filename)

    print len(S)
    clusL,clusP = [],[]
    clusP.append(0)
    for i in xrange(len(S)):
        if o == 1:
            clusL.append(len(S[i]))
            clusP.append(sum([clusP[-1],len(S[i])]))
    clusP.pop()

    p['F'][o].clear()
    #try:
        #frequency = Te.loadTendance(dir3[0][0].replace('X', 'S'))[:,1:].sum(0)
        #frequency,S = sort2by1(frequency,S, inverse=True)
        #if p['hide'].get():
            #frequency,S = frequency[:p['hide'].get()],S[:p['hide'].get()]

        #gs = gridspec.GridSpec(2, 1, height_ratios=[1, 7], hspace=0.)
        #p['A'][o] = p['F'][o].add_subplot(gs[0])
        #p['A'][o].bar(arange(len(frequency)), frequency, width=1)
        #p['A'][o].set_xticks([])
        #p['A'][o].set_ylim((0,frequency.max()))
        #p['A'][o].set_yticks((0,frequency.max()))
        #p['A'][o] = p['F'][o].add_subplot(gs[1])
    #except:
    if p['hide'].get():
        S = S[:p['hide'].get()]
        
    p['A'][o] = p['F'][o].add_subplot(111)

    if not p['switch'].get():
        if p['max'].get():
            p['A'][o].imshow(S, interpolation="nearest", \
                                aspect='auto', cmap=get_cmap(p['cmap'].get()+p['c_rev'].get()*'_r'), \
                                vmin=0, vmax=p['max'].get())
        else:
            p['A'][o].imshow(S, interpolation="nearest", aspect='auto', \
                                cmap=get_cmap(p['cmap'].get()+p['c_rev'].get()*'_r'))
        p['A'][o].xaxis.set_visible(False)
        if o ==1:
            p['A'][o].yaxis.set_visible(True)
            p['A'][o].set_yticks(clusP)
            p['A'][o].set_yticklabels(clusL)
        else:
            p['A'][o].yaxis.set_visible(False)
    else:
        if p['max'].get():
            p['A'][o].imshow(S.T, interpolation="nearest",aspect='auto', \
                                  cmap=get_cmap(p['cmap'].get()+p['c_rev'].get()*'_r'), \
                                  vmin=0, vmax=p['max'].get())
        else:
            p['A'][o].imshow(S.T, interpolation="nearest",aspect='auto',  \
                                  cmap=get_cmap(p['cmap'].get()+p['c_rev'].get()*'_r'))
        if   len(S[0]) == 66 : p['A'][o].set_yticks(xrange(66))
        elif len(S[0]) == 998: p['A'][o].set_yticks(a998i)
        p['A'][o].set_yticklabels(abre_ind, fontsize=9)
        if o ==1:
            p['A'][o].xaxis.set_visible(True)
            p['A'][o].set_xticks(clusP)
            p['A'][o].set_xticklabels(clusL)
        else:
            p['A'][o].xaxis.set_visible(False)
        p['A'][o].yaxis.set_visible(True)
        p['A'][o].tick_params(labelright=True)
    p['C'][o].show()
    if p['max'].get():
        permitSwitches(S, cm_x, maxs=(0, p['max'].get()))
    else:
        permitSwitches(S, cm_x)

def showCorrel():                 # 0     FCE / FCP Correlation
    global G_corr
    
    dFC = dir_+'FC/'
    if not os.path.exists(dFC):
        os.mkdir(dFC)
        
    correlation()
    G = G_corr

    if p['norm'].get():
        G /= G.max(0)

    p['F'][o].clear()
    p['A'][o] = p['F'][o].add_subplot(111)
    recurrentAxes()
    if p['seuil'].get():
        b = p['A'][o].contourf(p1lst, p2lst, where(G.T < p['seuil'].get(), G.T, p['seuil'].get()), CNTRF, \
                                cmap=get_cmap(p['cmap'].get()+p['c_rev'].get()*'_r'), \
                                antialiased=False)
    else:
        b = p['A'][o].contourf(p1lst, p2lst, G.T, CNTRF,  \
                                cmap=get_cmap(p['cmap'].get()+p['c_rev'].get()*'_r'), \
                                antialiased=False)
    p['F'][o].colorbar(b)

    #if not G_patt == None:
        #b = p['A'][o].contour(p1lst, p2lst, G_patt.T, p['Ncont'].get(), colors=p['cclbl'].get(), fontsize=15)
        #p['A'][o].clabel(b, inline=True, fmt='%i', fontsize=15, rightside_up=False)

    p['C'][o].show()

def showClusters():                 # 1     Clustering
    #TODO check why the first time it work and not after
    S0 = loadPatterns(dir3[0][0])
    S = []
    clusters = miniClust(S0, miea)
    for i in range(len(clusters)):
        S.extend(S0[clusters[i]].tolist())
    S = array(S)

    print len(S)


    p['F'][o].clear()
    p['A'][o] = p['F'][o].add_subplot(111)

    if not p['switch'].get():
        if p['max'].get():
            p['A'][o].imshow(S, interpolation="nearest", \
                                aspect='auto', cmap=get_cmap(p['cmap'].get()+p['c_rev'].get()*'_r'), \
                                vmin=0, vmax=p['max'].get())
        else:
            p['A'][o].imshow(S, interpolation="nearest", aspect='auto', \
                                cmap=get_cmap(p['cmap'].get()+p['c_rev'].get()*'_r'))
        p['A'][o].xaxis.set_visible(False)
        if o ==1:
            p['A'][o].yaxis.set_visible(True)
        else:
            p['A'][o].yaxis.set_visible(False)
    else:
        if p['max'].get():
            p['A'][o].imshow(S.T, interpolation="nearest",aspect='auto', \
                                  cmap=get_cmap(p['cmap'].get()+p['c_rev'].get()*'_r'), \
                                  vmin=0, vmax=p['max'].get())
        else:
            p['A'][o].imshow(S.T, interpolation="nearest",aspect='auto',  \
                                  cmap=get_cmap(p['cmap'].get()+p['c_rev'].get()*'_r'))
        if   len(S[0]) == 66 : p['A'][o].set_yticks(xrange(66))
        elif len(S[0]) == 998: p['A'][o].set_yticks(a998i)
        p['A'][o].set_yticklabels(abre_ind, fontsize=9)
        if o ==1:
            p['A'][o].xaxis.set_visible(True)
        else:
            p['A'][o].xaxis.set_visible(False)
        p['A'][o].yaxis.set_visible(True)
        p['A'][o].tick_params(labelright=True)
    p['C'][o].show()
    permitSwitches(S, cm_x)

def showOccurrences0D():            # 2     Occurrences 0D
    ''' Affiche les tendances des patterns pour les paramètres sélectionnés sous forme d'un bar plot
    avec des rectangles correspondants aux tendances des patterns en partant du bas.
    '''
    tendance = Te.loadTendance(dir3[0][0].replace('X','S'))

    # Calcul l'ordre de taux d'apparition en fonction du parametre P
    freqTot = []
    for i in range(1,len(tendance[0])):
        freqTot.append(array(tendance[:,i]).sum())
    ordre = sort2by1(freqTot, arange(1,len(tendance[0])), inverse=True)[1]

    x     = tendance[:,0]
    p['A'][o].clear()
    p['A'][o].set_xlim(min(x),max(x))
    p['A'][o].set_ylim(0,1)
    p['A'][o].set_xlabel(r'$Initial\ density$', fontsize=18)
    p['A'][o].set_ylabel(r'$Attractors\ basin\ size$', fontsize=20)
    cusu0 = zeros(len(tendance))
    cusu1 = cusu0.copy()
    if p['hide'].get():
        notHide = int(min(len(tendance[0]), p['hide'].get() + 1))
    else:
        notHide = len(tendance[0])
    for i in range(1,notHide):
        cusu1 += array(tendance[:,ordre[i-1]]) / 100.
        #cusu1 += array(tendance[:,i]) / 100.
        p['A'][o].fill_between(x, cusu0, cusu1, color=rand(3), alpha=.6)
        cusu0 = cusu1.copy()
    p['C'][o].show()

def showOccurrences1D():            # 3/4   Predominances 1D
    # Extraction des tendances
    frequencies = []

    if len(p1lst) > len(p2lst):
        paramX = p1lst
    else:
        paramX = p2lst

    for i in xrange(len(paramX)):
        if len(p1lst) > len(p2lst):
            params = dir3[i][0].replace('X','S')
        else:
            params = dir3[0][i].replace('X','S')
        if os.path.exists(params):
            tendance = Te.loadTendance(params)
            frequencies.append( (array(tendance[:,1:]).sum(0) / 3300.).tolist() )
        else:
            print "%s doesn't exists" %params

    # Remplissage par des zeros afin d'obtenir une matrice
    Nmax = 0
    Plen = len(paramX)
    for i in range(Plen):
        Nmax = max(Nmax, len(frequencies[i]))
    for i in range(Plen):
        for j in range(len(frequencies[i]), Nmax):
            frequencies[i].append(0)
    frequencies = array(frequencies)

    # Création de la bibliothèque des patterns selon leur apparition
    predoName = 'Predominances_%.2f_' %p['simcoef'].get()
    if p['similType'].get():
        predoName += '_sim'
    if len(p1lst) > len(p2lst):
        dicoFileName = dir_ + predoName + '%s.eva' %p['p2'].get()
    else:
        predoName += '_T'
        dicoFileName = dir_ + predoName + '%s.eva' %p['p1'].get()
    if os.path.exists(dicoFileName):
        dicoFile = open(dicoFileName, "rb")
        dico = pk_load(dicoFile)
    else:
        t1   = time()
        dico = {}
        # dico is 3D: number of futur stacked bars x appearing depending on P/G x (parameter P/G, pattern i, frequency)

        for q in xrange(len(paramX)):
            t0=time()
            ''' Purposes:
            params →
            S0/1    →
            clu0/1  →
            pos0/1  →
            arg0/1  →
            '''
            if len(p1lst) > len(p2lst):
                params = dir3[q][0]
            else:
                params = dir3[0][q]
            if q == 0:
                S1 = loadPatterns(params)
                clu1 = miniClust(S1, miea)
                pos1 = arange(len(clu1))
                for pos in pos1:
                    dico[pos] = [],[],[]
                    dico[pos][0].append(q)
                    dico[pos][1].append(frequencies[q,clu1[pos]].sum())
            else:
                S0   = array(S1)
                clu0 = list(clu1)
                del S1,clu1
                S1   = loadPatterns(params)
                clu1 = miniClust(S1, miea)

                arg0 = range(len(clu0))
                arg1 = range(len(clu1))
                pos0 = pos1.copy()
                pos1 = -ones(len(clu1), dtype=int)
                while arg0 and arg1:
                    try:    del corr
                    except: pass
                    corr = zeros((len(arg0), len(arg1)))
                    for i0,a0 in enumerate(arg0):
                        for i1,a1 in enumerate(arg1):
                            corr[i0,i1] = miea(fProximity(S0[clu0[a0]], S1[clu1[a1]]))
                        #if p['similType'].get():
                            #corr[i0] = similAndNotCorr(S0[a0], S1[arg1])
                        #else:
                            #corr[i0] = fastPearsonCorrelation1D(S0[a0], S1[arg1])

                    # To consider if it's a new pattern or not
                    if corr.max() > 0.8: # there is a same pattern
                        indM = unravel_index(corr.argmax(), corr.shape)
                        pos  = pos1[arg1[indM[1]]] = pos0[arg0[indM[0]]]
                        dico[pos][0].append(q)
                        dico[pos][1].append(frequencies[q, clu1[arg1[indM[1]]]].sum())
                        arg0.pop(indM[0])
                        arg1.pop(indM[1])
                    else: # New pattern
                        break
                    print len(arg0),'/',len(S0)
                if arg1: # Add news patterns in dico
                    news = arange(len(arg1))
                    lenD = len(dico)
                    for pos in news:
                        dico[lenD + pos] = [],[],[]
                        dico[lenD + pos][0].append(q)
                        dico[lenD + pos][1].append(frequencies[q, clu1[arg1[pos]]].sum())
                    pos1[pos1 == -1] = lenD + news
            print 'Etape %i/%i' %(q+1,len(paramX)),time()-t0
        print 'Temps de calcul:', time()-t1
        dicoFile = open(dicoFileName, "wb")
        pk_dump(dico, dicoFile)

        #t1   = time()

        #''' dico is 3D: number of futur stacked bars x appearing depending on P/G x (parameter P/G, pattern i, frequency)'''
        #dico = {}

        #for q in xrange(len(paramX)):
            #t0=time()
            #if len(p1lst) > len(p2lst):
                #params = dir3[q][0]
            #else:
                #params = dir3[0][q]
            #if q == 0:
                #S1 = loadPatterns(params)
                #pos0 = range(len(S1))
                #for pos in pos0:
                    #dico[pos] = [],[],[]
                    #dico[pos][0].append(q)
                    #dico[pos][1].append(pos)
                    #dico[pos][2].append(frequencies[q,pos])
            #else:
                #S0 = array(S1)
                #del S1
                #S1 = loadPatterns(params)

                #arg0 = range(len(S0))
                #arg1 = range(len(S1))
                #pos1 = -ones(len(S1), dtype=int)
                #while arg0 and arg1:
                    #try:    del corr
                    #except: pass
                    #corr = zeros((len(arg0), len(arg1)))
                    #for i0,a0 in enumerate(arg0):
                        #if p['similType'].get():
                            #corr[i0] = similAndNotCorr(S0[a0], S1[arg1])
                        #else:
                            #corr[i0] = fastPearsonCorrelation1D(S0[a0], S1[arg1])

                    ## To consider if it's a new pattern or not
                    #if corr.max() > 0.8: # there is a same pattern
                        #indM = unravel_index(corr.argmax(), corr.shape)
                        #pos  = pos1[arg1[indM[1]]] = pos0[arg0[indM[0]]]
                        #dico[pos][0].append(q)
                        #dico[pos][1].append(arg1[indM[1]])
                        #dico[pos][2].append(frequencies[q, arg1[indM[1]]])
                        #arg0.pop(indM[0])
                        #arg1.pop(indM[1])
                    #else: # New pattern
                        #break
                    #print len(arg0),'/',len(S0)
                #if arg1: # Add news patterns in dico
                    #news = arange(len(arg1))
                    #lenD = len(dico)
                    #for pos in news:
                        #dico[lenD + pos] = [],[],[]
                        #dico[lenD + pos][0].append(q)
                        #dico[lenD + pos][1].append(arg1[pos])
                        #dico[lenD + pos][2].append(frequencies[q, arg1[pos]])
                    #pos1[pos1 == -1] = lenD + news
                #pos0 = pos1
            #print 'Etape %i/%i' %(q,len(paramX)),time()-t0
        #print 'Temps de calcul:', time()-t1
        #dicoFile = open(dicoFileName, "wb")
        #pk_dump(dico, dicoFile)
    print 'Longueur du dico:', len(dico)
    dicoFile.close()

    # Calcul l'ordre de taux d'apparition en fonction du parametre P
    freqTot = []
    for i in range(len(dico)):
        freqTot.append(array(dico[i][1]).sum())
    ordre = sort2by1(freqTot, arange(len(dico)), inverse=True)[1]

    # Affiche les patterns selon leur ordre inverse de taux d'apparition
    p['A'][o].clear()
    if p['POG_PG'].get():
        if len(p1lst) < len(p2lst):
            paramX = array(paramX) * 2.
            p['A'][o].set_xlabel(r'$G$', fontsize=20)
        else:
            paramX = array(paramX) / 2.
            p['A'][o].set_xlabel(r'$P$', fontsize=20)
    else:
        if len(p1lst) < len(p2lst):
            p['A'][o].set_xlabel(r'$%s$' %p2str, fontsize=20)
        else:
            p['A'][o].set_xlabel(r'$%s$' %p1str, fontsize=20)
    p['A'][o].set_xlim(min(paramX),max(paramX))
    p['A'][o].set_ylim(0,1)
    p['A'][o].set_ylabel(r'$Attractors\ basin\ size$', fontsize=20)
    cusu0 = zeros(Plen)
    cusu1 = cusu0.copy()
    facHide = ones(Plen)
    if p['hide'].get():
        notHide = int(min(len(dico), p['hide'].get()))
        for i in range(notHide):
            cusu1[dico[ordre[i]][0]] += array(dico[ordre[i]][1])
            cusu0 = cusu1.copy()
        facHide = cusu1
        cusu0 = zeros(Plen)
        cusu1 = cusu0.copy()
    else:
        notHide = len(dico)
    for i in range(notHide):
        cusu1[dico[ordre[i]][0]] += array(dico[ordre[i]][1])
        p['A'][o].fill_between(paramX, cusu0/facHide, cusu1/facHide, color=rand(3), alpha=.6)
        cusu0 = cusu1.copy()
    p['C'][o].show()

def showOccurrences2D():            # 5     Occurrences 2D
    global G_tend
    if G_tend == None:
        tendances()

    p['F'][o].clear()
    p['A'][o] = p['F'][o].add_subplot(111)
    recurrentAxes()
    #if len(p['F'][o].axes) > 1:  # delete old colorbar
        #p['F'][o].delaxes(p['F'][o].axes[1])
        #p['F'][o].subplots_adjust(right=0.90)

    #p['A'][o].set_xlabel(r'$\Psi$', fontsize=20)
    #p['A'][o].set_ylabel(r'$\Omega$', fontsize=20)
    b = p['A'][o].contourf(p1lst, p2lst, G_tend.T, CNTRF,  \
                                cmap=get_cmap(p['cmap'].get()+p['c_rev'].get()*'_r'), \
                                antialiased=False)
    p['F'][o].colorbar(b)
    p['C'][o].show()

def showQuantities(clust=False):    # 6/12  Quantities (Clusterized or not)
    global G_patt
    if not clust:
        quantities()
        G = G_patt
    else:
        clusterizedQuant()
        G = G_CluPatt

    if p['norm'].get():
        G /= G.max(0)

    p['F'][o].clear()
    p['A'][o] = p['F'][o].add_subplot(111)
    recurrentAxes()
    if p['seuil'].get():
        b = p['A'][o].contourf(p1lst, p2lst, where(G.T < p['seuil'].get(), G.T, p['seuil'].get()), CNTRF, \
                                cmap=get_cmap(p['cmap'].get()+p['c_rev'].get()*'_r'), \
                                antialiased=False)
    else:
        b = p['A'][o].contourf(p1lst, p2lst, G.T, CNTRF,  \
                                cmap=get_cmap(p['cmap'].get()+p['c_rev'].get()*'_r'), \
                                antialiased=False)
    p['F'][o].colorbar(b, format='%i')

    if not G_dens == None:
        b = p['A'][o].contour(p1lst, p2lst, G_dens.T, p['Ncont'].get(), colors=p['cclbl'].get(), fontsize=15)
        p['A'][o].clabel(b, inline=True, fmt='%.2f', fontsize=15, rightside_up=False)

    p['C'][o].show()

def showDensities():                # 7     Densities
    global G_dens
    densities()

    p['F'][o].clear()
    p['A'][o] = p['F'][o].add_subplot(111)
    recurrentAxes()
    b = p['A'][o].contourf(p1lst, p2lst, G_dens.T, CNTRF,  \
                                cmap=get_cmap(p['cmap'].get()+p['c_rev'].get()*'_r'), \
                                antialiased=False)
    p['F'][o].colorbar(b)
    p['C'][o].show()

def showDensHistogram():            # 8     Densities Histogram
    #TODO Add polar axes for polar histograms
    params = dir3[0][0]
    nb_patt = Te.numberOfPatterns(params)
    dens = zeros(nb_patt)
    for k in xrange(nb_patt):
        patty   = np_load(params + '/pattern_%.3i.npy' %k)
        dens[k] = sum(patty) / (1.*len(patty))


    p['F'][o].clear()
    p['A'][o] = p['F'][o].add_subplot(211, polar=False)
    p['A'][o].hist(dens, bins=linspace(0, 1, 42))
    p['A'][o].set_xlabel(r'$Densities$', fontsize=15)
    #p['A'][o].set_xscale('log')
    #p['A'][o].set_yscale('log')

    p['A'][o] = p['F'][o].add_subplot(212, polar=True)
    p['A'][o].set_theta_zero_location("S")
    p['A'][o].set_theta_direction(-1)
    p['A'][o].set_axis_off()
    #p['A'][o].set_xscale('log')
    #p['A'][o].set_yscale('log')
    p['A'][o].hist(dens *2.*pi, bins=42, rwidth=1, range=(0, 2*pi))
    p['C'][o].show()

def showHistogram():                # 9     Histograms ?
    #TODO checker tendance_wout
    ''' Affiche l'histogramme des tendances des patterns pour les paramètres sélectionnés
    '''
    p['F'][o].clear()

    tendance = Te.loadTendance(dir3[0][0])
    tendance_resh = []
    for i in xrange(1, len(tendance[0])):
        tendance_resh.append(array(tendance[:,i]).sum(0) / 3300.)
    tendance_resh.sort(reverse=True)
    tendance_resh = array(tendance_resh)
    Kst = tendance_resh * arange(1,len(tendance_resh)+1)
    tendance_wout = where(tendance_resh < 0.0, -1., tendance_resh)
    tendance_wout2 = where(tendance_resh < 0.05, -1., tendance_resh)

    if size(tendance) == 1:
        return 0
    p['A'][o] = p['F'][o].add_subplot(311, polar=False)
    p['A'][o].hist(tendance_resh, bins=linspace(0, 1, 42))
    p['A'][o].set_ylabel(r'$Distribution$', fontsize=15)

    p['A'][o] = p['F'][o].add_subplot(312, polar=False)
    p['A'][o].hist(tendance_resh, bins=200, range=(0, 1))
    p['A'][o].set_ylabel(r'$Log\ Distribution$', fontsize=15)
    p['A'][o].set_xscale('log')
    p['A'][o].set_yscale('log')

    p['A'][o] = p['F'][o].add_subplot(313, polar=False)
    p['A'][o].plot(Kst,'-+')
    p['A'][o].set_ylabel(r'$K$', fontsize=15)

    p['A'][o] = p['F'][o].add_subplot(322, polar=True)
    p['A'][o].set_theta_zero_location("S")
    p['A'][o].set_theta_direction(-1)
    p['A'][o].set_axis_off()
    #p['A'][o].set_xscale('log')
    #p['A'][o].set_yscale('log')
    p['A'][o].hist(tendance_wout2 *2.*pi, bins=42, rwidth=1, range=(0, 2*pi))
    p['C'][o].show()

def showOccHistograms():            # 10    Occurrences Distributions
    #TODO find a faster way to subplot
    global G_Ohist

    if G_Ohist == None:
        occurrencesHistograms()
        recurrentAxes()
        if not G_patt == None:
            p['A'][o].contourf(p1lst, p2lst, G_patt.T, CNTRF,  \
                                aspect='auto', cmap=get_cmap(p['cmap'].get()+p['c_rev'].get()*'_r'), \
                                antialiased=False, alpha=0.5)

        for i in xrange(len(p1lst)):
            for j in xrange(len(p2lst)):
                a = p['F'][o].add_subplot(len(p1lst), len(p2lst), i+len(p1lst)*(len(p2lst)-j-1)+1, polar=True)
                a.set_theta_zero_location("S")
                a.set_theta_direction(-1)
                a.set_axis_off()
                a.hist(G_Ohist[i][j] *2.*pi, bins=20, rwidth=1, range=(0, 2*pi))
                #N, bins, patches = a.hist(G_Ohist[i][j] *2.*pi, bins=20, rwidth=1, range=(0, 2*pi))
                #for thisfrac, thispatch in zip(linspace(1,0,len(N)), patches):
                    #color = cm.get_cmap(p['cmap'].get())(thisfrac)
                    #thispatch.set_facecolor(color)
        p['C'][o].show()

def showDensHistograms():           # 11    Densities Distributions
    #TODO find a faster way to subplot
    global G_Dhist

    if G_Dhist == None:
        densitiesHistograms()
        recurrentAxes()
        if not G_patt == None:
            p['A'][o].contourf(p1lst, p2lst, G_patt.T, CNTRF,  \
                                aspect='auto', cmap=get_cmap(p['cmap'].get()+p['c_rev'].get()*'_r'), \
                                antialiased=False, alpha=0.5)

        l1 = len(p1lst)
        l2 = len(p2lst)
        for i in xrange(l1):
            for j in xrange(l2):
                a = p['F'][o].add_subplot(l2, l1, i+l1*(l2-j-1)+1, polar=True)
                a.set_theta_zero_location("S")
                a.set_theta_direction(-1)
                a.set_axis_off()
                a.hist(G_Dhist[i][j] *2.*pi, bins=42, rwidth=1, range=(0, 2*pi))
                #N, bins, patches = a.hist(G_Dhist[i][j] *2.*pi, bins=20, rwidth=1, range=(0, 2*pi))
                #for thisfrac, thispatch in zip(linspace(1,0,len(N)), patches):
                    #color = cm.get_cmap(p['cmap'].get())(thisfrac)
                    #thispatch.set_facecolor(color)
        p['C'][o].show()

def showAvgGOF():                   # 13    Average Goodness of Fit
    global p1lst, p2lst
    noise = [00,10]
    for ii in range(5):
        for nn in range(len(noise)):
            grid2[ii+nn*5].clear()

    if p['POG_PG'].get():
        p1lst /= 2.
        p2lst *= 2.
        grid2.axes_llc.set_xlabel(r'$P$', fontsize=20)
        grid2.axes_llc.set_ylabel(r'$G$', fontsize=20)
    else:
        grid2.axes_llc.set_xlabel(r'$%s$' %p1str, fontsize=20)
        grid2.axes_llc.set_ylabel(r'$%s$' %p2str, fontsize=20)

    for ii in range(5):
        for nn in range(len(noise)):
            avgGOF = np_load(p['dir_pri'].get() + '/At%.2i_S/avgGOF_p%i_fD1_c%i.npy' \
                                                   %(noise[nn],ii+1,not p['similType'].get()))
            im = grid2[ii+nn*5].contourf(p1lst, p2lst, avgGOF.T, CNTRF, \
                                cmap=get_cmap(p['cmap'].get()+p['c_rev'].get()*'_r'), \
                                antialiased=False)

    grid2.cbar_axes[0].colorbar(im)
    grid2.set_label_mode('L')
    p['C'][o].draw()

def recurrentAxes():
    global p1lst, p2lst
    if p['POG_PG'].get():
        p1lst /= 2.
        p2lst *= 2.
        p['A'][o].set_xlabel(r'$P$', fontsize=20)
        p['A'][o].set_ylabel(r'$G$', fontsize=20)
    #else:
        #p['A'][o].set_xlabel(r'$\Psi$', fontsize=20)
        #p['A'][o].set_ylabel(r'$\Omega$', fontsize=20)
    else:
        p['A'][o].set_xlabel(r'$%s$' %p1str, fontsize=20)
        p['A'][o].set_ylabel(r'$%s$' %p2str, fontsize=20)

def Show():
    global o, similarity, miea
    if   p['similType'].get() == 0: similarity = fCosineSimilarity
    elif p['similType'].get() == 1: similarity = fProximity
    elif p['similType'].get() == 2: similarity = fOverDistance
    if   p['mmm'].get() == 0: miea = type(rand(int())).min
    elif p['mmm'].get() == 1: miea = type(rand(int())).mean
    elif p['mmm'].get() == 2: miea = type(rand(int())).max

    o = p['PatTen'].get()
    mapCanvas(o)
    paramAndDir(o)
    if   o ==  0: showPatterns()                # Patterns
    elif o ==  1: showClusters()                # Clustering
    elif o ==  2: showOccurrences0D()           # Occ. 0D
    elif o ==  3: showOccurrences1D()           # Occ. 1D Psi
    elif o ==  4: showOccurrences1D()           # Occ. 1D Psi
    elif o ==  5: showOccurrences2D()           # Occ .2D
    elif o ==  6: showQuantities()              # Quantity
    elif o ==  7: showDensities()               # Densities
    elif o ==  8: showDensHistogram()           # Densities Hist
    elif o ==  9: showHistogram()               # Histogra ??
    elif o == 10: showOccHistograms()           # Occurrences Distributions
    elif o == 11: showDensHistograms()          # Densities Distributions
    elif o == 12: showQuantities(clust=True)    # Clusterized Quantities
    elif o == 13: showAvgGOF()                  # Average Goodness of Fit
    elif o == 14: showCorrel()                  # FCE / FCP correlation


'#################################################  PARAMETERS  #################################################'
global int_0, ext_0, int_N, ext_N, abre_ind, a998i
root = Tk()
root.title("Analyseur")
row_counter= 0
cm_x = cm.RdBu_r
D998_D66 = array([ \
     0,  0,  1,  1,  1,  1,  1,  1,  2,  2,  2,  3,  3,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,\
     4,  5,  5,  5,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,\
     7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,\
     8,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, 10, 10, 10, 10, 10,\
    10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,\
    12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 14, 14,\
    14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15,\
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 16, 16, 16,\
    16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 19,\
    19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,\
    21, 21, 21, 21, 21, 21, 21, 22, 22, 22, 22, 23, 23, 23, 23, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24,\
    24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 25, 25, 25, 25, 25, 25, 25,\
    25, 25, 25, 25, 25, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 28, 28,\
    28, 28, 28, 28, 28, 28, 28, 28, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 30, 30, 30, 30, 30, 30, 30, 30, 31, 31, 31, 31, 31, 31,\
    31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 32, 32, 32, 32, 32, 32, 32, 33, 33, 33, 33, 33, 33, 33, 34, 34, 34,\
    34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 35, 35, 35, 35, 35, 35, 35, 35, 36, 36, 36, 36, 36, 36,\
    36, 36, 36, 36, 36, 37, 37, 37, 37, 37, 37, 37, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39,\
    39, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41,\
    41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 42, 42,\
    42, 42, 43, 43, 43, 43, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 45, 45, 45, 45, 45, 45, 46, 46,\
    46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 47, 47, 47, 47, 47, 47, 47, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48,\
    48, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,\
    50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51,\
    51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52,\
    52, 52, 52, 52, 52, 52, 52, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 54, 54, 54, 54, 54, 55, 55, 55,\
    55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56,\
    56, 56, 56, 56, 56, 56, 56, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58,\
    58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59,\
    59, 59, 59, 59, 59, 60, 60, 60, 60, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 62, 62, 63,\
    63, 63, 63, 64, 64, 64, 64, 64, 64, 65, 65, 65])
AvgN_D66 = array([ \
    2,6,3,2,22,3,19,27,19,28,16,7,20,28,31,36,13,10,8,22,6,19,4,4,46,12,17,10,10,12,8,23,7,7,\
    23,8,11,8,9,16,12,50,4,4,20,6,19,7,11,13,36,30,29,19,5,19,25,17,27,22,4,22,2,4,6,3])
abrev = ['ENT','PARH','TP','FP','FUS','TT','LOCC','SP','IT','IP',\
         'SMAR','BTST','MT','ST','PSTC','PREC','CMF','POPE','PTRI','RMF',\
         'PORB','LOF','CAC','RAC','SF','MOF','LING','PCAL','CUN','PARC',\
         'ISTC','PCUN','PC']
abre_ind = Abreviations(abrev)
a998i = [1, 7, 10, 12, 34, 37, 56, 83, 102, 130, 146,\
         153, 173, 201, 232, 268, 281, 291, 299, 321, 327, 346,\
         350, 354, 400, 412, 429, 439, 449, 461, 469, 492, 499,\
         506, 529, 537, 548, 556, 565, 581, 593, 643, 647, 651,\
         671, 677, 696, 703, 714, 727, 763, 793, 822, 841, 846,\
         865, 890, 907, 934, 956, 960, 982, 984, 988, 994, 997]
region = Te.Pdir('Main/Regions/Images')
int_0 = imread(region + '/Int.png')      # Contours of internal part
ext_0 = imread(region + '/Ext.png')      # Contours of external part
int_N = zeros((33,513,837))              # Internal part
ext_N = zeros((33,513,837))              # External part
for i in range(33):
    try:    int_N[i] = imread(region + '/Int%s.png' %abrev[i])
    except: pass
    try:    ext_N[i] = imread(region + '/Ext%s.png' %abrev[i])
    except: pass

p = {}
p['dir_pri'] = StringVar();         p['dir_pri'].set("/home/golos/Main/Simulations/Patterns/SC0_998_DG")       # Dossier de simulation
p['dir_sub'] = StringVar();         p['dir_sub'].set("/Cand_S")                     # Sous-dossier pour les types de patterns
p['attr'] = StringVar();            p['attr'].set('00')                             # Valeur du bruit initial pour dir_sub = Attr_*
p['PatTen'] = IntVar();             p['PatTen'].set(0)                              # Type de Graphique

p['simcoef'] = DoubleVar();         p['simcoef'].set(0.75)                          # Coefficient de correlation de Clustering
p['rangeD'] = DoubleVar();          p['rangeD'].set(0.00)                           # Etendue du Clustering sur les dossiers
p['mmm'] = IntVar();                p['mmm'].set(0)                                 # Min / Moy / Max pour la similarite du Clustering

p['cmap'] = StringVar();            p['cmap'].set("Spectral")                       # Couleur des graphiques
p['max'] = DoubleVar();             p['max'].set(0)                                 # Maximum pour l'affichage des patterns
p['seuil'] = DoubleVar();           p['seuil'].set(0)                               # Seuil pour l'espace des patametres
p['hide'] = DoubleVar();            p['hide'].set(0)                                # Masque les patterns jusqu'à
p['Ncont'] = IntVar();              p['Ncont'].set(10)                              # Nombre de contour affiché
p['cclbl'] = StringVar();           p['cclbl'].set('k')                             # Color of the contour label
p['axis'] = IntVar();               p['axis'].set(0)                                # Affichage des abréviations
p['switch'] = IntVar();             p['switch'].set(1)                              # Patterns verticaux / horizontaux
p['c_rev'] = IntVar();              p['c_rev'].set(1)                               # Reverse colormap or not
p['POG_PG'] = IntVar();             p['POG_PG'].set(0)                              # Switch Psi/Omega/gamma to P/G
p['similType'] = IntVar();          p['similType'].set(0)                           # Similarity type (cosine, spe1, spe2)
p['norm'] = IntVar();               p['norm'].set(0)                                # Normalize the figure by the maximum on oen parameter

p['p1']    = DoubleVar();           p['p1'].set(1)                                  # Valeur du parametre 1
p['p2']    = DoubleVar();           p['p2'].set(1)                                  # Valeur du parametre 2


'##################################################  DISPLAY  ###################################################'
global At01,At02,At03,At04,At05,At06,At07,At08,At09,At10,At15,At20,At25,At30,At35,At40,At45,At50

W_frame   = LabelFrame(root, bd=0);                     W_frame.pack(side='left')
E_frame   = LabelFrame(root, bd=0);                     E_frame.pack(side='right', padx=5)
E1_frame  = LabelFrame(E_frame, text= "Clustering");    E1_frame.pack(side='top', fill=X)
E1a_frame = LabelFrame(E1_frame, bd=0);                 E1a_frame.pack(side='top', pady=4)
E1b_frame = LabelFrame(E1_frame, bd=0);                 E1b_frame.pack(side='top')
E1c_frame = LabelFrame(E1_frame, bd=0);                 E1c_frame.pack(side='top')
E2_frame  = LabelFrame(E_frame, text= "Directory");     E2_frame.pack(side='top', pady=4, fill=X)
E2a_frame = LabelFrame(E2_frame, bd=0);                 E2a_frame.pack(side='top')
E2b_frame = LabelFrame(E2_frame, bd=0);                 E2b_frame.pack(side='top')
E2c_frame = LabelFrame(E2_frame, bd=0);                 E2c_frame.pack(side='top')
E2d_frame = LabelFrame(E2_frame, bd=0);                 E2d_frame.pack(side='top')
E3_frame  = LabelFrame(E_frame, relief=SUNKEN);         E3_frame.pack(side='top', pady=8)
E4_frame  = LabelFrame(E_frame, text= "Displaying");    E4_frame.pack(side='top', pady=4, fill=X)
E4a_frame = LabelFrame(E4_frame, bd=0);                 E4a_frame.pack(side='top')
E4aL_frame= LabelFrame(E4a_frame, bd=0);                E4aL_frame.pack(side='left')
E4aR_frame= LabelFrame(E4a_frame, bd=0);                E4aR_frame.pack(side='right')
E4b_frame = LabelFrame(E4_frame, bd=0);                 E4b_frame.pack(side='top')
E4c_frame = LabelFrame(E4_frame, bd=0);                 E4c_frame.pack(side='top')


'##################################################  CANVAS  #################################################'
p['C'],p['F'],p['A'],p['T'] = [],[],[],[]
for i in range(16):
    p['F'].append(figure(figsize=(10.3, 9), dpi=100))
    #p['F'].append(figure(figsize=(10.3, 5), dpi=100))
    p['F'][i].patch.set_facecolor('white')
    p['A'].append(p['F'][i].add_subplot(111))
    p['C'].append(FigureCanvasTkAgg(p['F'][i], master=W_frame))
    p['C'][i].show()
    if i in [0,1,15]:
        p['C'][i].mpl_connect('motion_notify_event', getStatusXY)
    p['T'].append(NavigationToolbar2TkAgg(p['C'][i], W_frame))
    p['T'][i].pack_forget()
p['C'][0].get_tk_widget().pack()
p['T'][0].pack()
p['F'][-1].clf()
p['F'][13].clf()
grid = AxesGrid(p['F'][-1], 111,
                nrows_ncols = (2,2),
                axes_pad = 0.0,
                share_all = True,
                label_mode = 'L',
                cbar_location = 'right',
                cbar_mode = 'single')
grid2 = AxesGrid(p['F'][13], 111,
                nrows_ncols = (2,5),
                axes_pad = 0.0,
                aspect = False,
                share_all = True,
                cbar_location = 'right',
                cbar_mode = 'single')


'##################################################  WIDGETS  #################################################'
''' RIGHT 1 Clustering'''
textentry(E1a_frame, p['simcoef'], 'Sim. Coef', 7)
textentry(E1a_frame, p['rangeD'], 'Dir. Range', 7)
Radiobutton(E1b_frame, text="Min", variable=p['mmm'], value="0", command = Show, indicatoron=0).pack(side='left')
Radiobutton(E1b_frame, text="Moy", variable=p['mmm'], value="1", command = Show, indicatoron=0).pack(side='left')
Radiobutton(E1b_frame, text="Max", variable=p['mmm'], value="2", command = Show, indicatoron=0).pack(side='left')
Radiobutton(E1c_frame, text="Cosine", variable=p['similType'], value="0", command = Show, indicatoron=0).pack(side='left')
Radiobutton(E1c_frame, text="Spe 1",  variable=p['similType'], value="1", command = Show, indicatoron=0).pack(side='left')
Radiobutton(E1c_frame, text="Spe 2",  variable=p['similType'], value="2", command = Show, indicatoron=0).pack(side='left')

''' RIGHT 2 Directory'''
textentry(E2a_frame, p['dir_pri'], 'dir_pri', 20)
Button(E2b_frame, text='Refresh Lists', command=refreshLists).pack(side='top')
lista= Listbox(E2b_frame, width=10, exportselection=0);   lista.pack(side='left', pady=5)
listb= Listbox(E2b_frame, width=10, exportselection=0);   listb.pack(side='left', pady=5)
CS=Radiobutton(E2c_frame, text="Cand S", variable=p['dir_sub'], value="/Cand_S", command = Show, indicatoron=0); CS.grid(row=0, column=0)
CX=Radiobutton(E2c_frame, text="Cand X", variable=p['dir_sub'], value="/Cand_X", command = Show, indicatoron=0); CX.grid(row=1, column=0)
PS=Radiobutton(E2c_frame, text="Patt S", variable=p['dir_sub'], value="/Patt_S", command = Show, indicatoron=0); PS.grid(row=0, column=1)
PX=Radiobutton(E2c_frame, text="Patt X", variable=p['dir_sub'], value="/Patt_X", command = Show, indicatoron=0); PX.grid(row=1, column=1)
AS=Radiobutton(E2c_frame, text="Attr S", variable=p['dir_sub'], value="/Attr_S", command = Show, indicatoron=0); AS.grid(row=0, column=2)
AX=Radiobutton(E2c_frame, text="Attr X", variable=p['dir_sub'], value="/Attr_X", command = Show, indicatoron=0); AX.grid(row=1, column=2)
At01=Radiobutton(E2d_frame,text="00", variable=p['attr'], value="00", command = ShowAtr, indicatoron=0); At01.grid(row=0, column=0)
At01=Radiobutton(E2d_frame,text="01", variable=p['attr'], value="01", command = ShowAtr, indicatoron=0); At01.grid(row=0, column=1)
At02=Radiobutton(E2d_frame,text="02", variable=p['attr'], value="02", command = ShowAtr, indicatoron=0); At02.grid(row=0, column=2)
At03=Radiobutton(E2d_frame,text="03", variable=p['attr'], value="03", command = ShowAtr, indicatoron=0); At03.grid(row=0, column=3)
At04=Radiobutton(E2d_frame,text="04", variable=p['attr'], value="04", command = ShowAtr, indicatoron=0); At04.grid(row=0, column=4)
At05=Radiobutton(E2d_frame,text="05", variable=p['attr'], value="05", command = ShowAtr, indicatoron=0); At05.grid(row=0, column=5)
At06=Radiobutton(E2d_frame,text="06", variable=p['attr'], value="06", command = ShowAtr, indicatoron=0); At06.grid(row=0, column=6)
At07=Radiobutton(E2d_frame,text="07", variable=p['attr'], value="07", command = ShowAtr, indicatoron=0); At07.grid(row=0, column=7)
At08=Radiobutton(E2d_frame,text="08", variable=p['attr'], value="08", command = ShowAtr, indicatoron=0); At08.grid(row=0, column=8)
At09=Radiobutton(E2d_frame,text="09", variable=p['attr'], value="09", command = ShowAtr, indicatoron=0); At09.grid(row=0, column=9)
At10=Radiobutton(E2d_frame,text="10", variable=p['attr'], value="10", command = ShowAtr, indicatoron=0); At10.grid(row=1, column=1)
At15=Radiobutton(E2d_frame,text="15", variable=p['attr'], value="15", command = ShowAtr, indicatoron=0); At15.grid(row=1, column=2)
At20=Radiobutton(E2d_frame,text="20", variable=p['attr'], value="20", command = ShowAtr, indicatoron=0); At20.grid(row=1, column=3)
At25=Radiobutton(E2d_frame,text="25", variable=p['attr'], value="25", command = ShowAtr, indicatoron=0); At25.grid(row=1, column=4)
At30=Radiobutton(E2d_frame,text="30", variable=p['attr'], value="30", command = ShowAtr, indicatoron=0); At30.grid(row=1, column=5)
At35=Radiobutton(E2d_frame,text="35", variable=p['attr'], value="35", command = ShowAtr, indicatoron=0); At35.grid(row=1, column=6)
At40=Radiobutton(E2d_frame,text="40", variable=p['attr'], value="40", command = ShowAtr, indicatoron=0); At40.grid(row=1, column=7)
At45=Radiobutton(E2d_frame,text="45", variable=p['attr'], value="45", command = ShowAtr, indicatoron=0); At45.grid(row=1, column=8)
At50=Radiobutton(E2d_frame,text="50", variable=p['attr'], value="50", command = ShowAtr, indicatoron=0); At50.grid(row=1, column=9)
refreshLists()

''' RIGHT 3 Graphical Type'''
Radiobutton(E3_frame, text="Patterns",    variable=p['PatTen'], value="0", command = Show, indicatoron=0).grid(column=0, row=0, pady=5, padx=5)
Radiobutton(E3_frame, text="Clustering",  variable=p['PatTen'], value="1", command = Show, indicatoron=0).grid(column=0, row=1, pady=5, padx=5)
Radiobutton(E3_frame, text="Occ. 0D",     variable=p['PatTen'], value="2", command = Show, indicatoron=0).grid(column=0, row=2, pady=5, padx=5)
Radiobutton(E3_frame, text="Occ. 1D P",   variable=p['PatTen'], value="3", command = Show, indicatoron=0).grid(column=0, row=3, pady=5, padx=5)
Radiobutton(E3_frame, text="Occ. 1D G",   variable=p['PatTen'], value="4", command = Show, indicatoron=0).grid(column=0, row=4, pady=5, padx=5)
Radiobutton(E3_frame, text="Occ. 2D",     variable=p['PatTen'], value="5", command = Show, indicatoron=0).grid(column=0, row=5, pady=5, padx=5)
Radiobutton(E3_frame, text="Clus Quant",  variable=p['PatTen'], value="12",command = Show, indicatoron=0).grid(column=0, row=6, pady=5, padx=5)
Radiobutton(E3_frame, text="Correl.",     variable=p['PatTen'], value="14",command = Show, indicatoron=0).grid(column=0, row=7, pady=5, padx=5)
Radiobutton(E3_frame, text="Quantity",    variable=p['PatTen'], value="6", command = Show, indicatoron=0).grid(column=1, row=0, pady=5, padx=5)
Radiobutton(E3_frame, text="Densities",   variable=p['PatTen'], value="7", command = Show, indicatoron=0).grid(column=1, row=1, pady=5, padx=5)
Radiobutton(E3_frame, text="Dens Hist",   variable=p['PatTen'], value="8", command = Show, indicatoron=0).grid(column=1, row=2, pady=5, padx=5)
Radiobutton(E3_frame, text="All Den Hist",variable=p['PatTen'], value="11",command = Show, indicatoron=0).grid(column=1, row=3, pady=5, padx=5)
Radiobutton(E3_frame, text="Histog??",    variable=p['PatTen'], value="9", command = Show, indicatoron=0).grid(column=1, row=4, pady=5, padx=5)
Radiobutton(E3_frame, text="All Occ Hist",variable=p['PatTen'], value="10",command = Show, indicatoron=0).grid(column=1, row=5, pady=5, padx=5)
Radiobutton(E3_frame, text="avg GOF",     variable=p['PatTen'], value="13",command = Show, indicatoron=0).grid(column=1, row=6, pady=5, padx=5)

''' RIGHT 4 Displaying'''
textentry(E4aL_frame, p['max'], 'max', 7)
textentry(E4aL_frame, p['seuil'], 'seuil', 7)
textentry(E4aR_frame, p['hide'], 'hide', 7)
textentry(E4aR_frame, p['Ncont'], 'nCont', 7)
textentry(E4aR_frame, p['cclbl'], 'cclbl', 7)
listc= Listbox(E4b_frame, width=20, height=7, exportselection=0); listc.pack(side='top')
for nW in ['Blues','BuGn','BuPu','GnBu','Greens','Greys','Oranges','OrRd','PuBu','PuBuGn','PuRd','Purples','RdPu','Reds','YlGn','YlGnBu','YlOrBr','YlOrRd','',\
           'gist_earth','gist_gray','gist_heat','gist_ncar','gist_rainbow','gist_stern','gist_yarg','',\
           'autumn','bone','cool','copper','flag','gray','hot','hsv','jet','pink','prism','spectral','spring','summer','winter','',\
           'BrBG','PiYG','PRGn','PuOr','RdBu','RdGy','RdYlBu','RdYlGn','Spectral']:
    listc.insert(END, '%s' %nW)
listc.bind('<Double-1>', listcDef)
Checkbutton(E4c_frame, text="axis", variable=p['axis']).grid(column=0, row=0)
Checkbutton(E4c_frame, text="switch", variable=p['switch']).grid(column=1, row=0)
Checkbutton(E4c_frame, text="c_rev", variable=p['c_rev']).grid(column=2, row=0)
Checkbutton(E4c_frame, text="POGtoPG", variable=p['POG_PG']).grid(column=0, row=1)
Checkbutton(E4c_frame, text="norm", variable=p['norm']).grid(column=1, row=1)

Button(E_frame, text='Exit', command=Quit).pack(side='bottom')
root.mainloop()
root.destroy()
