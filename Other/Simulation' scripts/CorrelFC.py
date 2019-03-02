#!/usr/bin/env python
#-*- coding:Utf-8 -*-
''' Affiche les patterns, leur densité, leur répartition ou leur tendance en fonction des jeux de paramètres
    @author: mathieu.golos@gmail.com
'''

'#################################################  LIBRARIES  ##################################################'
import os
import copy
from pylab import figure, array, zeros, arange, get_cmap, corrcoef, where, average, size, linspace, pi, cm, log
from pylab import log, imread, subplot, rand, unravel_index, ones, sqrt, norm, newaxis, isfinite, diag, log10
from pylab import isreal, concatenate, colorbar, subplots
from numpy import load as np_load, save as np_save
from pickle import load as pk_load, dump as pk_dump
from matplotlib import colors, rcParams, gridspec
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.ticker import FormatStrFormatter
from scipy.stats import ss
from Tkinter import *
from time import time
from Tools.ext import data2array, Pdir
from Tools.matrices import triSupPearsonCor
rcParams.update({'font.size': 9})


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

def listdDef(e):
    p['dir_'].set(p['dir_ref'].get() + listd.get(listd.curselection()))
    refreshLists()
    Show()

def findFolderParameters(adress):

    l0 = os.listdir(adress)
    for i in range(len(l0))[::-1]:
        if not '.npy' in l0[i]:
            l0.remove(l0[i])

    lp1,lp2 = [],[]
    for i in range(len(l0)):
        splited = l0[i].replace('.npy','').rsplit('_')
        if splited[3] not in lp1:
            lp1.append( splited[3] )
        if splited[5] not in lp2:
            lp2.append( splited[5] )
    lp1 = array(lp1, dtype=float);  lp1.sort()
    lp2 = array(lp2, dtype=float);  lp2.sort()
    return lp1,lp2

def repertoire(CF=None, N='66'):
    dir_ = p['dir_'].get()

    if CF == 'CO':
        dir_ += '/CO'
    if CF == 'FC':
        dir_ += '/FC_%s' %N
    return dir_

def nameBoolTest(word, *args):
    mot = ''
    for i in range(len(args)):
        mot += args[i] * (args[i] in word)
    return mot

def paramNames():
    param= {}
    param['P']      = r'$P$'
    param['G']      = r'$G$'
    param['Psi']    = r'$\Psi$'
    param['Omega']  = r'$\Omega$'
    param['tauT']   = r'$\tau_\theta$'
    param['sigma2x']= r'$\sigma^2_x$'
    param['sigmax'] = r'$\sigma^2_x$'
    param['sigma2T']= r'$\sigma^2_\theta$'
    param['sigmaT'] = r'$\sigma^2_\theta$'
    param['sigmaA'] = r'$\sigma^2_A$'

    for test in os.listdir( repertoire(CF='FC') ):
        if '.npy' in test:
            p['p1str'].set(test.replace('.npy','').rsplit('_')[2])
            p['p2str'].set(test.replace('.npy','').rsplit('_')[4])
            break
    p['xlab'].set(param[p['p1str'].get()])
    p['ylab'].set(param[p['p2str'].get()])

def getOtherParam():
    param           = {}
    param['P']      = ["'-P','",        None, r'$P$']
    param['Psi']    = ["'-Psi','",      None, r'$\Psi$']
    param['G']      = ["'-G','",        None, r'$G$']
    param['Omega']  = ["'-Omega','",    None, r'$\Omega$']
    param['tauT']   = ["'-tauT','",     None, r'$\tau_\theta$']
    param['sigma2x']= ["'-sigma2x','",  None, r'$\sigma^2_x$']
    param['sigmax'] = ["'-sigmax','",   None, r'$\sigma^2_x$']
    param['sigma2T']= ["'-sigma2T','",  None, r'$\sigma^2_\theta$']
    param['sigmaT'] = ["'-sigmaT','",   None, r'$\sigma^2_\theta$']
    param['Cnnctm'] = ["'-Cnnctm','",   None]

    del param[p['p1str'].get()]
    del param[p['p2str'].get()]

    try:
        ofile = open(p['dir_'].get() + '/A_Tc.py', 'r')
        tfile = ofile.readlines()
        ofile.close()

        for line in tfile:
            for par in param.keys():
                if param[par][0] in line:
                    if par == 'Cnnctm':
                        param[par][1] = line.rsplit('/')[-1].rsplit("'")[0]
                        print 'Connectome:', param[par][1]
                    else:
                        param[par][1] = float(line.replace("'","").rsplit(",")[1])
    except:
        pass

    for par in param.keys():
        if param[par][1] == None:
            del param[par]
        #else:
            #param[par] = param[par][1]
    return param

def param2Darray(lp1,lp2):
    dir_    = repertoire(CF='FC',N=N)
    p2Dlist = []
    p1str   = p['p1str'].get()
    p2str   = p['p2str'].get()
    for i in xrange(len(lp1)):
        p2Dlist.append([])
        for j in xrange(len(lp2)):
            p2Dlist[i].append([])
            if len(str(lp2[1]-lp2[0]))-2 > 2:
                p2Dlist[i][j].append( dir_ + '/FC_A998_%s_%.2f_%s_%.3f' %(p1str,lp1[i],p2str,lp2[j]) )
                p2Dlist[i][j].append( dir_ +  '/FC_A66_%s_%.2f_%s_%.3f' %(p1str,lp1[i],p2str,lp2[j]) )
                p2Dlist[i][j].append( dir_ + '/FC_B998_%s_%.2f_%s_%.3f' %(p1str,lp1[i],p2str,lp2[j]) )
                p2Dlist[i][j].append( dir_ +  '/FC_B66_%s_%.2f_%s_%.3f' %(p1str,lp1[i],p2str,lp2[j]) )
            else:
                p2Dlist[i][j].append( dir_ + '/FC_A998_%s_%.2f_%s_%.2f' %(p1str,lp1[i],p2str,lp2[j]) )
                p2Dlist[i][j].append( dir_ +  '/FC_A66_%s_%.2f_%s_%.2f' %(p1str,lp1[i],p2str,lp2[j]) )
                p2Dlist[i][j].append( dir_ + '/FC_B998_%s_%.2f_%s_%.2f' %(p1str,lp1[i],p2str,lp2[j]) )
                p2Dlist[i][j].append( dir_ +  '/FC_B66_%s_%.2f_%s_%.2f' %(p1str,lp1[i],p2str,lp2[j]) )
    return p2Dlist

def Abreviations(abrev):
    abre_ind = ['r'+R for R in abrev]
    for i in xrange(33):
        abre_ind.append('l'+ abrev[::-1][i])
    return abre_ind

def mkdir(adress):
    ''' Crée un dossier s'il n'existe pas
    '''
    if not os.path.exists(adress):
        os.mkdir(adress)

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

def parcellation(vec, o, color, grid, maxs=None):
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

def clickParcellation(event, mat, o, color, grid, first=False, **kwargs):
    mapCanvas(o)
    if first:
        parcellation(mat, o, color, grid, **kwargs)
    else:
        print 'indice %i' %(xdata + 0.5)
        parcellation(mat[int(xdata + 0.5)], o, color, grid, **kwargs)
    p['C'][o].draw()

def permitSwitches(mat, color, o1, o2, **kwargs):
    p['C'][o1].get_tk_widget().bind("<Button-2>", \
                lambda event: clickParcellation(event, mat, o2, color, **kwargs))
    p['C'][o2].get_tk_widget().bind("<Button-2>", lambda event: Show())

def getStatusXY(event):
    global xdata, ydata
    xdata = event.xdata
    ydata = event.ydata

def mapCanvas(o):
    if not p['C'][o].get_tk_widget().winfo_ismapped():
        for i in p['C'].keys():
            p['C'][i].get_tk_widget().pack_forget()
            p['T'][i].pack_forget()
        p['C'][o].get_tk_widget().pack()
        p['T'][o].pack()

def fastCovariance(TC):
    TCm = (TC - TC.mean(axis=0)).T
    N   = TC.shape[1]
    cov = zeros((N, N))
    for i in xrange(N):
        cov[i,i:]  = TCm[i:].dot(TCm[i].T)
        cov[i:,i]  = cov[i,i:]
        #print i,'/',N
    return cov, TCm

def fastPearsonCorrelation(TC):
    N        = TC.shape[1]
    corr,TCm = fastCovariance(TC)
    TCss     = sqrt(ss(TCm, axis=1))
    for i in xrange(N):
        corr[i,i:] /= TCss[i:] * TCss[i]
        corr[i:,i]  = corr[i,i:]
    return where(isfinite(corr), corr, 0.)


'#############################################  MAIN FUNCTIONS  ##############################################'
def refreshDir():
    global listd
    listd.destroy()
    listd= Listbox(E1a_frame, width=20, exportselection=0);   listd.pack(side='left', pady=5)
    listoslistdir = os.listdir(p['dir_ref'].get())
    listoslistdir.sort()
    for nW in listoslistdir:
        if os.path.isdir(p['dir_ref'].get() + nW):
            listd.insert(END, nW)
    listd.bind('<Double-1>', listdDef)

def refreshLists():
    global lista, listb, lp1, lp2
    lp1,lp2 = findFolderParameters(p['dir_'].get() + '/FC_998')

    lista.destroy()
    listb.destroy()
    lista= Listbox(E1b_frame, width=10, exportselection=0);   lista.pack(side='left', pady=5)
    listb= Listbox(E1b_frame, width=10, exportselection=0);   listb.pack(side='left', pady=5)
    for nW in lp1:
        lista.insert(END, '%.3f' %nW)
    for nW in lp2:
        listb.insert(END, '%.3f' %nW)
    lista.bind('<Double-1>', listaDef)
    listb.bind('<Double-1>', listbDef)

    p['p1'].set(lp1[0])
    p['p2'].set(lp2[0])

    paramNames()
    p2Dlist = param2Darray(lp1,lp2)

def matricesCorrelation(FC1, FC2, posit, avg=True):
    ''' Return the correlation between the two matrices, defined by the mean across
    correlations for each columns. If posit=True, the correlation is done for positive matrices.
    '''
    from pylab import isfinite
    if FC1.shape != FC2.shape:
        print 'Functional Connectomes shape are not the same'
        return 0
    l = FC1.shape[0]
    corr = zeros(l)
    if posit:
        for c in range(l):
            corr[c] = corrcoef(where(FC1[c]<0,0,FC1[c]), where(FC2[c]<0,0,FC2[c]))[0,1]
    else:
        for c in range(l):
            corr[c] = corrcoef(FC1[c], FC2[c])[0,1]
    corr = where(isfinite(corr), corr, 0.)

    if avg: return corr.mean()
    else:   return corr

def compareMatrices(fig, lMat,lTitl, aspect=None, switch=0, color='jet'):
    nb = len(lMat)
    nc = int(ceil(sqrt(nb)))
    nl = int(floor(nb / float(nc)))
    fig.clf()
    for i in xrange(nb):
        if switch: p['A'][0] = fig.add_subplot(nc,nl,i+1)
        else:      p['A'][0] = fig.add_subplot(nl,nc,i+1)
        if aspect:
            p['A'][0].imshow(lMat[i], interpolation = 'nearest', aspect='auto', cmap=color)
        else:
            p['A'][0].imshow(lMat[i], interpolation = 'nearest', aspect = 1, cmap=color)
        p['A'][0].title(lTitl[i])
        p['A'][0].colorbar()
    fig.draw()

def TwoTriangleFunCon(FC1, FC2, posit=False, diago=False, uniHemi=False):
    ''' Show two functional connectivity matrices using only one square matrix
    FC1 is upper right and FC2 bottom left
    '''
    if FC1.shape != FC2.shape:
        print 'Functional Connectomes shape are not the same'
        return 0
    FCC = zeros(FC1.shape)
    if p['norm'].get():
        no1 = norm(FC1)
        no2 = norm(FC2)
    else:
        no1 = 1.
        no2 = 1.

    if uniHemi:
        for i in range(FCC.shape[0]):
            FCC[i,i:]   = FC1[i,i:] / no1
            FCC[i+1:,i] = FC2[i+1:,i] / no2
    else:
        for i in range(FCC.shape[0]):
            FCC[i,i:]   = FC1[i,i:] / no1
            FCC[i+1:,i] = FC2[i+1:,i] / no2
        if not diago:
            #FCC -= diag(diag(FCC))
            FCC[diag(FCC.shape) == 1] = FCC.min()

    if posit: return where(FCC < 0.0, 0., FCC)
    else:     return FCC

def loadEFC():
    global p
    dir_ = Pdir('Main/Connectomes/Hagmann/')
    for i in range(6):
        p['E%i_998'%i] = data2array(dir_+ 'FC_D_998_%s.npy'%i)
        p['E%i_66' %i] = data2array(dir_+ 'FC_D_66_%s.npy' %i)
        p['SC%i_998'%i] = data2array(dir_+ 'SC_D_998_%s.npy'%i)
        p['SC%i_66' %i] = data2array(dir_+ 'SC_D_66_%s.npy' %i)

def getSlices():
    if p['LR'].get() == 0:
        sl1 = slice(0,998), slice(0,998), slice(0,998), slice(0,998)
        sl2 = slice(0,66), slice(0,66), slice(0,66), slice(0,66)
    elif p['LR'].get() == 1:
        sl1 = slice(0,499), slice(0,499), slice(997,498,-1), slice(997,498,-1)
        sl2 = slice(0,33), slice(0,33), slice(65,32,-1), slice(65,32,-1)
    elif p['LR'].get() == 2:
        sl1 = slice(499,998), slice(0,499), slice(-500,-999,-1), slice(997,498,-1)
        sl2 = slice(33,66), slice(0,33), slice(-34,-67,-1), slice(65,32,-1)
    return sl1, sl2

def getFCCorrelation():
    simul = [   p['dir_'].get() + '/FC_66/FC_A66_', \
                p['dir_'].get() + '/FC_998/FC_A998_', \
                p['dir_'].get() + '/FC_66/FC_B66_', \
                p['dir_'].get() + '/FC_998/FC_B998_']
    dirCorr = p['dir_'].get() + '/CO'
    heminame = '' + (p['LR'].get() != 0) * '_uni' \
                  + (p['LR'].get() == 2) * '_ant' \
                  + (p['SC'].get() == 1) * '_SC' \
                  + (p['log'].get() == 1) * '_log'
    subChaine = '%i_%s_%s%s' %(p['patient'].get(), p['p1str'].get(), p['p2str'].get(), heminame)
    seuil     = p['seuil'].get()

    mkdir(dirCorr)
    CO_A66, CO_A998, CO_B66, CO_B998 = None, None, None, None

    if p['AB'].get() != 1:
        try:
            if p['force'].get(): goToExcept
            CO_B66  = data2array(dirCorr + '/CO_B66_%s.npy'  %subChaine)
            CO_B998 = data2array(dirCorr + '/CO_B998_%s.npy' %subChaine)
            print dirCorr + '/CO_B998_%s.npy' %subChaine
        except:
            shapeTS = len(lp1), len(lp2)
            CO_B66  = zeros(shapeTS)
            CO_B998 = zeros(shapeTS)
            sl1,sl2 = getSlices()
            if p['SC'].get():
                FC_E998 = array(meanHemisph(p['SC%i_998'%p['patient'].get()], sl1))
                FC_E66  = array(meanHemisph(p['SC%i_66' %p['patient'].get()], sl2))
            else:
                FC_E998 = array(meanHemisph(p['E%i_998'%p['patient'].get()], sl1))
                FC_E66  = array(meanHemisph(p['E%i_66' %p['patient'].get()], sl2))

            FC_E998 /= norm(FC_E998)
            FC_E66  /= norm(FC_E66)

            if p['log'].get() == 1:
                FC_E998 = log10(FC_E998)
                FC_E66  = log10(FC_E66)
                FC_E998[isfinite(FC_E998) == False] = seuil
                FC_E66 [isfinite(FC_E66 ) == False] = seuil
                FC_E998[FC_E998 < seuil] = seuil
                FC_E66 [FC_E66  < seuil] = seuil

            for p1 in range(len(lp1)):
                for p2 in range(len(lp2)):
                    if p['p1str'].get() == 'sigmaA':
                        name    = '%s_%.3f_%s_%.2f.npy' %(p['p1str'].get(), lp1[p1], p['p2str'].get(), lp2[p2])
                    else:
                        name    = '%s_%.2f_%s_%.2f.npy' %(p['p1str'].get(), lp1[p1], p['p2str'].get(), lp2[p2])
                    try:
                        if p['forceFC'].get():
                            adress = simul[2].replace('FC_B','BO_').replace('FC','BO') + name
                            B66    = data2array(adress)[10:]
                            FC_B66 = meanHemisph(fastPearsonCorrelation(B66), sl2)
                            np_save(adress, B66)
                            del B66
                        else:
                            FC_B66 = meanHemisph(data2array(simul[2] + name), sl2)

                        FC_B66 /= norm(FC_B66)
                        if p['log'].get() == 1:
                            FC_B66 = log10(FC_B66)
                            FC_B66[isfinite(FC_B66) == False] = seuil
                            FC_B66[FC_B66 < seuil] = seuil
                        if p['posit'].get():
                            CO_B66 [p1,p2] = triSupPearsonCor(where(FC_B66<0,0,FC_B66), where(FC_E66<0,0,FC_E66))
                        else:
                            CO_B66 [p1,p2] = triSupPearsonCor(FC_B66 ,FC_E66)
                            
                    except: print 'no file'

                    try:
                        if p['forceFC'].get():
                            adress  = simul[3].replace('FC_B','BO_').replace('FC','BO') + name
                            B998    = data2array(adress)[10:]
                            FC_B998 = meanHemisph(fastPearsonCorrelation(B998), sl1)
                            np_save(adress, B998)
                            del B998
                        else:
                            FC_B998 = meanHemisph(data2array(simul[3] + name), sl1)

                        FC_B998 /= norm(FC_B998)
                        if p['log'].get() == 1:
                            FC_B998 = log10(FC_B998)
                            FC_B998[isfinite(FC_B998) == False] = seuil
                            FC_B998[FC_B998 < seuil] = seuil
                        if p['posit'].get():
                            CO_B998 [p1,p2] = triSupPearsonCor(where(FC_B998<0,0,FC_B998), where(FC_E998<0,0,FC_E998))
                        else:
                            CO_B998 [p1,p2] = triSupPearsonCor(FC_B998 ,FC_E998)
                    except: print 'no file'
                    print '%i/%i ; %i/%i' %(p1+1, len(lp1), p2+1, len(lp2))
            print dirCorr + '/CO_B998_%s.npy' %subChaine
            np_save(dirCorr + '/CO_B66_%s.npy'  %subChaine, CO_B66)
            np_save(dirCorr + '/CO_B998_%s.npy' %subChaine, CO_B998)
            del FC_B998, FC_B66, FC_E998, FC_E66

    if p['AB'].get() != 2:
        try:
            if p['force'].get(): goToExcept
            CO_A66  = data2array(dirCorr + '/CO_A66_%s.npy'  %subChaine)
            CO_A998 = data2array(dirCorr + '/CO_A998_%s.npy' %subChaine)
        except:
            shapeTS = len(lp1), len(lp2)
            CO_A66  = zeros(shapeTS)
            CO_A998 = zeros(shapeTS)
            sl1,sl2 = getSlices()
            if p['SC'].get():
                FC_E998 = array(meanHemisph(p['SC%i_998'%p['patient'].get()], sl1))
                FC_E66  = array(meanHemisph(p['SC%i_66' %p['patient'].get()], sl2))
            else:
                FC_E998 = array(meanHemisph(p['E%i_998'%p['patient'].get()], sl1))
                FC_E66  = array(meanHemisph(p['E%i_66' %p['patient'].get()], sl2))

            FC_E998 /= norm(FC_E998)
            FC_E66  /= norm(FC_E66)

            if p['log'].get() == 1:
                FC_E998 = log10(FC_E998)
                FC_E66  = log10(FC_E66)
                FC_E998[isfinite(FC_E998) == False] = seuil
                FC_E66 [isfinite(FC_E66 ) == False] = seuil
                FC_E998[FC_E998 < seuil] = seuil
                FC_E66 [FC_E66  < seuil] = seuil

            for p1 in range(len(lp1)):
                for p2 in range(len(lp2)):
                    if p['p1str'].get() == 'sigmaA':
                        name    = '%s_%.3f_%s_%.2f.npy' %(p['p1str'].get(), lp1[p1], p['p2str'].get(), lp2[p2])
                    else:
                        name    = '%s_%.2f_%s_%.2f.npy' %(p['p1str'].get(), lp1[p1], p['p2str'].get(), lp2[p2])
                    try:
                        FC_A66  = meanHemisph(data2array(simul[0] + name), sl2)
                        FC_A66 /= norm(FC_A66)
                        if p['log'].get() == 1:
                            FC_A66 = log10(FC_A66)
                            FC_A66[isfinite(FC_A66) == False] = seuil
                            FC_A66[FC_A66 < seuil] = seuil
                        if p['posit'].get():
                            CO_A66 [p1,p2] = triSupPearsonCor(where(FC_A66<0,0,FC_A66), where(FC_E66<0,0,FC_E66))
                        else:
                            CO_A66 [p1,p2] = triSupPearsonCor(FC_A66 ,FC_E66)
                    except: print 'no file'
                    try:
                        FC_A998 = meanHemisph(data2array(simul[1] + name), sl1)
                        FC_A998 /= norm(FC_A998)
                        if p['log'].get() == 1:
                            FC_A998 = log10(FC_A998)
                            FC_A998[isfinite(FC_A998) == False] = seuil
                            FC_A998[FC_A998 < seuil] = seuil
                        if p['posit'].get():
                            CO_A998 [p1,p2] = triSupPearsonCor(where(FC_A998<0,0,FC_A998), where(FC_E998<0,0,FC_E998))
                        else:
                            CO_A998 [p1,p2] = triSupPearsonCor(FC_A998 ,FC_E998)
                    except: print 'no file'
                    print '%i/%i ; %i/%i' %(p1+1, len(lp1), p2+1, len(lp2))
            CO_A998[0,0] = (CO_A998[0,1]+CO_A998[1,0])/2.
            CO_A998[4,4] = (CO_A998[4,2]+CO_A998[4,3])/2.

            np_save(dirCorr + '/CO_A66_%s.npy'  %subChaine, CO_A66)
            np_save(dirCorr + '/CO_A998_%s.npy' %subChaine, CO_A998)
            del FC_A998, FC_A66, FC_E998, FC_E66

    return CO_A66, CO_A998, CO_B66, CO_B998

def meanHemisph(FC,sli):
    if sli[0:2] == sli[2:4]:
        return FC
    else:
        return (FC[sli[0:2]] + FC[sli[2:4]]) *.5

def showCorr():
    CO_A66, CO_A998, CO_B66, CO_B998 = getFCCorrelation()

    minarg = []
    maxarg = []
    if p['AB'].get() != 2 and p['All99866'].get() != 2: # A998
        minarg.append(CO_A998.min())
        maxarg.append(CO_A998.max())
    if p['AB'].get() != 2 and p['All99866'].get() != 1: # A66
        minarg.append(CO_A66.min())
        maxarg.append(CO_A66.max())
    if p['AB'].get() != 1 and p['All99866'].get() != 2: # B998
        minarg.append(CO_B998.min())
        maxarg.append(CO_B998.max())
    if p['AB'].get() != 1 and p['All99866'].get() != 1: # B66
        minarg.append(CO_B66.min())
        maxarg.append(CO_B66.max())
    minarg.append(p['min'].get())
    maxarg.append(p['max'].get())
    vmin = min(minarg)
    vmax = max(maxarg)

    if p['AB'].get() != 1:
        argm_CO_B998 = unravel_index(CO_B998.argmax(), CO_B998.shape)
        argm_CO_B66  = unravel_index(CO_B66.argmax(),  CO_B66.shape)
        print lp1[argm_CO_B998[0]], lp2[argm_CO_B998[1]], 'CO_B998: %.3f' %CO_B998.max()
        print lp1[argm_CO_B66[0]],  lp2[argm_CO_B66[1]],  'CO_B66: %.3f' %CO_B66.max()

    p['F'][o].clear()
    nbl = 1 + (p['AB'].get() == 0)
    nbc = 1 + (p['All99866'].get() == 0)
    grid = AxesGrid(p['F'][o], 111,
                nrows_ncols = (nbl,nbc),
                axes_pad = 0.3,
                aspect = True,
                share_all = True,
                label_mode = 'L',
                cbar_location = 'right',
                cbar_mode = 'single')
    cmap = get_cmap(p['cmap'].get()+p['c_rev'].get()*'_r')
    kwargs = {'vmin': vmin, \
              'vmax': vmax, \
              'cmap': cmap, \
              #'aspect': '1', \
              'origin': 'lower', \
              'interpolation': p['interp'].get()}
    i = 0
    while i < nbl * nbc:
        if p['AB'].get() != 2 and p['All99866'].get() != 2:
            im = grid[i].imshow(CO_A998.T, **kwargs)
            grid[i].set_title(r'$Activity\ 998\ Nodes$')
            i += 1
        if p['AB'].get() != 2 and p['All99866'].get() != 1:
            im = grid[i].imshow(CO_A66.T,  **kwargs)
            grid[i].set_title(r'$Activity\ 66\ Nodes$')
            i += 1
        if p['AB'].get() != 1 and p['All99866'].get() != 2:
            im = grid[i].imshow(CO_B998.T, **kwargs)
            grid[i].set_title(r'$BOLD\ 998\ Nodes$')
            i += 1
        if p['AB'].get() != 1 and p['All99866'].get() != 1:
            im = grid[i].imshow(CO_B66.T,  **kwargs)
            grid[i].set_title(r'$BOLD\ 66\ Nodes$')
            i += 1
    grid.axes_llc.set_xlabel(p['xlab'].get(), fontsize=20)
    grid.axes_llc.set_ylabel(p['ylab'].get(), fontsize=20)
    grid.axes_llc.set_xticks(range(len(lp1))[::2])
    grid.axes_llc.set_yticks(range(len(lp2))[::2])
    grid.axes_llc.set_xticklabels(lp1[::2])
    grid.axes_llc.set_yticklabels(lp2[::2])
    grid.cbar_axes[0].colorbar(im)
    del CO_A66, CO_A998, CO_B66, CO_B998

def getFC():
    sl1,sl2 = getSlices()
    name    = '%s_%.2f_%s_%.2f.npy' %(p['p1str'].get(), p['p1'].get(), p['p2str'].get(), p['p2'].get())
    seuil   = p['seuil'].get()

    # Simulated
    if p['forceFC'].get() == 1:
        B9   = data2array(p['dir_'].get() + '/BO_998/BO_998_' + name)[10:]
        B6   = data2array(p['dir_'].get() + '/BO_66/BO_66_' + name)[10:]
        print B9.shape, B6.shape
        FCS9 = fastPearsonCorrelation(B9)
        FCS6 = fastPearsonCorrelation(B6)
        del B9, B6
    else:
        FCS9 = data2array(p['dir_'].get() + '/FC_998/FC_B998_' + name)
        FCS6 = data2array(p['dir_'].get() + '/FC_66/FC_B66_' + name)
    FCS9 = meanHemisph(FCS9, sl1)
    FCS6 = meanHemisph(FCS6, sl2)

    # Empirical
    if p['SC'].get():
        FCE9 = array(meanHemisph(p['SC%i_998'%p['patient'].get()], sl1))
        FCE6 = array(meanHemisph(p['SC%i_66' %p['patient'].get()], sl2))
    else:
        FCE9 = array(meanHemisph(p['E%i_998'%p['patient'].get()], sl1))
        FCE6 = array(meanHemisph(p['E%i_66' %p['patient'].get()], sl2))


    FCS9 /= norm(FCS9)
    FCS6 /= norm(FCS6)
    FCE9 /= norm(FCE9)
    FCE6 /= norm(FCE6)

    # Other
    if p['posit'].get():
        FCS9[FCS9 < 0] = 0
        FCS6[FCS6 < 0] = 0
        FCE9[FCE9 < 0] = 0
        FCE6[FCE6 < 0] = 0

    if p['log'].get() == 1:
        FCS9 = log10(FCS9)
        FCS6 = log10(FCS6)
        FCE9 = log10(FCE9)
        FCE6 = log10(FCE6)
        FCS9[isfinite(FCS9) == False] = seuil
        FCS6[isfinite(FCS6) == False] = seuil
        FCE9[isfinite(FCE9) == False] = seuil
        FCE6[isfinite(FCE6) == False] = seuil
        FCS9[FCS9 < seuil] = seuil
        FCS6[FCS6 < seuil] = seuil
        FCE9[FCE9 < seuil] = seuil
        FCE6[FCE6 < seuil]

    return FCS9, FCS6, FCE9, FCE6

def showFC():
    FCS9, FCS6, FCE9, FCE6 = getFC()
    if p['posit'].get():
        corr1 = triSupPearsonCor(where(FCS9<0,0,FCS9), where(FCE9<0,0,FCE9))
        corr2 = triSupPearsonCor(where(FCS6<0,0,FCS6), where(FCE6<0,0,FCE6))
    else:
        corr1 = triSupPearsonCor(FCS9 ,FCE9)
        corr2 = triSupPearsonCor(FCS6 ,FCE6)

    if p['diag'].get() == 0:
        FCS9[diag(ones(FCS9.shape[0])) == 1] = FCS9.min()
        FCS6[diag(ones(FCS6.shape[0])) == 1] = FCS6.min()
        FCE9[diag(ones(FCE9.shape[0])) == 1] = FCE9.min()
        FCE6[diag(ones(FCE6.shape[0])) == 1] = FCE6.min()

    if p['SqTri'].get() == 1:
        FCS9 = TwoTriangleFunCon(FCS9, FCE9, posit=False, diago=p['diag'].get(), uniHemi=False)
        FCS6 = TwoTriangleFunCon(FCS6, FCE6, posit=False, diago=p['diag'].get(), uniHemi=False)

    p['F'][o].clear()
    AX  = [[],[],[],[]]
    nbl = 2 - p['SqTri'].get()
    nbc = (p['All99866'].get() == 0) * 1 + 1
    for i in range(nbl * nbc):
        AX[i] = p['F'][o].add_subplot(nbl,nbc,i+1)
    cmap = get_cmap(p['cmap'].get()+p['c_rev'].get()*'_r')
    cmap.set_bad('w',1.)

    if p['All99866'].get() == 0: # FCS998 + FCS66
        im = AX[0].imshow(FCS9, cmap=cmap, interpolation=p['interp'].get())
        im = AX[1].imshow(FCS6, cmap=cmap, interpolation=p['interp'].get())
        AX[0].set_title(r'$Simulated\ FC\ 998\ Nodes\ (corr:\ %.2f)$' %corr1)
        AX[1].set_title(r'$Simulated\ FC\ 66\ Nodes\ (corr:\ %.2f)$' %corr2)

        if p['SqTri'].get() == 0: # FCS998 + FCS66 + FCE998 + FCE66
            im = AX[2].imshow(FCE9, cmap=cmap, interpolation=p['interp'].get())
            im = AX[3].imshow(FCE6, cmap=cmap, interpolation=p['interp'].get())
            AX[2].set_title(r'$Empirical\ FC\ 998\ Nodes$')
            AX[3].set_title(r'$Empirical\ FC\ 66\ Nodes$')

    elif p['All99866'].get() == 1: # FCS998
        im = AX[0].imshow(FCS9, cmap=cmap, interpolation=p['interp'].get())
        AX[0].set_title(r'$Simulated\ FC\ 998\ Nodes\ (corr:\ %.2f)$' %corr1)

        if p['SqTri'].get() == 0: # FCS998 + FCE998
            im = AX[1].imshow(FCE9, cmap=cmap, interpolation=p['interp'].get())
            AX[1].set_title(r'$Empirical\ FC\ 998\ Nodes$')

    elif p['All99866'].get() == 2: # FCS66
        im = AX[0].imshow(FCS6, cmap=cmap, interpolation=p['interp'].get())
        AX[0].set_title(r'$Simulated\ FC\ 66\ Nodes\ (corr:\ %.2f)$' %corr2)

        if p['SqTri'].get() == 0: # FCS66 + FCE66
            im = AX[1].imshow(FCE6, cmap=cmap, interpolation=p['interp'].get())
            AX[1].set_title(r'$Empirical\ FC\ 66\ Nodes$')

    cbar_ax = p['F'][o].add_axes([0.85, 0.15, 0.05, 0.7])
    p['F'][o].subplots_adjust(right=0.8)
    p['F'][o].colorbar(im, cax=cbar_ax)

def showBOLDSIM():
    p['F'][o].clear()

    if p['All99866'].get() in [0,1]: N = 998
    else:                            N = 66

    name = '%s_%.2f_%s_%.2f.npy' %(p['p1str'].get(), p['p1'].get(), p['p2str'].get(), p['p2'].get())
    BO   = data2array(p['dir_'].get() + '/BO_%i/BO_%i_%s' %(N,N,name))
    if p['forceFC'].get():
        BO = BO[10:]
    BO  -= BO.mean(0)
    if p['AA'].get():
        BO /= abs(BO).max(0)

    ax = p['F'][o].add_subplot(111)
    ax.set_yticklabels(abre_ind)
    if p['All99866'].get() in [0,1]: ax.set_yticks(a998i)
    else:                            ax.set_yticks(xrange(66))

    im = ax.imshow(BO.T, cmap=get_cmap(p['cmap'].get()+p['c_rev'].get()*'_r'), \
                         aspect = 'auto', origin='lower', interpolation=p['interp'].get())
    ax.tick_params(labelright=True)
    ax.set_xlabel(r'$Time\ (s)$', fontsize = 12)
    permitSwitches(BO, cm_x, o, 'v', grid=p['A']['v'])
    del BO

def showBOLDEMP():
    p['F'][o].clear()

    if p['patient'].get():
        if p['All99866'].get() in [0,1]: N = 998
        else:                            N = 66

        adrs = Pdir('Main/TimeCourses')
        if p['AB'].get() == 0:
            BOa  = data2array(adrs + '/Hagmann/TC_%i_%ia.dat' %(N,p['patient'].get())).T
            BOb  = data2array(adrs + '/Hagmann/TC_%i_%ib.dat' %(N,p['patient'].get())).T
            BOa -= BOa.mean(0)
            BOb -= BOb.mean(0)
            BO   = concatenate((BOa, BOb), axis=0)
            del BOa, BOb
        elif p['AB'].get() == 1:
            BO  = data2array(adrs + '/Hagmann/TC_%i_%ia.dat' %(N,p['patient'].get())).T
            BO -= BO.mean(0)
        else:
            BO  = data2array(adrs + '/Hagmann/TC_%i_%ib.dat' %(N,p['patient'].get())).T
            BO -= BO.mean(0)

        if p['AA'].get():
            BO /= abs(BO).max(0)

        ax = p['F'][o].add_subplot(111)
        ax.set_yticklabels(abre_ind)
        if p['All99866'].get() in [0,1]: ax.set_yticks(a998i)
        else:                            ax.set_yticks(xrange(66))

        im = ax.imshow(BO.T, cmap=get_cmap(p['cmap'].get()+p['c_rev'].get()*'_r'), \
                            aspect = 'auto', origin='lower', interpolation=p['interp'].get())
        ax.tick_params(labelright=True)
        ax.set_xlabel(r'$Time\ (s)$', fontsize = 12)
        permitSwitches(BO, cm_x, o, 'v', grid=p['A']['v'])
        del BO
    else:
        print 'select a patient'

def showROICorr():
    p['F'][o].clear()
    ax   = p['F'][o].add_subplot(111)
    name = '%s_%.2f_%s_%.2f.npy' %(p['p1str'].get(), p['p1'].get(), p['p2str'].get(), p['p2'].get())
    if p['All99866'].get() <= 1: # FCS998
        if p['SC'].get():
            FCE = p['SC%i_998' %p['patient'].get()]
        else:
            FCE = p['E%i_998' %p['patient'].get()]
        FCS = data2array(p['dir_'].get() + '/FC_998/FC_B998_' + name)
        sli = getSlices()[0]
    else: # FCS66
        if p['SC'].get():
            FCE = p['SC%i_66' %p['patient'].get()]
        else:
            FCE = p['E%i_66' %p['patient'].get()]
        FCS = data2array(p['dir_'].get() + '/FC_66/FC_B66_' + name)
        sli = getSlices()[1]

    FCS  = meanHemisph(FCS, sli)
    FCS /= norm(FCS)
    FCE  = meanHemisph(FCE, sli)
    FCE /= norm(FCE)
    corr = matricesCorrelation(FCS,FCE,p['posit'].get(), avg=False)

    if p['All99866'].get() <= 1: # FCS998
        ax.set_xlim(1,998)
        ax.bar(range(1,999), corr)
    else:
        ax.set_xticks(range(66))
        ax.set_xticklabels(abre_ind, rotation=80)
        ax.bar(range(66), corr)
    ax.set_ylim(0,1)
    if p['SC'].get():
        ax.set_title(r'$Simulated\ /\ Empirical\ Functional\ Connectivity\ Correlation\ (mean:\ %.2f)$' %corr.mean(), fontsize=16)
    else:
        ax.set_title(r'$Simulated\ Functional\ /\ Structural\ Connectivity\ Correlation\ (mean:\ %.2f)$' %corr.mean(), fontsize=16)

    permitSwitches(corr, cm_x, o, 'v', grid=p['A']['v'], first=True, maxs=(0,1))

def Show():
    global o
    o = p['ShowTyp'].get()
    mapCanvas(o)
    othPar = getOtherParam()
    if   o == 0: showCorr()
    elif o == 1: showFC()
    elif o == 2: showBOLDSIM()
    elif o == 3: showROICorr()
    elif o == 4: showBOLDEMP()

    subtitle = ''
    for par in othPar.keys():
        if par != 'Cnnctm':
            subtitle += '%s (%s)  ' %(othPar[par][2], othPar[par][1])
    for i in range(len(p['F'])-1):
        p['F'][i].suptitle(subtitle, y=0.03, fontsize = 12)
    p['C'][o].show()


'#################################################  PARAMETERS  #################################################'
global abre_ind, a998i
root = Tk()
root.title("CorrelFC")
row_counter= 0
cm_x = cm.PRGn

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
abrev    = ['ENT','PARH','TP','FP','FUS','TT','LOCC','SP','IT','IP',\
         'SMAR','BTST','MT','ST','PSTC','PREC','CMF','POPE','PTRI','RMF',\
         'PORB','LOF','CAC','RAC','SF','MOF','LING','PCAL','CUN','PARC',\
         'ISTC','PCUN','PC']
abre_ind = Abreviations(abrev)
a998i    = [1, 7, 10, 12, 34, 37, 56, 83, 102, 130, 146,\
         153, 173, 201, 232, 268, 281, 291, 299, 321, 327, 346,\
         350, 354, 400, 412, 429, 439, 449, 461, 469, 492, 499,\
         506, 529, 537, 548, 556, 565, 581, 593, 643, 647, 651,\
         671, 677, 696, 703, 714, 727, 763, 793, 822, 841, 846,\
         865, 890, 907, 934, 956, 960, 982, 984, 988, 994, 997]
region= Pdir('Main/Regions/Images')
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
p['dir_ref'] = StringVar();         p['dir_ref'].set("/home/golos/Main/Simulations/TcHrfDowFc/")          # Dossier racine
p['dir_']    = StringVar();         p['dir_'].set("/home/golos/Main/Simulations/TcHrfDowFc/SC0_P_sx")     # Dossier de simulation
p['p1']    = DoubleVar();           p['p1'].set(1)                              # Valeur du parametre 1
p['p2']    = DoubleVar();           p['p2'].set(1)                              # Valeur du parametre 2
p['p1str'] = StringVar();           p['p1str'].set('')                          # x file label
p['p2str'] = StringVar();           p['p2str'].set('')                          # y file label
p['xlab'] = StringVar();            p['xlab'].set('')                           # x label
p['ylab'] = StringVar();            p['ylab'].set('')                           # y label

p['ShowTyp'] = IntVar();            p['ShowTyp'].set(0)                         # show type Corr/FC/BOLD
p['All99866']= IntVar();            p['All99866'].set(0)                        # 998 and/or 66
p['LR']      = IntVar();            p['LR'].set(0)                              # 2 hemispheres(0), LL/RR(1), LR/RL(2)
p['SqTri']   = IntVar();            p['SqTri'].set(1)                           # Total square matrix or mixed triangles matrices
p['AB']      = IntVar();            p['AB'].set(2)                              # from Activity / Bold signals
p['log']     = IntVar();            p['log'].set(0)                             # log10(FC)
p['patient'] = IntVar();            p['patient'].set(0)                         # Patient
p['SC']      = IntVar();            p['SC'].set(0)                              # Structural connectivity and not Functional

p['diag']    = IntVar();            p['diag'].set(0)                            # diagonal or not
p['posit']   = IntVar();            p['posit'].set(1)                           # positive values only for correlations
p['force']   = IntVar();            p['force'].set(0)                           # force to write on existing file
p['forceFC'] = IntVar();            p['forceFC'].set(0)                         # force calculation of FC without first 10s
p['c_rev']   = IntVar();            p['c_rev'].set(0)                           # Reverse colormap or not
p['norm']    = IntVar();            p['norm'].set(0)                            # Normalize the figure by the maximum on oen parameter
p['seuil']   = DoubleVar();         p['seuil'].set(-3)                          # log seuil for FC

p['min']     = DoubleVar();         p['min'].set(1)                             # Minimum for colorbar
p['max']     = DoubleVar();         p['max'].set(0)                             # Maximum for colorbar
p['cmap']    = StringVar();         p['cmap'].set("Spectral")                   # Couleur des graphiques
p['interp']  = StringVar();         p['interp'].set('nearest')                  # Interpolation for imshow 'nearest' or 'bicubic'
p['AA']      = IntVar();            p['AA'].set(1)                              # Average over time and divide by max over space

loadEFC()


'##################################################  DISPLAY  ###################################################'
W_frame   = LabelFrame(root, bd=0);                     W_frame.pack(side='left')

E_frame   = LabelFrame(root, bd=0);                     E_frame.pack(side='right', padx=5)
E1_frame  = LabelFrame(E_frame, text= "Directory");     E1_frame.pack(side='top', pady=4, fill=X)
E1a_frame = LabelFrame(E1_frame, bd=0);                 E1a_frame.pack(side='top')
E1b_frame = LabelFrame(E1_frame, bd=0);                 E1b_frame.pack(side='top')

E12_frame = LabelFrame(E_frame, relief=SUNKEN);         E12_frame.pack(side='top', pady=4)

E2_frame  = LabelFrame(E_frame, text= "Options");       E2_frame.pack(side='top', pady=4, fill=X)
E2b_frame = LabelFrame(E2_frame, bd=0);                 E2b_frame.pack(side='top')
E2c_frame = LabelFrame(E2_frame, bd=0);                 E2c_frame.pack(side='top')
E2d_frame = LabelFrame(E2_frame, bd=0);                 E2d_frame.pack(side='top')
E2e_frame = LabelFrame(E2_frame, bd=0);                 E2e_frame.pack(side='top')
E2f_frame = LabelFrame(E2_frame, bd=0);                 E2f_frame.pack(side='top')
E2g_frame = LabelFrame(E2_frame, bd=0);                 E2g_frame.pack(side='top')
E2h_frame = LabelFrame(E2_frame, bd=0);                 E2h_frame.pack(side='top')

E3_frame  = LabelFrame(E_frame);                        E3_frame.pack(side='top', pady=4, fill=X)
E3a_frame = LabelFrame(E3_frame, bd=0);                 E3a_frame.pack(side='top')
E3b_frame = LabelFrame(E3_frame, bd=0);                 E3b_frame.pack(side='top')
E3c_frame = LabelFrame(E3_frame, bd=0);                 E3c_frame.pack(side='top')
E3d_frame = LabelFrame(E3_frame, bd=0);                 E3d_frame.pack(side='top')


'##################################################  CANVAS  #################################################'
p['C'],p['F'],p['A'],p['T'] = {},{},{},{}
for i in [0,1,2,3,4,'v']:
    p['F'][i] = figure(figsize=(10.8, 9), dpi=100)
    p['F'][i].patch.set_facecolor('white')
    p['C'][i] = FigureCanvasTkAgg(p['F'][i], master=W_frame)
    p['C'][i].show()
    if i in [2,3]:
        p['C'][i].mpl_connect('motion_notify_event', getStatusXY)
    p['T'][i] = NavigationToolbar2TkAgg(p['C'][i], W_frame)
    p['T'][i].pack_forget()

p['C'][0].get_tk_widget().pack()
p['T'][0].pack()
p['F']['v'].clf()
p['A']['v'] = AxesGrid(p['F']['v'], 111,
                    nrows_ncols = (2,2),
                    axes_pad = 0.0,
                    share_all = True,
                    label_mode = 'L',
                    cbar_location = 'right',
                    cbar_mode = 'single')


'##################################################  WIDGETS  #################################################'
listd= Listbox(E1a_frame, width=10, exportselection=0);   listd.pack(side='left', pady=5)
refreshDir()
lista= Listbox(E1b_frame, width=10, exportselection=0);   lista.pack(side='left', pady=5)
listb= Listbox(E1b_frame, width=10, exportselection=0);   listb.pack(side='left', pady=5)
refreshLists()

Radiobutton(E12_frame, text="Corr.", variable=p['ShowTyp'], value=0, command = Show, indicatoron=0).grid(column=0, row=0)
Radiobutton(E12_frame, text="F.Co.", variable=p['ShowTyp'], value=1, command = Show, indicatoron=0).grid(column=1, row=0)
Radiobutton(E12_frame, text="ROIc",  variable=p['ShowTyp'], value=3, command = Show, indicatoron=0).grid(column=2, row=0)
Radiobutton(E12_frame, text="BOS",   variable=p['ShowTyp'], value=2, command = Show, indicatoron=0).grid(column=3, row=0)
Radiobutton(E12_frame, text="BOE",   variable=p['ShowTyp'], value=4, command = Show, indicatoron=0).grid(column=4, row=0)

Radiobutton(E2b_frame, text="All", variable=p['All99866'], value=0, command = Show, indicatoron=0).grid(column=0, row=0)
Radiobutton(E2b_frame, text="998", variable=p['All99866'], value=1, command = Show, indicatoron=0).grid(column=1, row=0)
Radiobutton(E2b_frame, text="66",  variable=p['All99866'], value=2, command = Show, indicatoron=0).grid(column=2, row=0)
Radiobutton(E2c_frame, text="All",   variable=p['LR'], value=0, command = Show, indicatoron=0).grid(column=0, row=0)
Radiobutton(E2c_frame, text="LL RR", variable=p['LR'], value=1, command = Show, indicatoron=0).grid(column=1, row=0)
Radiobutton(E2c_frame, text="LR RL", variable=p['LR'], value=2, command = Show, indicatoron=0).grid(column=2, row=0)
Radiobutton(E2d_frame, text="Square",   variable=p['SqTri'], value=0, command = Show, indicatoron=0).grid(column=0, row=0)
Radiobutton(E2d_frame, text="Triangle", variable=p['SqTri'], value=1, command = Show, indicatoron=0).grid(column=1, row=0)
Radiobutton(E2e_frame, text="AB", variable=p['AB'], value=0, command = Show, indicatoron=0).grid(column=0, row=0)
Radiobutton(E2e_frame, text="A",  variable=p['AB'], value=1, command = Show, indicatoron=0).grid(column=1, row=0)
Radiobutton(E2e_frame, text="B",  variable=p['AB'], value=2, command = Show, indicatoron=0).grid(column=2, row=0)
Radiobutton(E2f_frame, text="Normal", variable=p['log'], value=0, command = Show, indicatoron=0).grid(column=0, row=0)
Radiobutton(E2f_frame, text="Log10",  variable=p['log'], value=1, command = Show, indicatoron=0).grid(column=1, row=0)
Radiobutton(E2g_frame, text="A", variable=p['patient'], value=0, command = Show, indicatoron=0).grid(column=0, row=0)
Radiobutton(E2g_frame, text="1", variable=p['patient'], value=1, command = Show, indicatoron=0).grid(column=1, row=0)
Radiobutton(E2g_frame, text="2", variable=p['patient'], value=2, command = Show, indicatoron=0).grid(column=2, row=0)
Radiobutton(E2g_frame, text="3", variable=p['patient'], value=3, command = Show, indicatoron=0).grid(column=3, row=0)
Radiobutton(E2g_frame, text="4", variable=p['patient'], value=4, command = Show, indicatoron=0).grid(column=4, row=0)
Radiobutton(E2g_frame, text="5", variable=p['patient'], value=5, command = Show, indicatoron=0).grid(column=5, row=0)
Radiobutton(E2h_frame, text="FCE", variable=p['SC'], value=0, command = Show, indicatoron=0).grid(column=0, row=0)
Radiobutton(E2h_frame, text="SC",  variable=p['SC'], value=1, command = Show, indicatoron=0).grid(column=1, row=0)

Checkbutton(E3a_frame, text="diag",  variable=p['diag']   ).grid(column=0, row=0)
Checkbutton(E3a_frame, text="posit", variable=p['posit']  ).grid(column=1, row=0)
Checkbutton(E3a_frame, text="force", variable=p['force']  ).grid(column=2, row=0)
Checkbutton(E3a_frame, text="c_rev", variable=p['c_rev']  ).grid(column=0, row=1)
Checkbutton(E3a_frame, text="norm",  variable=p['norm']   ).grid(column=1, row=1)
Checkbutton(E3a_frame, text="AA",    variable=p['AA']     ).grid(column=2, row=1)
Checkbutton(E3a_frame, text="foFC",  variable=p['forceFC']).grid(column=1, row=2)
Radiobutton(E3b_frame, text="Interp:None", variable=p['interp'], value='bicubic', command = Show, indicatoron=0).grid(column=0, row=0)
Radiobutton(E3b_frame, text="Interp:Near", variable=p['interp'], value='nearest', command = Show, indicatoron=0).grid(column=1, row=0)
textentry(E3c_frame, p['min'], 'min', 7)
textentry(E3c_frame, p['max'], 'max', 7)
textentry(E3c_frame, p['seuil'], 'log seuil', 7)
listc= Listbox(E3d_frame, width=20, height=7, exportselection=0); listc.pack(side='top')
for nW in ['Blues','BuGn','BuPu','GnBu','Greens','Greys','Oranges','OrRd','PuBu','PuBuGn','PuRd','Purples','RdPu','Reds','YlGn','YlGnBu','YlOrBr','YlOrRd','',\
           'gist_earth','gist_gray','gist_heat','gist_ncar','gist_rainbow','gist_stern','gist_yarg','',\
           'autumn','bone','cool','copper','flag','gray','hot','hsv','jet','pink','prism','spectral','spring','summer','winter','',\
           'BrBG','PiYG','PRGn','PuOr','RdBu','RdGy','RdYlBu','RdYlGn','Spectral']:
    listc.insert(END, '%s' %nW)
listc.bind('<Double-1>', listcDef)

Button(E_frame, text='Exit', command=Quit).pack(side='bottom')
root.mainloop()
root.destroy()
