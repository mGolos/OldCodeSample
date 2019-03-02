#!/usr/bin/env python
#-*- coding:Utf-8 -*-
''' Affiche les patterns en fonction de gamma, Psi et la densité
@author: gmoaltos@hotmail.com
'''

import os
from pylab import *
from numpy import load as LOAD
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from Tkinter import *


def loadTendance(Folder):
    ofi = open(Folder + '/tendance.txt', 'r')
    tendance = array([line.split() for line in ofi]).astype(float)
    ofi.close()
    return tendance

def textentry(root, variable, label):
    global row_counter
    l= Label(root, text=label)
    l.grid(column=0, row=row_counter, sticky='w')
    widget = Entry(root, textvariable=variable, width=10)
    widget.grid(column=1, row=row_counter)
    row_counter += 1
    return widget

def Quit():
    root.quit()
    W_frame.destroy()

def Direct():
    if p['PatTen'].get():
        letMeSee( clustering( giveMeMyPatterns() ) )
    else:
        if p['rangeD'].get():
            letMeTendanceAll()
        else:
            letMeTendance()

def listaDef(e):
    p['NormW'].set(float(lista.get(lista.curselection())))
    Direct()


def giveMeMyPatterns():
    global nombre_patterns
    S = []
    NormW = p['NormW'].get()
    dir_pri = p['dir_pri'].get()
    rangeD = p['rangeD'].get()
    distrib_patt = []

    for NormW_i in arange(NormW - rangeD, NormW + rangeD + 0.001, 0.05):
        dir_sub = dir_pri + '/Norm_%.2f' %(NormW_i)
        if os.path.exists(dir_sub):
            nb_patt = len(os.listdir(dir_sub)) / 2
            for i in range(nb_patt):
                S.append([ LOAD(dir_sub + '/pattern_%.3i.npy' %i) ])
            distrib_patt.append(nb_patt)
        else:
            print "dir_sub doesn't exists"
    print distrib_patt
    nombre_patterns = len(S)
    return S

def similarity(E,F):
    sim = []
    for i in range(len(E)):
        for j in range(len(F)):
            sim.append(corrcoef(E[i], F[j])[0,1])
    return where(p['mmm'].get() == 0, min(sim), \
           where(p['mmm'].get() == 1, average(sim), \
                                      max(sim)))

def clustering(S):
    sim_final = 1
    while len(S) > 1 and sim_final > p['sim'].get():
        sim = -1
        for i in range(len(S)):
            for j in range(i + 1, len(S)):
                simEF = similarity(S[i],S[j])
                if simEF > sim:
                    indE = i
                    indF = j
                    sim = simEF
        sim_final = sim
        for i in range(len(S[indF])):
            S[indE].append(S[indF].pop(0))
        S.remove([])
    return S

def letMeSee(S):
    fig.clf()
    indice = 1
    abreviations = ['ENT','PARH','TP','FP','FUS','TT','LOCC','SP','IT','IP',\
                    'SMAR','BTST','MT','ST','PSTC','PREC','CMF','POPE','PTRI','RMF',\
                    'PORB','LOF','CAC','RAC','SF','MOF','LING','PCAL','CUN','PARC',\
                    'ISTC','PCUN','PC']
    abre_ind = []
    for i in range(len(abreviations)): abre_ind.append('r'+ abreviations[i])
    abreviations.reverse()
    for i in range(len(abreviations)): abre_ind.append('l'+ abreviations[i])

    for i in range(len(S)):
        for j in range(len(S[i])):

            if not p['axis'].get():
                a = fig.add_subplot(nombre_patterns, 1, indice)
                if p['max'].get():
                    a.imshow([S[i][j]], interpolation="nearest", cmap=get_cmap(p['cmap'].get()), vmin=0, vmax=p['max'].get())
                else:
                    a.imshow([S[i][j]], interpolation="nearest", cmap=get_cmap(p['cmap'].get()))
                a.axis('off')
            else:
                a = fig.add_subplot(1, nombre_patterns, indice)
                if p['max'].get():
                    a.imshow(array([S[i][j]]).T, interpolation="nearest", cmap=get_cmap(p['cmap'].get()), vmin=0, vmax=p['max'].get())
                else:
                    a.imshow(array([S[i][j]]).T, interpolation="nearest", cmap=get_cmap(p['cmap'].get()))
                if nombre_patterns < 9 or indice == 1:
                    a.set_yticks(range(66))
                    a.set_yticklabels(abre_ind, fontsize=9)
                else:
                    a.yaxis.set_visible(False)
                    a.axis('off')
                a.xaxis.set_visible(False)
            if j == 0:
                a.set_title(len(S[i]))
            indice += 1
    canvas.show()

def letMeTendance():
    fig.clf()
    tendance = loadTendance(p['dir_pri'].get() + '/Norm_%.2f' %(p['NormW'].get()))
    a = fig.add_subplot(111)

    width = 1./ len(tendance)
    color = linspace(0,1,len(tendance[0])-1)

    a.bar(tendance[:,0], tendance[:,1], width, color=str(color[0]))
    soum=tendance[:,1]
    for i in range(2,len(tendance[0])) :
        a.bar(tendance[:,0], tendance[:,i], width, bottom=soum, color=str(color[i-1]))
        soum = soum + tendance[:,i]

    a.set_xlim(0,1)
    a.set_xlabel(r'$Initial\ density$', fontsize=15)
    a.set_ylabel(r'$Tendance$', fontsize=15)
    #print [sum(tendance[k,1:]) for k in range(len(tendance))]
    canvas.show()

def letMeTendanceAll():
    fig.clf()
    dir_pri = p['dir_pri'].get()
    rangefory0 = arange(1.4, 2.6 + 0.001, 0.05)

    tab = []
    for NormW_i in rangefory0:
        dir_sub = dir_pri + '/Norm_%.2f' %(NormW_i)
        if os.path.exists(dir_sub):
            tendance = loadTendance(dir_sub)
            tendance[:,1] = [sum(array(tendance[k,1:])!=0) for k in range(len(tendance))]
            tab.append(tendance[:,1])
        else:
            print "dir_sub doesn't exists"
    Tendances = array(tab)

    fig.add_subplot(111)
    x0 = tendance[:,0]
    y0 = rangefory0
    contourf(y0, x0, Tendances.T, int(floor(Tendances.max())), cmap=get_cmap(p['cmap'].get()))
    #imshow(Tendances.T, interpolation='nearest', cmap=get_cmap(p['cmap'].get()))
    colorbar()
    ylabel(r'$Initial\ density$', fontsize=20)
    xlabel(r'$||W||$', fontsize=20)
    canvas.show()


root = Tk()
root.title("Clustering_0")
W_frame   = Frame(root, height=20);             W_frame.pack(side='left', pady=5, padx=8)
E_frame   = LabelFrame(root, bd=0);             E_frame.pack(side='right', pady=5, padx=8)
EN_frame  = LabelFrame(E_frame, text= "Model"); EN_frame.pack(side='top', pady=5, padx=8)
ES_frame  = LabelFrame(E_frame, bd=0);          ES_frame.pack(side='top', pady=5, padx=8)
ES1_frame = LabelFrame(ES_frame, bd=0);         ES1_frame.pack(side='top', pady=5, padx=8)
ES2_frame = LabelFrame(ES_frame, bd=0);         ES2_frame.pack(side='top', pady=5, padx=8)
ES4_frame = LabelFrame(ES_frame, bd=0);         ES4_frame.pack(side='top', pady=5, padx=8)
ES3_frame = LabelFrame(ES_frame, bd=0);         ES3_frame.pack(side='top', pady=5, padx=8)


''' LEFT '''
fig = figure(figsize=(10,8), dpi=100)
canvas = FigureCanvasTkAgg(fig, master=W_frame)
canvas.show()
canvas.get_tk_widget().pack(side="top", fill=BOTH, expand=1)
canvas._tkcanvas.pack(side="top", fill=BOTH, expand=1)


''' RIGHT '''
row_counter= 0
p = {}
p['dir_pri'] = StringVar();         p['dir_pri'].set("../DtiCluster/Premieres/VarNorm_Discret_66_theta_fix")
p['NormW'] = DoubleVar();           p['NormW'].set(0.5)
p['rangeD'] = DoubleVar();          p['rangeD'].set(0.00)
p['algo'] = IntVar();               p['algo'].set(1)
p['sim'] = DoubleVar();             p['sim'].set(0.7)
p['mmm'] = IntVar();                p['mmm'].set(0)
p['cmap'] = StringVar();            p['cmap'].set("bone")
p['max'] = DoubleVar();             p['max'].set(0)
p['PatTen'] = IntVar();             p['PatTen'].set(0)

p['axis'] = IntVar();               p['axis'].set(0)
p['theta_dyn'] = IntVar();          p['theta_dyn'].set(0)

textentry(EN_frame, p['rangeD'], 'range_dens_moy')
textentry(EN_frame, p['sim'], 'Sim_coef')
textentry(EN_frame, p['dir_pri'], 'dir_pri')
textentry(EN_frame, p['max'], 'max')
textentry(EN_frame, p['cmap'], 'cmap')

Checkbutton(E_frame, text="axis + switch", variable=p['axis']).pack()
Button(E_frame, text='Exit', command=Quit).pack(side='bottom')

Radiobutton(ES1_frame, text="Dis Tfix", variable=p['dir_pri'], value="../DtiCluster/Premieres/VarNorm_Discret_66_theta_fix", command = Direct, indicatoron=0).grid(row=0, column=0)
Radiobutton(ES1_frame, text="Dis Tvar", variable=p['dir_pri'], value="../DtiCluster/Premieres/VarNorm_Discret_66_theta_var", command = Direct, indicatoron=0).grid(row=1, column=0)
Radiobutton(ES1_frame, text="G20 Tfix", variable=p['dir_pri'], value="../DtiCluster/Premieres/VarNorm_Gamma_20_ss_vartheta", command = Direct, indicatoron=0).grid(row=0, column=1)
Radiobutton(ES1_frame, text="G30 Tvar", variable=p['dir_pri'], value="../DtiCluster/Premieres/VarNorm_Gamma_30", command = Direct, indicatoron=0).grid(row=1, column=1)
Radiobutton(ES1_frame, text="G30 Tfix", variable=p['dir_pri'], value="../DtiCluster/Premieres/VarNorm_Gamma_30_ss_vartheta", command = Direct, indicatoron=0).grid(row=0, column=2)
Radiobutton(ES2_frame, text="Min", variable=p['mmm'], value="0", command = Direct, indicatoron=0).pack(side='left')
Radiobutton(ES2_frame, text="Moy", variable=p['mmm'], value="1", command = Direct, indicatoron=0).pack(side='left')
Radiobutton(ES2_frame, text="Max", variable=p['mmm'], value="2", command = Direct, indicatoron=0).pack(side='left')
Radiobutton(ES4_frame, text="Tendance", variable=p['PatTen'], value="0", command = Direct, indicatoron=0).pack(side='left')
Radiobutton(ES4_frame, text="Patterns", variable=p['PatTen'], value="1", command = Direct, indicatoron=0).pack(side='left')

lista= Listbox(ES3_frame, width=10, height=33, exportselection=0);   lista.pack(side='left')
for nW in arange(1.4,2.601,0.05):
    lista.insert(END, '%.2f' %nW)
lista.bind('<Double-1>', listaDef)


root.mainloop()
root.destroy()