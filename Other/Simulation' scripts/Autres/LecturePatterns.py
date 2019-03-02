#!/usr/bin/env python
#-*- coding:Utf-8 -*-
''' Affiche les patterns en fonction de gamma, normW et la densité
@author: gmoaltos@hotmail.com
'''

from os import mkdir, path, listdir
from pylab import *
from numpy import load as LOAD
from random import randrange
from Tkinter import*
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

def textentry(root, variable, label):
    ''' Description '''
    global row_counter
    l= Label(root, text=label)
    l.grid(column=0, row=row_counter, sticky='w')
    widget = Entry(root, textvariable=variable, width=10)
    widget.grid(column=1, row=row_counter)
    row_counter += 1
    return widget

def Quit():
    ''' Description '''
    root.quit()
    W_frame.destroy()

def listaDef(e):
    global p
    p['NormW'].set(float(lista.get(lista.curselection())))
    letMeSee()

def listbDef(e):
    global p
    p['dens_moy'].set(float(listb.get(listb.curselection())))
    letMeSee()

def letMeSee():
    global p
    dir_pri = p['dir_pri'].get()
    dir_sub = dir_pri + '/Norm_%.2f' %(p['NormW'].get())

    '''FIGURE 1'''
    W_f_norm = zeros(len(listdir(dir_pri)) - 1)
    for i in range(len(W_f_norm)):
        tmp = listdir(dir_pri)
        tmp.sort()
        W_f_norm[i] = (len(listdir(dir_pri + '/' + tmp[i+1])) -1)/2

    fig1.clf()
    b= fig1.add_subplot(111)
    b.plot(W_f_norm)
    canvas1.show()

    '''FIGURE 2'''
    nb_patt = len(listdir(dir_sub)) / 2
    fig2.clf()
    for i in range(nb_patt):
        patt = LOAD(dir_sub + '/pattern_%.3i.npy' %i)
        if p['N'].get()==47 and p['to66'].get():
            patt = (patt.reshape(1,len(patt))).tolist()[0]
            aRemettre = [0,1,2,3,5,11,21,22,25,27,32,40,55,60,61,62,63,64,65]
            for j in aRemettre:
                patt.insert(j,0)

        a = fig2.add_subplot(nb_patt,1,i+1)
        a.imshow([patt], interpolation="nearest", cmap=get_cmap("bone"), vmin=0, vmax=1)
        if not p['axis'].get():
            a.axis('off')
    canvas2.show()


fig1 = Figure(figsize=(10,4), dpi=50)
fig2 = Figure(figsize=(10,6), dpi=100)

root = Tk()
root.title("Lecture Patterns")
W_frame  = Frame(root, height=20);              W_frame.pack(side='left', pady=5, padx=8)
E_frame = LabelFrame(root);                     E_frame.pack(side='right', pady=5, padx=8)
EN_frame = LabelFrame(E_frame, text= "Model");  EN_frame.pack(side='top', pady=5, padx=8)
ES_frame = LabelFrame(E_frame);                 ES_frame.pack(side='top', pady=5, padx=8)

''' LEFT '''
canvas1 = FigureCanvasTkAgg(fig1, master=W_frame)
canvas1.show()
canvas1.get_tk_widget().pack(side="top", fill=BOTH, expand=1)
canvas1._tkcanvas.pack(side="top", fill=BOTH, expand=1)
canvas2 = FigureCanvasTkAgg(fig2, master=W_frame)
canvas2.show()
canvas2.get_tk_widget().pack(side="top", fill=BOTH, expand=1)
canvas2._tkcanvas.pack(side="top", fill=BOTH, expand=1)

''' RIGHT '''
row_counter= 0
p = {}
p['dir_pri'] = StringVar();         p['dir_pri'].set("../DtiCluster/Premieres/VarNorm_Discret_66_theta fix")
p['NormW'] = DoubleVar();           p['NormW'].set(0.5)
p['dens_moy'] = DoubleVar();        p['dens_moy'].set(0.02)
p['gamma'] = IntVar();              p['gamma'].set(10)
p['N'] = IntVar();                  p['N'].set(66)
p['to66'] = IntVar();               p['to66'].set(0)
p['axis'] = IntVar();               p['axis'].set(0)
p['theta_dyn'] = IntVar();          p['theta_dyn'].set(0)

textentry(EN_frame, p['N'], 'N')
textentry(EN_frame, p['gamma'], 'gamma')

Entry(E_frame, textvariable=p['dir_pri']).pack()
Checkbutton(E_frame, text="47 to 66", variable=p['to66']).pack()
Checkbutton(E_frame, text="axis", variable=p['axis']).pack()
#Checkbutton(E_frame, text="theta dynamique", variable=p['theta_dyn']).pack()


lista= Listbox(ES_frame, width=10, height=33, exportselection=0);   lista.pack(side='left')
listb= Listbox(ES_frame, width=10, height=33, exportselection=0);   listb.pack(side='left')
for nW in arange(1.4,2.601,0.05):
    lista.insert(END, '%.2f' %nW)
for nW in arange(0.02,0.981,0.03):
    listb.insert(END, '%.2f' %nW)
lista.bind('<Double-1>', listaDef)
listb.bind('<Double-1>', listbDef)

Button(E_frame, text='Exit', command=Quit).pack(side='bottom')

root.mainloop()
root.destroy()

