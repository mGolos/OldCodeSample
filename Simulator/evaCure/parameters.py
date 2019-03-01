#!/usr/bin/env python
#-*- coding:Utf-8 -*-

''' @author: mathieu.golos@gmail.com
'''

from attribute import attributeClass


class evaPar(attributeClass):
    '''All parameters needed to run evaCure:
    main:
    - dur     : simulation lenght (steps)
    - rperiod : step period for monitoring
    - out     : monitoring variables
    - init    : 'rand', 'patt', 'ext', 'file' initialization
    - dens    : density for random initialization
    - normX   : norm of the initial vector
    - stdD_0  : initial noise
    - patt    : pattern number for initialization
    - m       : number of patterns
    
    connectivity:
    - connAd   : adress of the structural connectivity matrix
    - distAd   : adress of the nodes distances matrix
    - vel      : defines the velocity for delays (0 disable delays)
    - tauC     : time constant for plasticity
    - normC    : defines the renormalization function (<0 disable it)
    - normT    : defines the normalization type (Frobenius by default)
    
    noise:
    - colors  : types of noise ('white', 'pink', 'brown')
    - stdD_x  : standard deviation for potentials
    - stdD_T  : standard deviation for inhibition
    - AonF    : parameter for brownian noise (1: pink noise)
        
    model:
    - model     : model
    - threshold : type of threshold ('global', 'local')
    - dt        : timestep
    - taux      : time constant for potentials
    - tauT      : time constant for global inhibition theta
    - P         : Excitatory factor on Inhibitory factor Psi/Omega
    - G         : Inhibitory factor dot thresholding function factor Omega.gamma
    - x         : Potential value for each nodes
    - theta     : inhibition value
    - A         : activity of neurons for each nodes
    
    other:
    - argv : use external parameters from sys.argv

    Functions:
    - __init__      
    - fromOutside
    - update ?.
    '''
    keys = ['dur', 'rperiod', 'out', 'init', 'dens', 'normX', 'stdD_0', 'patt', 'm', 
            'connAd', 'distAd', 'vel', 'tauC', 'normC', 'normT',
            'colors', 'stdD_x', 'stdD_T', 'AonF',
            'model', 'threshold', 'dt', 'taux', 'tauT', 'P', 'G', 
            'x', 'A', 'theta',
            'argv']

    # main
    dur = 10000
    rperiod = 10
    out = {'rtime':[], 'x': []}
    init = 'rand'
    dens = 0.5
    normX = None
    stdD_0 = 0.
    patt = 0
    m = 0
    
    # connectivity
    connAd = 'Connectomes/SC_H_998_0.npy'
    distAd = 'Connectomes/FL_H_998_0.npy'
    vel = 0.
    tauC = 0.
    normC = 1.
    normT = 'fro'
    
    # noise
    colors = 'white'     
    stdD_x = 0.
    stdD_T = 0.
    AonF = 1
    
    # model
    model = 'HopfieldBasedDynamic' 
    threshold = 'global'
    dt = 0.1
    taux = 10.
    tauT = 50.
    P = 1.
    G = 60.
    x = None
    A = None
    theta = 0.5
    
    #other
    argv = False

    
    def __init__(self, **kw):
        '''Construction from attributeClass class.'''
        attributeClass.__init__(self, **kw)
        if self.argv:
            self.fromOutside()
            
    
    def fromOutside(self):
        int_var = ['dur', 'rtime', 'rperiod', 'patt', 'm', 'AonF']   
        dbl_var = ['dens', 'normX', 'stdD_0', 'vel', 'tauC', 'normC', 'stdD_x', 'stdD_T', 
                   'dt', 'taux', 'tauT', 'P', 'G', 'theta']
        str_var = ['init', 'connAd', 'distAd', 'model', 'threshold']
        lst_var = ['out', 'colors']
        
        from sys import argv as sys_argv
        for i in range(1, len(sys_argv), 2):
            try:
                option = sys_argv[i]
                if   option in int_var:
                    setattr(self, option, int(sys_argv[i+1]))
                    
                elif option in dbl_var:
                    setattr(self, option, float(sys_argv[i+1]))
                    
                elif option in str_var:
                    setattr(self, option, str(sys_argv[i+1]))
                    
                elif option in lst_var:
                    # list of var format have to be splitted by '|' like : 't|x|A|theta'
                    list_var = str(sys_argv[i+1]).split('|')                    
                    setattr(self, option, list_var)
                    
                else:
                    print 'Options invalides :',option,'->',sys_argv[i+1]
            except:
                pass


#    #def Take_Tk():
#    #TODO
#        #''' Initialization of the Tkinter window in GUI.py
#        #'''
#        #global root
#        #root = Tk()
#        #root.title("evaCure")
#
#    #def TakeParameters_GUI():
#    #TODO
#        #''' Send parameters in a library for Tkinter button and all in GUI.py
#        #'''
#        #global p
#        #p = {}
#        #p['N']          = IntVar();     p['N'].set(N)
#        #p['m']          = IntVar();     p['m'].set(m)
#        #p['priT']       = IntVar();     p['priT'].set(priT)
#        #p['nper']       = IntVar();     p['nper'].set(nper)
#        #p['GUI']        = IntVar();     p['GUI'].set(GUI)
#        #p['netw']       = IntVar();     p['netw'].set(netw)
#        #p['patt']       = IntVar();     p['patt'].set(patt)
#        #p['local']      = IntVar();     p['local'].set(local)
#        #p['init']       = IntVar();     p['init'].set(init)
#        #p['who']        = IntVar();     p['who'].set(who)
#        #p['TtS']        = IntVar();     p['TtS'].set(TtS)
#        #p['DSTC']       = IntVar();     p['DSTC'].set(DSTC)
#        #p['dt']         = DoubleVar();  p['dt'].set(dt)
#        #p['dur']        = DoubleVar();  p['dur'].set(dur)
#        #p['theta']      = DoubleVar();  p['theta'].set(theta)
#        #p['err']        = DoubleVar();  p['err'].set(err)
#        #p['dens']       = DoubleVar();  p['dens'].set(dens)
#        #p['sigma20']    = DoubleVar();  p['sigma20'].set(sigma20)
#        #p['sigma2x']    = DoubleVar();  p['sigma2x'].set(sigma2x)
#        #p['sigma2T']    = DoubleVar();  p['sigma2T'].set(sigma2T)
#        #p['sigma2A']    = DoubleVar();  p['sigma2A'].set(sigma2A)
#        #p['P']          = DoubleVar();  p['P'].set(P)
#        #p['G']          = DoubleVar();  p['G'].set(G)
#        #p['taux']       = DoubleVar();  p['taux'].set(taux)
#        #p['tauT']       = DoubleVar();  p['tauT'].set(tauT)
#        #p['tauW']       = DoubleVar();  p['tauW'].set(tauW)
#        #p['fileIN']     = StringVar();  p['fileIN'].set(fileIN)
#        #p['direIN']     = StringVar();  p['direIN'].set(direIN)
#        #p['Cnnctm']     = StringVar();  p['Cnnctm'].set(Cnnctm)
#        #p['TCOUT']      = StringVar();  p['TCOUT'].set(TCOUT)
#
#    #def Give_InitialParameters_GUI():
#    #TODO
#        #''' Send parameters from the library used for Tkinter in GUI.py to variables used in Main.py
#        #'''
#        #global N,m,dt,dur,theta,nper,err,GUI,test,netw,patt,local,init,who
#        #global fileIN,direIN,Cnnctm,dens,sigma20,sigma2x,sigma2T,sigma2A,P,G
#        #global tauT,tauW,priT,taux,DSTC,TtS,TCOUT
#
#        #N         = p['N'].get()
#        #m         = p['m'].get()
#        #dt        = p['dt'].get()
#        #dur       = p['dur'].get()
#        #theta     = p['theta'].get()
#        #nper      = p['nper'].get()
#        #err       = p['err'].get()
#        #GUI       = p['GUI'].get()
#        #netw      = p['netw'].get()
#        #patt      = p['patt'].get()
#        #local     = p['local'].get()
#        #init      = p['init'].get()
#        #who       = p['who'].get()
#        #fileIN    = p['fileIN'].get()
#        #direIN    = p['direIN'].get()
#        #Cnnctm    = p['Cnnctm'].get()
#        #dens      = p['dens'].get()
#        #sigma20   = p['sigma20'].get()
#        #sigma2x   = p['sigma2x'].get()
#        #sigma2T   = p['sigma2T'].get()
#        #sigma2A   = p['sigma2A'].get()
#        #P         = p['P'].get()
#        #G         = p['G'].get()
#        #tauT      = p['tauT'].get()
#        #tauW      = p['tauW'].get()
#        #priT      = p['priT'].get()
#        #taux      = p['taux'].get()
#        #DSTC      = p['DSTC'].get()
#        #TtS       = p['TtS'].get()
#        #TCOUT     = p['TCOUT'].get()
#
#    #def Give_RunParameters_GUI(Network, Itechek):
#    #TODO
#        #''' Permit some parameters to be updated in GUI.py
#        #'''
#        #global p, priT, dur, dt
#        #priT            = p['priT'].get()
#        #dur             = p['dur'].get()
#        #dt              = p['dt'].get()
#        #Itechek.nmax    = int(p['dur'].get() / p['dt'].get())
#        #Network.dt      = p['dt'].get()
#        #Network.dur     = p['dur'].get()
#        #Network.err     = p['err'].get()
#        #Network.taux    = p['taux'].get()
#        #Network.tauT    = p['tauT'].get()
#        #Network.tauW    = p['tauW'].get()
#        #Network.P       = p['P'].get()
#        #Network.G       = p['G'].get()
#        #Network.sigma2x = p['sigma2x'].get()
#        #Network.sigma2T = p['sigma2T'].get()
#        #Network.sigma2A = p['sigma2A'].get()
#        #Network.sigma20 = p['sigma20'].get()
#        #TCOUT           = p['TCOUT'].get()
#        #DSTC            = p['DSTC'].get()
