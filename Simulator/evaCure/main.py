#!/usr/bin/env python
#-*- coding:Utf-8 -*-

'''@author: mathieu.golos@gmail.com
'''

from attribute import attributeClass
from tools import loadPatterns, array2data, estimAndPercent, tic
from time import time
from pylab import zeros, copy
from numpy import ndarray
try:
    from pyprind import ProgBar
except:
    pass

class evaCure(attributeClass):
    '''Construction from the basic attributeClass.
    Parameters:
    - dur     : simulation lenght (steps)
    - rperiod : step period for monitoring
    - out     : monitoring variables
    - init    : 'rand', 'patt', 'ext', 'file' initialization
    - dens    : density for random initialization
    - normX   : norm of the initial vector
    - stdD_0  : initial noise
    - patt    : pattern number for initialization
    - m       : nuber of patterns
    - kwa     : any other parameters is favored
    
    Build or Given:
    - evaPar : Parameters class
    - evaCon : Connectivity class
    - evaNoi : Noise class
    - evaMod : Model class
    
    Functions:
    - __init__
    - getParameters
    - getConnectivity
    - getModel
    - getPatterns
    - initialization
    - defineOut
    - updateOut
    - update
    - updateTillEnd
    - toEquilibrium.
    '''
    keys = ['dur', 'rtime', 'rperiod', 'out', 'init', 'dens', 'normX', 'stdD_0', 'patt', 'm']
           
    dur = 10000
    rtime = 0
    rperiod = 10
    out = {'rtime':[], 'x': []}
    init = 'rand'
    dens = 0.5
    normX = None
    stdD_0 = 0.
    patt = 0
    m = 0
    iFile = 0
    

    def __init__(self, **kwa):
        '''Initialize the system depending on the different classes or dictionnaries inputed.'''
        attributeClass.__init__(self, **kwa)
        
        self.setAttr(**kwa)
        kwa = self.getParameters(**kwa)
        self.setAttr(**kwa)
        kwa = self.getConnectivity(**kwa)
        kwa = self.getModel(**kwa)
        kwa = self.getNoise(**kwa)
        #self.keys.extend(['evaPar', 'evaCon', 'evaNoi', 'evaMod'])
        self.evaPar.setAttr(**kwa)
        
        self.getPatterns(**kwa)
        self.initialization(**kwa)
        self.initialization(**kwa)
        self.defineOut(**kwa)
        
        
    def getParameters(self, evaPar=None, **kwa):
        if evaPar is None:
            from evaCure import parameters
            self.evaPar = parameters.evaPar(**kwa)
        elif evaPar.__class__ == dict:
            from evaCure import parameters
            self.evaPar = parameters.evaPar(**dict(kwa, **evaPar))
        else:
            self.evaPar = evaPar
            
        kwa.update(**self.evaPar.getAttr())
        return kwa

            
    def getConnectivity(self, evaCon=None, **kwa):
        if evaCon is None:
            from evaCure import connectivity
            self.evaCon = connectivity.evaCon(**kwa)
        elif evaCon.__class__ == dict:
            from evaCure import connectivity
            self.evaCon = connectivity.evaCon(**dict(kwa, **evaCon))
        else:
            self.evaCon = evaCon
            
        kwa.update(self.evaCon.getAttr())
        return kwa
           
           
    def getNoise(self, evaNoi=None, **kwa):
        if evaNoi is None:
            from evaCure import noise
            self.evaNoi = noise.evaNoi(**kwa)
        elif evaNoi.__class__ == dict:
            from evaCure import noise
            self.evaNoi = noise.evaNoi(**dict(kwa, **evaNoi))
        else:
            self.evaNoi = evaNoi
            
        kwa.update(self.evaNoi.getAttr())
        return kwa
           
           
    def getModel(self, evaMod=None, **kwa):
        if evaMod is None:
            from evaCure import model
            self.evaMod = model.choice(**kwa)
        elif evaMod.__class__ == dict:
            from evaCure import model
            self.evaMod = model.choice(**dict(kwa, **evaMod))
        else:
            self.evaMod = evaMod
            
        kwa.update(self.evaMod.getAttr())
        return kwa


    def getPatterns(self, mmax=None, **kwa):
        ''' Define A and x patterns from *.npy files in PatternsIN folder.
        The value mmax permits to load only the mmax first patterns.'''
        if self.m:
            self.ksix = loadPatterns(self.pattD.replace('_S','_X'), mmax=self.m)
            self.ksiA = loadPatterns(self.pattD.replace('_X','_S'), mmax=self.m)
            self.m = self.ksix.shape[0]
           
            
    def initialization(self, **kwa):
        '''Diverse initializations.'''
        if self.init == 'rand':
            self.evaMod.initialRand(self.dens, self.normX)
            
        elif self.init == 'patt':
            if self.patt < self.m:
                self.evaMod.initialPatt(self.ksix[self.patt],
                                        self.ksiA[self.patt])
            elif self.m != 0:
                print 'There is no pattern %i -> 1st choosen on %i patterns' %(self.patt, self.m)
                self.patt = 0
                self.initialization(**kwa)
            else:
                print 'Not any patterns loaded -> random initialization.'
                self.init = 'rand'
                self.initialization(**kwa)
                
        elif self.init == 'ext':
            self.evaMod.initialExt()
            
        elif self.init == 'file':
            self.evaMod.initialFile()
            
        if self.stdD_0 != None:
            from tools import whiteNoise
            noise_0 = whiteNoise(dim1=self.evaCon.N, stdD=self.stdD_0)
            self.evaMod.initialNoise(noise_0)
            
        self.evaMod.initTheta()
        self.iteCheck = IterationChecker(**kwa)


    def defineOut(self, out=None, **kwa):
        '''Define witch variables will be followed and stocked.
        Allowed variables: x, A, ...
        '''
        if out != None:
            self.out = {}
            for p in out:
                self.out[p] = []
        self.updateOut()
        

    def updateOut(self):
        for p in self.out.keys():
            if p in self.keys:
                tmp = getattr(self, p)
            elif p in self.evaMod.keys:
                tmp = getattr(self.evaMod, p)
            elif p in self.evaCon.keys:
                tmp = getattr(self.evaCon, p)
            else:
                print 'Error in updateOut with %s.' %p
                pass
                
            self.out[p].append(copy(tmp))
            #if type(tmp) == ndarray:
                #self.out[p].append(tmp.copy())
            #else:
                #self.out[p].append(tmp)
        

    def update(self):
        '''Update the system for a timestep.'''
        self.evaMod.update( self.evaNoi.update() )
        self.evaCon.update()
        self.rtime += 1
        if not self.rtime % self.rperiod:
            self.updateOut()
            
            
    def updateTillEnd(self, estimOver=1000, t0=None, Tsave=None, Fname=None):
        if t0 is None:
            t0 = tic(ret=True)
        
        if type(self.dur) == str:
            self.dur = timeToSteps(self.dur, self.evaMod.dt)     
            
        for i in xrange(estimOver):
            self.update()
            estimAndPercent(self.rtime, self.dur, avg=estimOver, t0=t0)           
               
        if self.rtime == estimOver:
            try:    bar = ProgBar(self.dur - estimOver, stream=1)
            except: pass
            
        for i in xrange(estimOver, self.dur):
            self.update()
            try:    bar.update()
            except: pass
            
            if (Tsave != None) and (Fname != None):
                if not self.rtime % Tsave:
                    array2data(self.out, Fname + '%i.npy'%self.iFile)
                    self.iFile += 1
                    for k in self.out.keys():
                        del self.out[k][:]
       
        
    def toEquilibrium(self):
        '''Run the simulation till it takes too much time or 
        if the system become stationnary.
        '''
        while self.iteCheck.dontstop(self.evaMod.x.mean()):
            self.update()
            

    #def printStop(self):
    #TODO
        #if Network.m:
            #Network.closers()
            #print "dHamming at start       :", Network.dHamm(S0)
            #print "dHamming at end         :", Network.dHamm(Network.x)
            #print "dCorrelation at start   :", ["%.1f"%k for k in Network.correlation(S0)]
            #print "dCorrelation at end     :", ["%.1f"%k for k in Network.correlation(Network.x)]
            #print "Closer patterns         :", Network.argpos
            #print "dCorrelation from them  :", Network.dpos

        #symmetry = ((Network.W.T - Network.W).todense() > 1e-14).sum() * 100 / float(N * N)
        #print "Theta                   :", Network.theta
        #print "Norm of connectivity    :", norm(Network.W.todense())
        #print "Density of connectivity :", Network.connectivityDensity()
        #print "Symmetry and minima     :", where(symmetry == 0, "Symmetric", "Asymmetric("+str(int(symmetry))+"%),")
        #print norm(Network.A), norm(Network.W.dot(Network.A))
        
        
        
    #def closers(self):
    #TODO
        #''' Calculate the position and the distance for the minimum Hamming distance
        #between neurons state and patterns ksi.
        #'''
        #D = self.correlation(self.x)
        #self.argpos = []
        #for i in xrange(self.m):
            #if D[i] == min(D):
                #self.argpos.append(i)
        #self.dpos = min(D)

    #def correlation(self, S):
    #TODO
        #''' Return the distance vector, between neurons state and patterns ksi,
        #extracted from the correlation matrix.
        #'''
        #C = corrcoef(S, self.ksix)[0,1:]
        #C = where(isnan(C), 0., C)
        #return where(C < 0., 0., C)

    #def dHamm(self, S):
    #TODO
        #''' Return the Hamming distance [m] between state S [N] and patterns ksi [N][m].'''
        #d = zeros(self.m)
        #for k in range(self.m):
            #d[k] = abs(S - self.ksix[k]).sum() / float(self.N)
        #return d


class IterationChecker(attributeClass):
    '''Class that stop iterations if the simulation takes too much time
    or if the system's state become stationnary during some time :
    Parameters:
    - nmax   : number of maximum iterations
    - nper   : number of iteration to consider a stationnary system (opt)
    
    Builded:
    - ncount : counter for iterations
    - system : past values of the system [nper]  (opt).
    
    Functions:
    .. __init__
    .. dontstop.
    '''
    keys = ['nmax', 'nper', 'err']
    nmax = 100000
    nper = 1000
    err  = 1e-6
    ncount = 0  
    
    def __init__(self, **kwa):
        '''Initialize the system depending on the different classes or dictionnaries inputed.'''
        attributeClass.__init__(self, **kwa)
        self.setAttr(**kwa)
        if self.nper:
            self.system = zeros(self.nper)


    def dontstop(self, variable):
        '''Always stop for ncount<nmax.
        Do not stop for the first nper iteration if nper<nmax.
        After nper, stop if the relative difference between average
        of 'system' and the last variable is lower than err.
        '''
        self.ncount += 1

        if self.ncount >= self.nmax:
            return False
        elif self.nper:
            self.system[self.ncount % self.nper] = variable
            if self.ncount < self.nper:
                return True
            else:
                sMean = self.system.mean()
                if abs((sMean - variable) / sMean) < self.err \
                or self.ncount >= self.nmax \
                or sMean < 1e-42:
                    return False
                else:
                    return True
        else:
            return True
           
           
def timeToSteps(dur, dt):
    '''Return the maximal number of step from a string value of time as (XXmXXs) and depending of the timestep.
    '''
    ste = 0
    if dur.endswith('ms'):
        dur = dur.rstrip('ms')
        ste += int(dur[-3:])
        dur = dur[:-3]
    if dur.endswith('s'):
        dur = dur.rstrip('s')
        ste += int(dur[-2:]) * 1000
        dur = dur[:-2]
    if dur.endswith('m'):
        dur = dur.rstrip('m')
        ste += int(dur) * 60000
    return int(ste / dt)

            
if __name__ == '__main__':
    eva = evaCure()
    eva.updateTillEnd()
