#!/usr/bin/env python
#-*- coding:Utf-8 -*-

''' @author: mathieu.golos@gmail.com
'''

from pylab import randint, rand, randn, dot, tanh, norm, exp, sqrt
from pylab import where, arctanh, mean, permutation, zeros
from tools import data2array
from attribute import attributeClass, getAttr


def choice(**kwa):
    '''Choose between:
    - HopfieldBasedDynamic
    - HopfieldBasedStatic
    - ReducedWongWang
    - DemiWilsonCowanNet
    - HopfieldNetDiscreet.
    '''
    try:
        mod = getAttr(kwa)['model']
    except IOError:
        print 'Model undefinied.'
        raise
    
    if mod == 'HopfieldBasedDynamic':
        return HopfieldBasedDynamic(**kwa)
    
    elif mod == 'HopfieldBasedStatic':
        return HopfieldBasedStatic(**kwa)
    
    elif mod == 'ReducedWongWang':
        return ReducedWongWang(**kwa)
    
    elif mod == 'DemiWilsonCowanNet':
        return DemiWilsonCowanNet(**kwa)
    
    elif mod == 'HopfieldNetDiscreet':
        return HopfieldNetDiscreet(**kwa)
    
    elif mod == 'Generic2DOscillator':
        return Generic2DOscillator(**kwa)
    
    else:
        print 'Wrong model name.'
        return 0    


class modelBase(attributeClass):
    '''Network models:
    - HopfieldNetDiscreet: for binary neurons
    - DemiWilsonCowanNet: based on
    - HopfieldBasedStatic: static threshold
    - HopfieldBasedDynamic: dynamic threshold,
    
    Parameters:
    - threshold : type of threshold (global or local)
    - dt        : timestep
    - taux      : time constant for potentials
    - tauT      : time constant for global inhibition theta
    - P         : Excitatory factor on Inhibitory factor Psi/Omega
    - G         : Inhibitory factor dot thresholding function factor Omega.gamma
    - x         : Potential value for each nodes
    - theta     : inhibition value
    - A         : activity of neurons for each nodes
    
    Given:
    - N      : number of neurons
    - conn   : connectivity matrix
    - delays : delay matrix
    
    Build:
    - x     : Potential value for each nodes
    - theta : inhibition value
    - A     : activity of neurons for each nodes
    - E     : energy of the system
    
    Functions:
    - __init__
    - initialRand
    - initialExt
    - initialPatt
    - initialFile
    - initialNoise
    - energy_ThetaLocal
    - energy_ThetaGlobal
    - thresholding.
    '''
    keys = ['threshold', 'dt', 'taux', 'tauT', 'P', 'G', 
            'x', 'A', 'theta', 'E',
            'N', 'conn', 'delays']
            
    threshold = 'global'
    dt = 0.1
    taux = 10.
    tauT = 50.
    P = 1.
    G = 60.
    
    x = None
    A = None
    theta = 0.5
    E = None
   
   
    def __init__(self, **kwa):
        '''Construction from attributeClass class.'''
        attributeClass.__init__(self, **kwa)
       

    def density(self, x):
        '''Return the density of a vector as...
        '''
        return where(x < mean(x), 0., 1.).mean()
    

    def initialRand(self, dens, normX=None):
        '''Initialize (x, A) with a certain density 'dens'.
        '''
        #self.x = rand(self.N)
        #self.x[self.x <  (1.- dens)] = 0
        #self.x[self.x >= (1.- dens)] = 1
        
        #if normX != None:
            #self.x *= normX / norm(self.x)
            
        #self.A = self.thresholding(self.x, dens)
        
       
        self.A = zeros(self.N)
        self.A[permutation(self.N)[range(int(dens * self.N))]] = 1.
        self.initTheta()
        self.x = (arctanh(1.99 * (self.A - 0.5)) / self.G + self.theta) / self.P
        #self.x = where(self.A == 1., (self.theta + 2./ self.G) / self.P \
                                   #, (self.theta - 2./ self.G) / self.P)
        #self.x = self.A.copy()


    def initialExt(self):
        if self.x is None:
            if self.A is None:
                print 'No variables loaded in initialExt.'
            else:
                self.x = self.thresholdingM1(self.A)
        elif self.A is None:
            self.A = self.thresholding(self.x, self.density(self.x)) #self.theta)


    def initialPatt(self, ksix, ksiA):
        '''Initialize the network (x, A).'''
        self.x = ksix.copy()
        self.A = ksiA.copy()
           

    def initialFile(self, adress):
        '''Initialize the network (x, A) from the file pattern in adress.'''
        self.x = data2array(adress.replace('_S','_X'))
        self.A = data2array(adress.replace('_X','_S'))
        
        
    def initialNoise(self, noise):
        '''Add noise to initialization.'''
        self.x += noise
        self.A  = self.thresholding(self.x, self.theta)


    def energy_ThetaLocal(self):
        '''Returns the energy of the system for a local theta.'''
        self.E = - self.P * 0.5 * dot(self.conn.dot(self.A), self.A)\
                 + (self.theta * self.A).sum()
        return self.E


    def energy_ThetaGlobal(self):
        '''Returns the energy of the system for a global theta.'''
        self.E = - self.P * 0.5 * dot(self.conn.dot(self.A), self.A)\
                 + self.theta * self.A.sum()
        return self.E


    def thresholding(self, x, theta):
        '''Thresholding definition depending on threshold_type: 0 for tanh ; 1 for Brunel's function.'''
        return 0.5 * (1. + tanh(self.G * (self.P * x - theta)))


class HopfieldBasedDynamic(modelBase):
    '''Based on the Continuous Hopfield Network model.'''

    def __init__(self, **kw):
        '''Construction from modelBase class.'''
        modelBase.__init__(self, **kw)
        if self.taux: self.v1 = self.dt / self.taux
        else:         self.v1 = 0.
        if self.tauT: self.v2 = self.dt / self.tauT
        else :        self.v2 = 0.

        if self.threshold == 'global':
            self.initTheta = self.initThetaGlobal
            self.energy    = self.energy_ThetaGlobal
            self.meanoupas = lambda A: A.mean()
        elif self.threshold == 'local':
            self.initTheta = self.initThetaLocal
            self.energy    = self.energy_ThetaLocal
            self.meanoupas = lambda A: A
        else:
            print 'Threshold type bad definied : %s' %self.threshold
            
        self.theta_0 = 0.5

    def initThetaLocal(self):
        #self.theta = (self.conn.sum(0) * self.theta_0).A1
        self.theta = self.P * self.conn.dot(self.A)
        #self.theta = self.A.copy()

    def initThetaGlobal(self):
        #self.theta = self.A.mean()
        #self.theta = 1./ self.A.mean() * self.P * 0.0042
        self.theta = self.P * self.conn.dot(self.A).mean() / self.A.mean()**0.5

    def update(self, noises):
        '''Dynamical update of activity procedure with the thresholding function
        and a probability for some neurone to change his state.
        '''
        self.x     += self.v1 * (- self.x     + self.conn.dot(self.A)  + noises[0])
        self.theta += self.v2 * (- self.theta + self.meanoupas(self.A) + noises[1])
        self.A      = self.thresholding(self.x, self.theta)
            
        
        '''DELAYS'''
        #TODO
        #IAN = take(self.TA, self.delays).reshape(self.N, self.N)
        ## For 'conn x IAN' the sum is done on the second axis to match the asymetry
        #self.x     += self.v1 * (- self.x     + self.conn.multiply(IAN).sum(1).A1 + s2mx)
        #self.theta += self.v2 * (- self.theta + self.sumoupas(self.TA[0])      + s2mT))
        #self.TA     = roll(self.TA, 1, axis=0)
        #self.TA[0]  = self.thresholding(self.x)


class HopfieldBasedStatic(modelBase):
    '''Continuous Hopfield Network.'''

    def __init__(self, **kw):
        '''Construction from modelBase class.'''
        modelBase.__init__(self, **kw)
        if self.taux: self.v1 = self.dt / self.taux
        else:         self.v1 = 0.
        if self.tauT: self.v2 = self.dt / self.tauT
        else :        self.v2 = 0.

        if self.threshold == 'local':
            self.initTheta = self.initThetaLocal
            self.energy = self.energy_ThetaLocal
        elif self.threshold == 'global':
            self.initTheta = self.initThetaGlobal
            self.energy = self.energy_ThetaGlobal
            
        self.magic_number = 0.5

    def initThetaLocal(self):
        self.theta = (self.conn.sum(0) * self.magic_number).A1
        self.theta_0 = (self.conn.sum(0) * self.magic_number).A1

    def initThetaGlobal(self):
        self.theta = ((self.conn.sum(0) * self.magic_number).A1).mean()
        self.theta_0 = ((self.conn.sum(0) * self.magic_number).A1).mean()

    def update(self, noises):
        '''Dynamical update of activity procedure with the thresholding function
        and a probability for some neurone to change his state.
        '''
        self.x += self.v1 * (- self.x + self.conn.dot(self.A) + noises[0])
        self.theta +=  self.v2 * (-(self.theta-self.theta_0)  + noises[1])
        #self.x += self.v1 * (- self.x + self.conn.dot(self.A) + self.theta_0 * noises[0] / self.theta_0.mean())
        #self.theta +=  self.v2 * (-(self.theta-self.theta_0)  + self.theta_0 * noises[1] / self.theta_0.mean())
        self.A  = self.thresholding(self.x, self.theta)


class ReducedWongWang(modelBase):
    '''Reduced Wong Wang.'''

    a     = 0.270
    b     = 0.108
    c     = 154
    gamma = 0.641
    tau_s = 100.
    w     = 0.6
    J_N   = 0.2609
    I_o   = 0.33
    theta = 0.5
    we    = 2.4
    err   = 1e-6
    S     = None
       
       
    def __init__(self, **kw):
        '''Construction from modelBase class.'''
        
        self.keys = ['dt', 'a', 'b', 'c', 'gamma', 'tau_s', 'w', 'tau_s', 'J_N', 'I_o', 
                'P', 'G', 'theta', 'we',
                'S', 'N', 'conn', 'delays']
        
        modelBase.__init__(self, keys=self.keys, **kw)
                    
        self.mn1  = self.w * self.J_N
        self.mn2  = self.dt / self.tau_s
        self.mn3  = self.J_N * self.we
        self.sqdt = sqrt(self.dt)
   
   
    def initialRand(self, dens, normX=None):
        '''...'''
        self.S = zeros(self.N)
        self.S[permutation(self.N)[range(int(dens * self.N))]] = 1.
        if normX != None:
            self.S *= normX / norm(self.S)


    def initialExt(self):
        pass
       
       
    def initialNoise(self, noise):
        '''Add noise to initialization.'''
        self.S += noise
        
        
    def initTheta(self):
        pass


    def thresholding(self, x):
        '''...'''
        #return 0.5 * (1. + tanh(self.G * (self.P * x - self.theta)))
        
        r = self.a * x - self.b
        
        ind = abs(self.c * r) < self.err  # Problematic indices (x too close to 0.4)
        if any(ind):
            nruter = zeros(len(x))
            nruter[ind] = 0.5 * r[ind] + 1./ self.c   # For numerical stability
            nruter[~ind] = r[~ind] / (1.- exp(- self.c * r[~ind]))
        else:
            nruter = r / (1. - exp(-self.c * r))
        
        return nruter
   

    def update(self, noises):    
        '''...'''        
        x = self.mn1 * self.S + self.I_o + self.mn3 * self.conn.dot(self.S)
        r = self.thresholding(x)
        self.S += self.mn2 * (- self.S / self.tau_s + (1.- self.S) * self.gamma * r) \
                + self.sqdt * noises[0]
        self.S[self.S < 0] = 0.
        self.S[self.S > 1] = 1.
        

class Generic2DOscillator(modelBase):
    '''Generic 2D Oscillator.'''

    tau = 1.
    I = 0.
    a = 0.
    b = -10.
    c = 0.
    d = 0.025
    e = 3.
    f = 1.
    g = 0.
    alpha = 1.
    beta = 1.
    gamma = 1.
    x = None
    A = None
       
       
    def __init__(self, **kw):
        '''Construction from modelBase class.'''
        
        self.keys = ['tau','I','a','b','c','d','e','f','g','alpha','beta','gamma',
                     'x','A',
                     'P', 'G', 'theta',
                     'N', 'conn', 'delays']
        
        modelBase.__init__(self, keys=self.keys, **kw)
                    
        self.v1  = self.dt / self.tau
   
   
    def initialRand(self, dens, normX=None):
        '''...'''
        self.x = zeros(self.N)
        self.A = zeros(self.N)
        self.x[permutation(self.N)[range(int(dens * self.N))]] = 1.
        self.A[permutation(self.N)[range(int(dens * self.N))]] = 1.
        if normX != None:
            self.x *= normX / norm(self.x)


    def initialExt(self):
        pass
        
        
    def initTheta(self):
        pass

       
       
    def initialNoise(self, noise):
        '''Add noise to initialization.'''
        self.x += noise
        
        
    def update(self, noises):    
        '''...'''       
        self.x += self.v1 * (self.d * (- self.x**3 + self.e * self.x**2 + self.A + self.conn.dot(self.x)) + noises[0])
        self.A += self.v1 * self.d * (self.a + self.b * self.x - self.A)
        
        
        
class DemiWilsonCowanNet(modelBase):
    '''Demi WilsonCowan Network.'''

    def __init__(self, **kw):
        '''Construction from modelBase class.'''
        modelBase.__init__(self, **kw)

    def initTheta(self):
        '''Initialize inputs from nodes activity using the inverse of thresholding definition.'''
        self.x0 = arctanh(1.999999999 * (self.A - 0.5 + 1e-15)) / self.gamma + self.theta
        self.x  = self.x0.copy()

    def thresholding(self, x):
        '''Thresholding definition depending on threshold_type: 0 for tanh ; 1 for Brunel's function.'''
        return 0.5 * (1. + tanh(self.gamma * (x - self.theta)))

    def update(self, noises):
        '''Dynamical update of activity procedure with the thresholding function
        and a probability for some neurone to change his state.
        '''
        if self.tauT:
          self.theta += self.dt / self.tauT * (- self.theta + self.Omega / self.N * self.A.sum())

        self.x = self.Psi * self.conn.dot(self.A)
        self.A += self.dt / self.taux * ( - self.A + noises[0] + self.thresholding(self.x) )


class HopfieldNetDiscreet(modelBase):
    '''Discreet Hopfield Network.'''

    def __init__(self, **kw):
        '''Construction from modelBase class.'''
        modelBase.__init__(self, **kw)

    def initTheta(self):
        '''Initialize inputs from nodes activity using the inverse of thresholding definition.'''
        self.x0 = where(self.A < 0.5, self.theta - 0.5, self.theta + 0.5)
        self.x = self.x0.copy()

    def thresholding(self, x):
        '''Thresholding definition.'''
        return where(x < self.theta, 0., where(x > self.theta, 1., x))

    def update(self, noises):
        '''Dynamical update of activity as Monte Carlo procedure with a probability for
        some neuron to change his state.
        '''
        if self.tauT:
            if self.test == 0:
                self.theta = self.A.sum() / self.N
            elif self.test == 1:
                self.theta = self.Omega * self.A.sum() / self.N
            elif self.test == 2:
                self.theta += self.dt / self.tauT * self.N * (- self.theta + self.Omega * self.A.sum() / self.N)
            elif self.test == 3:
                self.theta += self.dt / self.tauT * (- self.theta + self.Omega * self.A.sum() / self.N)

        k1 = randint(self.N)
        self.x[k1] = self.Psi * self.conn.getcol(k1).T.dot(self.A)[0]
        self.A[k1] = self.thresholding(self.x[k1])

        for nothing in xrange(self.N):
            k2 = randint(self.N)
            if rand() < self.sigma2x / self.N:
                self.A[k2] = self.A[k2] == 0
