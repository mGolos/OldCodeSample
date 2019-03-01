#!/usr/bin/env python
#-*- coding:Utf-8 -*-

''' @author: mathieu.golos@gmail.com
'''

from tools import data2array
from pylab import norm, dot
from scipy import sparse
from attribute import attributeClass


class evaCon(attributeClass):
    '''Structural connectivity.
    Parameters:
    - connAd   : adress of the structural connectivity matrix
    - distAd   : adress of the nodes distances matrix
    - vel      : defines the velocity for delays (0 disable delays)
    - tauC     : time constant for plasticity
    - normC    : defines the renormalization function (<0 disable it)
    - normT    : defines the normalization type (Frobenius by default)
    
    Given:
    - dt : time step
    
    Build:
    - conn   : structural density matrix
    - N      : number of nodes
    - delays : delays between nodes (if vel != 0)
    - conn0  : initial structural connectivity matrix for plasticity (if tauC != 0)
    
    Functions:
    - __init__
    - connectivity
    - initPlasticity
    - sparsity
    - delays
    - updatePlasticity
    - updateNo.
    '''
    keys = ['connAd', 'distAd', 'vel', 'tauC', 'normC', 'normT',
            'dt',
            'conn', 'N', 'delays', 'conn0']
    
    connAd = 'Connectomes/SC_H_998_0.npy'
    distAd = 'Connectomes/FL_H_998_0.npy'
    vel = 0
    tauC = 0
    normC = 1.
    normT = 'fro'
    
    conn = None
    N = None
    delays= None
    conn0 = None
    
    
    def __init__(self, **kw):
        '''Construction from attributeClass class.'''
        attributeClass.__init__(self, **kw)
        
        # Load the matrix
        self.connectivity(self.conn)
        
        # Generate delays
        if self.vel:
            self.delays()
        
        #Initial connectivity copy for plasticity
        if self.tauC:
            self.initPlasticity()
            self.v3 = self.dt / self.tauC
            self.update = self.updatePlasticity
        else:
            self.update = self.updateNo


    def connectivity(self, conn=None):
        '''Initialize the connectivity matrix from a data file,
        renormalizes and transform it into a sparse matrix.
        If plasticity is enable, generate also a copy.
        '''
        #Loading
        if conn is None:
            self.conn = data2array(self.connAd)
        else:
            self.conn = conn
        self.N = len(self.conn)
        
        #Renormalization
        if self.normC:
            self.conn *= self.normC / norm(self.conn, ord=self.normT)
            
        #Transformation to sparse matrix
        self.conn = sparse.csr_matrix(self.conn)


    updateConnectivity = connectivity
            
            
    def initPlasticity(self):
        self.conn0 = self.conn.copy()
        #self.connB = 1 * (self.conn != 0)


    def updatePlasticity(self, A):
        '''If there is plasticity, the connectivity is update here.'''
        self.conn = self.conn + self.v3 * (- self.conn + self.conn0 - self.conn0 * sparse.diags(A,0))
        
        
    def updateNo(self, *la, **kwa):
        pass
        

    def sparsity(self):
        '''Return the sparsity of the structural connectivity matrix.'''
        return self.conn.nnz / self.conn.shape[0] * self.conn.shape[1]
        

    def delays(self):
        #TODO verify and generalize
        ''' Create the indices matrix for the delays depending on
        nodes distances, timestep and velocity
        '''
        
        indices = np.ma.masked_array(mask)
        temporalvalues[indices]
        
        try:
            Dist = data2array(self.distAd)
        except ImportError:
            print 'Unable to load distance matrix from %s' %self.distAd
            raise

        if self._debug:
            print 'Maximal distance between nodes:', Dist.max(),'mm'
            
        if self.vel: self.delays = (Dist / self.vel / self.dt / 10.).astype(int)
        else:        self.delays = zeros(Dist.shape)
        vec = zeros(self.N * self.N)
        for i in xrange(self.N):
            for j in xrange(self.N):
                vec[i*self.N+j] = self.N * self.delays[i,j] + j
        self.delays = vec.astype(int)

        #Ij   = dot(ones((self.N,1)), [range(self.N)]).astype(int)
        #self.delays = [Ii, Ij]
        #self.delays = sparse.csr_matrix(self.delays)
