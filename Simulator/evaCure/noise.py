#!/usr/bin/env python
#-*- coding:Utf-8 -*-

''' @author: mathieu.golos@gmail.com
'''

from tools import noiseClass
from attribute import attributeClass


class evaNoi(attributeClass):
    '''Return different number of white &/or pink noise:
    Parameters:
    - colors  : types of noise (None, white', 'pink', 'brown')
    - stdD_x  : standard deviation for potentials
    - stdD_T  : standard deviation for inhibition
    - AonF    : parameter for brownian noise (1: pink noise)
    
    Given:
    - threshold : type of threshold (global or local)
    - N : number of neurons
    
    Build:
    - noises : vector of different noises
    
    Functions:
    - __init__
    - configure
    - update
    - updateNoise
    - noNoise.
    '''
    keys = ['colors', 'stdD_x', 'stdD_T', 'AonF', 
            'N', 'threshold']

    colors = 'white'     
    stdD_x = 0.
    stdD_T = 0.
    AonF = 1
        
    
    def __init__(self, **kwa):
        '''Construction from attributeClass class.'''
        attributeClass.__init__(self, **kwa)
        self.noises = [[], []]
        self.configure()
        
        
    def configure(self):
        if type(self.colors) == str:
            self.colors = [self.colors] * 2
            
        if self.colors is None:
            self.colors = [self.colors] * 2
                   
        if self.threshold == 'global':
            dimT= 1
        elif self.threshold == 'local':
            dimT = self.N
        else:
            print "Wrong threshold: %s" %self.threshold
                
        self.noises[0] = noiseClass(self.colors[0], dim1=self.N, stdD=self.stdD_x)
        self.noises[1] = noiseClass(self.colors[1], dim1=dimT,   stdD=self.stdD_T)
        
        
    def update(self):
        return self.noises[0].update(), self.noises[1].update()
    
    
    def updateNoise(self, stdD_x=None, stdD_T=None):
        if stdD_x != None:
            self.noises[0].stdD = stdD_x
        if stdD_T != None:
            self.noises[1].stdD = stdD_T
