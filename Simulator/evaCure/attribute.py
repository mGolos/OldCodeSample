#!/usr/bin/env python
#-*- coding:Utf-8 -*-

''' @author: mathieu.golos@gmail.com
'''

class attributeClass:
    '''Base class for parameters attribution.
    Parameters:
    - _debug : return warnings for bad inputs if True
    
    Functions:
    - __init__
    - setAttr
    - getAttr
    - printAttr
    - verAttr.
    '''    
        
    def __init__(self, _debug=False, **kw):
        self._debug = _debug
        self.setAttr(**kw)
        self.verAttr()

    def setAttr(self, **kw):
        '''Change specific parameters from a dictionary.'''
        for p in kw.keys():
            if p in self.keys:
                setattr(self, p, kw[p])
            elif self._debug:
                print 'Bad attribute %s' %p         

    def getAttr(self, args=None):
        return getAttr(self, args)

    def printAttr(self, args=None):
        return printAttr(self, args)
    
    def verAttr(self):
        for k in self.keys:
            try: 
                self.getAttr([k])
            except:
                raise ValueError('%s not defined.'%k)
        

def printAttr(kwa, args=None):
    '''Print specific or all the parameters of a dictionnary or
    an attribute class.'''
    if type(kwa) == dict:
        if args is None:
            args = kwa.keys()    
        for p in args:
            print '%s: '%p, kwa[p]
    else:
        if args is None:
            args = kwa.keys
        for p in args:
            print '%s: '%p, getattr(kwa, p)
            
            
def getAttr(kwa, args=None):            
    '''Return specific or all the parameters in a dictionary.'''
    nruter = {}
    if type(kwa) == dict:
        if args is None:
            args = kwa.keys()
        for p in args:
            tmp = kwa[p]
            if tmp is not None:
                nruter[p] = tmp
    else:
        if args is None:
            args = kwa.keys
        for p in args:
            tmp = getattr(kwa, p)
            if tmp is not None:
                nruter[p] = tmp
    return nruter
