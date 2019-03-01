#!/usr/bin/env python
#-*- coding:Utf-8 -*-

from pylab import linspace, sqrt, exp, sin, norm, einsum, dot, tanh, arctanh
from scipy.signal import butter, lfilter


def fft(TS, tmax=None, sr=1, nf=1):
    '''Return the absolute of the real part of the fast fourrier transformation 
    from a time serie 'TS' with a sampling rate 'sr' and a fft resolution 'nf'.
    '''
    from numpy import linspace
    from numpy.fft import fft
    
    N = len(TS)
    if not tmax:
        tmax = N
    #frequencies = sr * linspace(0, tmax/2., N/2., endpoint=True) / nf
    frequencies = linspace(0, tmax/2., N/2., endpoint=True) / nf
    amplitudes = abs(fft(TS, n=nf*N)[:N/2] / N) 
    return frequencies, amplitudes    


def sigmoid(x, H=.5, Q=1, G=60, P=1, T=0.5):
    '''Sigmoidal thresholding function.'''
    return H * (Q + tanh(G * (P * x - T)))


def sigmoidM1(A, H=.5, Q=1, G=60, P=1, T=0.5):
    '''Inverse of sigmoidal thresholding function.'''
    #TODO Add epsilon error, Q and H
    return (arctanh(1.99 * (A - 0.5)) / G + T) / P


def HRF(ms = 1.):
    ''' The Heamodynamic response function
    '''
    T     = 10
    tau_s = 0.8
    tau_f = 0.4
    scale = 1000. / ms
    v_time= linspace(0., T, scale * T)
    sqrt_tfts = sqrt(1./ tau_f - 1./ (4.* tau_s ** 2))
    exp_ts    = exp(-0.5 * (v_time / tau_s))
    h         = exp_ts * sin(v_time * sqrt_tfts) / sqrt_tfts
    return h


def Lyapunov(v, P, W):
    '''W is a connectome
    '''
    W /= norm(W)
    theta = W.sum(0) * 0.5
    if v.ndim == 1:
        return -P * 0.5 * v.dot(W).dot(v) + (theta * v).sum()
    else:
        return -P * 0.5 * einsum('ij,ij->i', dot(v, W), v) + (theta * v).sum(1)


def partition(V, P, beta, **kwa):
    return exp(- beta * Lyapunov(V, P, **kwa)).sum()


def probaLya(V, P=1, beta=1, **kwa):
    return exp(- beta * Lyapunov(V, P, **kwa)) / partition(V, P, beta)


def butterBandPass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butterBandPassFilter(data, lowcut, highcut, fs, order=5):
    '''Time has to be the last dimension.
    fs = sample rate.
    '''
    b, a = butterBandPass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y
