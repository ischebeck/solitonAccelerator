# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d, interp1d
import scipy.constants as const
from fieldIntegration import *
#plt.close('all')

def loadData(path):
#    data = open(path).read().split('\n\n')
#    x = np.array(data[1].split('\n')[1:], dtype = float)
#    y = np.array(data[2].split('\n')[1:], dtype = float)
#    e = np.array([np.array(l.split(' ')[1:], dtype = float) for l in data[3].split('\n')[1:]])
    x = np.loadtxt(path+'_x.txt')
    y = np.loadtxt(path+'_y.txt')
    real = np.loadtxt(path+'_REx.txt', dtype = np.complex)
    imag = np.loadtxt(path+'_IEx.txt', dtype = np.complex)
    e = real+1j*imag
    return x*1e6, y*1e6, e

def plotField(x, y, e):
    fig = plt.figure(figsize = (9,6), dpi = 200)
    ax = fig.add_subplot(111)
    ax.imshow(e.T, aspect = 'auto', extent = [x[0], x[-1], y[0], y[-1]])

def plotSlice(ax, x, y, e, yslice=0.):

    yidx = np.argmin(np.abs(yslice-y))
    ax.plot(x, e[:, yidx])

def plotTrace(ax, trace, beta0 = 1., frameSpeed = 0.):
    ax.plot(trace[:,1]-trace[:,0]*const.c*frameSpeed*1e6*1e-15, (rel_E(trace[:,2]) - rel_E(trace[0,2]))*1e-3)
    ax.set_xlabel(r'$z-vt$ / um')
    ax.set_ylabel(r'$\Delta E_{kin}$ / keV')
    ax.title.set_text(r'$\beta$ = ' + str(beta0) + ' c, frame speed: ' + str(frameSpeed) + ' c')

path = 'C:\\Users\\hermann_b\\switchdrive\\Lumerical\\Geometry\\cross'
widths = [700, 750, 800, 850]
ydist = 0.05

fig = plt.figure(1, figsize = (9,6), dpi = 200)
ax = fig.add_subplot(111)

for i in range(len(widths)):
    x, y, Ex = loadData(path+str(widths[i]))
    
    yslice = widths[i]/2000. + ydist

    yidx = np.argmin(np.abs(yslice-y))
    f = interp1d(x, Ex[:, yidx])
    
    w0 = const.c/1.55e-6*1e-15 # in 1e15 Hz
    
    def field(t, s):
        if s<x[0] or s>x[-1]:
            return 0.
        else:
            return np.real(f(s)*np.exp(-2.*const.pi*1j*w0*t))
    
    '''
    fig = plt.figure(figsize = (9,6), dpi = 200)
    ax = fig.add_subplot(111)    
    
    g= np.linspace(x[0], x[-1], 1000)
    for t in np.linspace(0, 2, 1):
        ax.plot(g, field(t,g))
    '''
    bInit = np.linspace(0.6, 0.99, 30)
    maxGain = np.zeros(len(bInit))
    for j, b in enumerate(bInit):
        print(j)
        #fig = plt.figure(figsize = (9,6), dpi = 200)
        #ax = fig.add_subplot(111)
        gain = []
        for phase in np.linspace(0, 1, 30):
            y0= [0,rel_p(b)]
            trace = trackParticle(field, phase/w0, y0, 100+phase/w0)
            #plotTrace(ax, trace, beta0 = b)
            gain.append(np.max((rel_E(trace[:,2]) - rel_E(trace[0,2]))*1e-3))
        maxGain[j] = np.max(gain)
        #ax.set_xlim(0,16)
    
    ax.plot(bInit, maxGain, marker='o', ls = '-', label = 'width = ' + str(widths[i]/1000.) + ' um')

ax.set_xlabel(r'$\beta_0$')
ax.set_ylabel(r'$\Delta E_{max}$ / keV')

ax.legend()