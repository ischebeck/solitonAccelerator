# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from fieldIntegration import *
from pulse import solitonField
import scipy.constants as const

field = solitonField
t0 = 0
tEnd = 2000

def plotTrace(ax, trace, frameSpeed = 0.5):
    
    #ax.plot(trace[:,0], rel_beta(trace[:,2]))
    ax.plot(trace[:,1]-trace[:,0]*const.c*frameSpeed*1e6*1e-15, (rel_E(trace[:,2]) - m0)*1e-3)
    #ax.plot(trace[:,1], (rel_E(trace[:,2]) - E0_e)*1e-3)
    ax.set_xlabel(r'$z-vt$ / um')
    ax.set_ylabel(r'$E_{kin}$ / keV')
    ax.title.set_text('e Propagation, frame speed: ' + str(frameSpeed) + ' c')

def lAccel(trace):
    #acceleration length in um
    maxPos = np.argmax(trace[:,2])
    minPos = np.argmin(trace[:,2])
    return abs(trace[maxPos,1]-trace[minPos,1])/2.

def maxEGain(trace):
    return (max(rel_E(trace[:, 2]))-min(rel_E(trace[:, 2])))/2.

fig = plt.figure(1, figsize = (9,6), dpi = 200)
ax = fig.add_subplot(111)

betaInit = np.linspace(0.48, 0.52, 31)
l = [] #for acceleration length
e = [] #for max e gain

for b in betaInit:
    y0= [0,rel_p(b)]
    trace = trackParticle(field, t0, y0, tEnd)
    plotTrace(ax, trace)
    l.append(lAccel(trace))
    e.append(maxEGain(trace))
    
fig = plt.figure(2, figsize = (9,6), dpi = 200)
ax = fig.add_subplot(111)
ax.plot(betaInit, l, label = 'fieldFactor: 0.005')
ax.set_xlabel(r'$\beta$')
ax.set_ylabel(r'$L_{accel}$ / um')
ax.legend()

fig = plt.figure(3, figsize = (9,6), dpi = 200)
ax = fig.add_subplot(111)
ax.plot(betaInit, np.array(e)*1e-3, label = 'fieldFactor: 0.005')
ax.set_xlabel(r'$\beta$')
ax.set_ylabel(r'$\Delta E $ / keV')
ax.legend()
plt.show()