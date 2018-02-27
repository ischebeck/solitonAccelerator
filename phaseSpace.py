# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from fieldIntegration import *
from pulse import solitonField
import scipy.constants as const

field = solitonField
t0 = 0
tEnd = 2000

def plotTrace(ax, trace, frameSpeed = 0.):
    
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

fig = plt.figure(1, figsize = (9,6), dpi = 200)
ax = fig.add_subplot(111)

pInit = np.linspace(0.45, 0.55, 101)
l = [] #for acceleration length
for p in pInit:
    y0= [0,rel_p(p)]
    trace = trackParticle(field, t0, y0, tEnd)
    plotTrace(ax, trace)
    l.append(lAccel(trace))

fig = plt.figure(2, figsize = (9,6), dpi = 200)
ax = fig.add_subplot(111)
ax.plot(pInit, l, label = 'fieldFactor: 0.005')
ax.axhline(1.5/2., ls = '--')
ax.set_xlabel(r'$\beta$')
ax.set_ylabel(r'$L_{accel}$ / um')
plt.show()