# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from fieldIntegration import *
from pulse import solitonField, wigglerField
import scipy.constants as const

field = wigglerField
t0 = 0
tEnd = 200

def plotTrace(ax, trace, frameSpeed = 0.):
    
    ax.plot(trace[:,0], trace[:,1])
    #ax.plot(trace[:,1]-trace[:,0]*const.c*frameSpeed*1e6*1e-15, (rel_E(trace[:,2]) - m0)*1e-3)
    #ax.plot(trace[:,1], (rel_E(trace[:,2]) - E0_e)*1e-3)
    ax.title.set_text('e Propagation, frame speed: ' + str(frameSpeed) + ' c')
    
    
fig = plt.figure(1, figsize = (9,6), dpi = 200)
ax = fig.add_subplot(111)

betaInit = [0.7]

for b in betaInit:
    y0= np.array([0, 0, 0, rel_p(b)])
    trace = trackParticle(field, t0, y0, tEnd)
    plotTrace(ax, trace)
    
plt.show()