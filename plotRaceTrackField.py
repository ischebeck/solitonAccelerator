# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d, interp1d
import scipy.constants as const
from fieldIntegration import *
from scipy.fftpack import fft
from matplotlib import animation, rc



def loadData(path):
    x = np.loadtxt(path+'x.txt')
    y = np.loadtxt(path+'y.txt')
    real = np.loadtxt(path+'REx.txt', dtype = np.complex)
    imag = np.loadtxt(path+'IEx.txt', dtype = np.complex)
    e = real+1j*imag
    return x*1e6, y*1e6, e

fig = plt.figure(figsize = (9,6), dpi = 200)
ax = fig.add_subplot(111)

def plotField(x, y, e):

    im = plt.imshow(e.T, cmap = 'bwr', aspect = 'equal', origin = 'lower', extent = [x[0], x[-1], y[0], y[-1]], animated = True)
    ax.set_xlabel(r'$x$ (um)')
    ax.set_ylabel(r'$y$ (um)')
    ax.set_title(r'$E_x$')
    return im,
    
path = 'C:\\Users\\hermann_b\\switchdrive\\Lumerical\\Ring Resonator\\'

x, y, e = loadData(path)
r = 9
f = np.roll(np.hstack((e, np.zeros(np.shape(e)))), r, axis=1)[::-1]
f2 = f[:,::-1]

nFrames = 30
ph = np.linspace(0,np.pi*2, nFrames, endpoint = False)



def animate(i):
    return plotField(x,2*y,np.real((f+f2)*np.exp(-1j*ph[i])))


anim = animation.FuncAnimation(fig, animate, frames=nFrames, interval=20, blit=True, repeat= True)
anim.save('field.mp4', writer="ffmpeg", fps=30)



