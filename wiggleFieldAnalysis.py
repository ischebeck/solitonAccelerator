# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d, interp1d
import scipy.constants as const
from fieldIntegration import *
from scipy.fftpack import fft
#plt.close('all')

def loadData(path):
    x = np.loadtxt(path+'_x.txt')
    y = np.loadtxt(path+'_y.txt')
    real = np.loadtxt(path+'_REx.txt', dtype = np.complex)
    imag = np.loadtxt(path+'_IEx.txt', dtype = np.complex)
    e = real+1j*imag
    return x*1e6, y*1e6, e

def plotField(x, y, e):
    fig = plt.figure(figsize = (9,6), dpi = 200)
    ax = fig.add_subplot(111)
    ax.imshow(e.T, aspect = 'auto', origin = 'lower', extent = [x[0], x[-1], y[0], y[-1]])

def plotSlice(ax, x, y, e, yslice=0.):

    yidx = np.argmin(np.abs(yslice-y))
    ax.plot(x, e[:, yidx])

def fourier(x,val, ax1, ax2):
    dx = x[1]-x[0]    
    N = len(val)
    f = val - np.average(val)
    ff = fft(f)
    xf = np.linspace(0.0, 1.0/(2.0*dx), N//2)

    ax1.plot(x, val)
    ax2.plot(2*np.pi*xf, 2.0/N * np.abs(ff[0:N//2]), label = 'width = ' + str(widths[i]/1000.) + ' um')
    
path = 'C:\\Users\\hermann_b\\switchdrive\\Lumerical\\Geometry\\cross'
widths = [700, 750, 800, 850]
ydist = 0.05

fig = plt.figure(figsize = (9,6), dpi = 200)
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

for i in range(4):
    

    yslice = ydist + widths[i]/2000.
    
    x, y, Ex = loadData(path+str(widths[i]))
    yidx = np.argmin(np.abs(yslice-y))
    xmin = np.argmin(np.abs(3.-x))
    ExSlice = Ex[:, yidx]
        
    #fourier(x[xmin:], np.real(ExSlice[xmin:]), ax1, ax2)
    fourier(x[xmin:], np.abs(ExSlice[xmin:]), ax1, ax2)
    
ax1.set_xlabel(r'$x$ / um')
ax1.set_ylabel(r'$|E_x|$ / GV/m')
ax1.set_title(r'Field Profile along x')
ax2.set_xlabel(r'$k_x$ / um$^{-1}$')
ax2.set_ylabel(r'Intensity / a.u.')
ax2.set_title(r'FFT')
ax2.set_xlim((0,12))
ax2.legend()
fig.tight_layout()
