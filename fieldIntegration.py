# -*- coding: utf-8 -*-
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.integrate import ode
import os
import scipy.constants as const
    
m0 = const.electron_mass*const.c**2/const.electron_volt #eV

def genOdeFunction(field):
    # field(t, z) electric field at position z, at time t
    def func(t,y):
        #fs, [um, eV]
        z = y[0] # um
        p = y[1] # ev/c
        E = np.sqrt(m0**2 + p**2) #eV
        dz = p/E*const.c*1e-15*1e6 #um/fs
        dp = np.real(field(t, z)*(1e9*1e-15)*const.c)
        
        return [dz, dp]
    
    return func

def genOdeFunction2D(field):
    # field(t, z) electric field at position z, at time t
    def func(t,y):
        #fs, [um, eV]
        x = y[0] # um
        z = y[1] # um
        px = y[2] # ev/c
        pz = y[3] # ev/c
        E = np.sqrt(m0**2 + px**2 + pz**2) #eV
        dx = px/E*const.c*1e-15*1e6 #um/fs
        dz = pz/E*const.c*1e-15*1e6 #um/fs
        dpx, dpz = np.real(field(t, z)*(1e9*1e-15)*const.c)
        
        return [dx, dz, dpx, dpz]
    
    return func

def trackParticle(field, t0, y0, tEnd):
    # field(t, z) electric field at position z, at time t
    # t0, y0 units: fs, [um, eV]
    r = ode(genOdeFunction(field))
    r.set_initial_value(y0, t0)
    dt = 0.05 #fs
    trace = []
    while r.successful() and r.t < tEnd:
        r.integrate(r.t+dt)
        trace.append(np.hstack((r.t, r.y)))

    return np.array(trace)

def rel_gamma(beta):
    return 1./np.sqrt(1-beta**2)
def rel_p(beta):
    return m0*rel_gamma(beta)*beta
def rel_E(p):
    return np.sqrt(m0**2 + p**2)
def rel_pE(E):
    return np.sqrt(E**2-m0**2)
def rel_beta(p):
    return p/rel_E(p)
def rel_betaE(E):
    return rel_beta(rel_pE(E))