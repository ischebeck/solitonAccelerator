# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const

waveguideN = 2 #SiN

solLambda0 = 1.5 #um (in vacuum)
solLambdaN = 1.5/waveguideN #um (in waveguide)
solW0 = 2*np.pi*const.c/solLambda0*(1e6*1e-15) #1e15 rad/s
solEnergy = 10 #nJ
solT = 20. #fs
solArea = 1.5**2 #um2
fieldFactor = 0.005 #fraction of soliton field in electron channel

wigglerLambda = 1.5 #um 
wigglerW0 = 2*np.pi*const.c/wigglerLambda*(1e6*1e-15) #1e15 rad/s
wigglerPeak = 1000.

def peakField():
    #in GV/m
    #estimation based on gaussian pulse
    power = solEnergy/solT #MW
    intensity = power/solArea #1e18 W/m2
    field = np.sqrt(2*intensity / (const.c*const.epsilon_0*waveguideN**2)) #GV/m
    return field

solPeakField = peakField()


def solitonField(t, z):
    #in GV/m
    xi = 2*np.pi/solLambdaN*z - solW0*t
    shape = 2./(np.exp(xi/(solW0*solT))+np.exp(-xi/(solW0*solT)))
    amp = shape*np.exp(1j*xi)*solPeakField*fieldFactor
    #return 1.
    return amp

def solitonEnvelope(t, z):
    return np.abs(solitonField(t,z))

def wigglerField(t, z):
    #in GV/m
    xi = wigglerW0*t
    amp = np.exp(1j*(xi))*wigglerPeak
    #return np.array([1,0])
    return np.array([amp,0])

"""
z= np.linspace(-25, 25, 5000)
fig = plt.figure(figsize = (9,6), dpi = 200)
ax = fig.add_subplot(111)
ax.plot(z, np.real(solitonField(0, z)))
ax.plot(z, np.real(solitonEnvelope(0, z)))
ax.set_xlabel('z / um')
ax.set_ylabel('E-Field / GV/m')
ax.title.set_text('Pulse Shape (20 fs)')
plt.show()
"""


