'''
Auxiliary functions
'''

import numpy as np
import matplotlib as mpl

def gauss(x, mu, sigma, A):
    return(A*np.exp(-((x-mu)/sigma)**2))

# function for linear fits
def linear(x, a, b):
    return a - b*x


def PowerLaw_exp(params, t):
    [A, B, C] = params
    P = (A*(t)**-C) * np.exp(-B*t)
    return P

# least squares fit
def lstsq_PowerLaw_exp(fit_params, longcorr_xdata, longcorr_ydata, func_name, fit_guess):
    return np.sum((PowerLaw_exp(fit_params, longcorr_xdata)-longcorr_ydata)**2)



def Isacolmap():
    x = np.linspace(0,1, num=256)
#    xnew = np.linspace(0,1, num=256)
#    xnew = np.linspace(0,1, num=256)**1.5
    xnew = np.linspace(0,1, num=256)*(np.sin(np.linspace(0,0.5*np.pi, num=256)))
    y = mpl.cm.magma(np.arange(256))
#    y = mpl.cm.turbo(np.arange(256))
    np.shape(x)
    np.shape(y)
    yinterp = np.ones(np.shape(y))
    for i in [0,1,2]:
        yinterp[:,i] = np.interp(x, xnew, y[:,i])
    return yinterp




def IsacolmapFDID():
    x = np.linspace(0,1, num=256)
#    xnew = np.linspace(0,1, num=256)
#    xnew = np.linspace(0,1, num=256)**1.5
    xnew = np.linspace(0,1, num=256)*(np.sin(np.linspace(0,0.5*np.pi, num=256)))
#    y = mpl.cm.magma(np.arange(256))
    y = mpl.cm.turbo(np.arange(256))
    np.shape(x)
    np.shape(y)
    yinterp = np.ones(np.shape(y))
    for i in [0,1,2]:
        yinterp[:,i] = np.interp(x, xnew, y[:,i])
    return yinterp