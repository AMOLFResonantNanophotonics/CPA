"""
Auxiliary functions for constructing FDID diagrams

Contains
makeGaussian() -  Gaussian entry into FDID from a given segment
MakeFDID() -  runs over an entire segment list to construct FDID
"""

import numpy as np
import math


def makeGaussian(xmin, xmax, xres, ymin, ymax, yres, xval, erx, yval, ery):
    """ Make a gaussian kernel.
    """

    x = np.arange(xmin, xmax, xres)
    y = np.arange(ymin, ymax, yres)

    x, y = np.meshgrid(x, y)
    
    # for very narrow Gaussians there is sometimes a pixelation and normalization issue
    # For very long segments relative errors in intensities can be small, catch it as follows: 

    if ery < 0.2*yres:
        ery=0.2*yres
    
    if erx < 0.2*xres:
        erx=0.2*xres
        
    func = np.exp(-np.square(x-xval)/(2.0*erx*erx)-np.square(y-yval)/(2.0*ery*ery))/(2*math.pi*erx*ery)
    
    return func


# def LoadFDID(loadFilepathFit):
#     Nheader = 0
#     with open(loadFilepathFit,"r") as f:
#         for nr, line in enumerate(f):
#             if line.startswith("# chunk nr, gamma fit"):
#                 Nheader = nr
    
#     gammaIntensitylst = np.loadtxt(loadFilepathFit, delimiter=',', skiprows=Nheader+1)
#     return gammaIntensitylst


def MakeFDID(gamma_lst, gamma_err_lst, I_lst, I_err_lst, countPriority=False, CPAlengthpriority=False, difftimelst=[], xmax=1):
    '''
    gamma_lst :      list of the decay rates for each segment/binned element
    gamma_err_lst :  list of the associated errors of these decay rates
    I_lst :          list of the intensities for each segment/binned element
    I_err_lst :      list of the associated errors of these intensities
    '''
    
    # =============================================================================
    # import preamble
    # =============================================================================
        
    
    if CPAlengthpriority and difftimelst==[]:
        print('You want to normalize the probability density by CPA length.\nYou need to give a value for difftimelst')
        return []
    
    xmin = 0 # in 1/ns
    ymin = 0 # in nr of counts
    ymax = 3*np.mean(I_lst)+2*np.mean(I_err_lst) # in nr of counts
    
    xres = (xmax-xmin)/70 # in 1/ns
    yres = (ymax-ymin)/70 # in nr of counts
    
    fdid = np.zeros((math.ceil((ymax-ymin)/yres), math.ceil((xmax-xmin)/xres)))
    
    
    for i in range(len(gamma_lst)):
    
        
        newGauss = makeGaussian(xmin, xmax, xres, ymin, ymax, yres, gamma_lst[i], gamma_err_lst[i], I_lst[i], I_err_lst[i])

        if countPriority:
            newGauss = newGauss*I_lst[i]
        if CPAlengthpriority:
            newGauss = newGauss*difftimelst[i]
        fdid += newGauss
    
    #fdid = fdid / len(gamma_lst)       #to normalize the FDID plot
    fdid = fdid / np.sum(np.sum(fdid))  #to normalize the FDID plot
  
    return fdid, (xmin, xmax,ymin, ymax)




