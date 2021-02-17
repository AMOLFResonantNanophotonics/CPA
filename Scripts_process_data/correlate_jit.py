'''
This script's function is to manage and test the 
Wahl algorithm for calculating correlations of 
time-tagged data sets

that means it is used for 
(1) g2  [correlate two detector channels]
(2) fluorescence decay traces  [correlate detector and laser channel]
(3) logarithmic time axis correlation [as you would use in FCS, or in blinking analysis Houel et al. ACS Nano 2015, 9, 1, 886–893]


Wahl, Gregor, Patting and Enderlein
Optics Express Vol. 11, Issue 26, pp. 3583-3591 (2003)

This implementation uses and needs jit accelleration,  and for good performance expects solely ints as inputs for time tags

'''

import numpy as np
import math
from numba import jit
# import time

 

# =============================================================================
# non-coarsened correlation of time traces 
# i.e., g(2) between any two channels, at resolution of timing card
# =============================================================================
@jit(nopython=True,parallel=False)   # use Numba just in time precompil. for speed up
def Corr(times1, times2, tau):  # Basic Wahl algorithm, with no time axis coarsening
    
        
    result=0*tau
    
    len1 = len(times1)
    len2 = len(times2)

    ## implementation of section 2.4 in Wahl paper     
    for lag in np.arange(len(tau)):
        
        swi = 1
        times2prime=times2+tau[lag]
        p1=0
        p2=0
  
        while (p1<len1 and p2<len2): 

                
            if swi == 1:
                if times1[p1] >= (times2prime[p2]):
                    if times1[p1] == (times2prime[p2]):
                        result[lag] = result[lag] + 1
                        p1 = p1+1
                    swi = -swi
                else:
                    p1 = p1 + 1
                    
            else:
                if (times2prime[p2]) >= times1[p1]:
                    if times1[p1] == (times2prime[p2]):
                        result[lag] = result[lag] + 1
                        p2 = p2+1
                    swi = -swi
                else:
                    p2 = p2 + 1
           
    return result


# =============================================================================
# logarithmic coarsening of time axis for long time g(2) correlations 
# [Wahl & Enderlein algorithm]
# =============================================================================
    
def NormalizedAutoCorrLog(times,coarseningbase):
    '''
    Wrapper that simply runs
    (1) time-axis and weight generator [2.5 in Wahl paper]
    (2) correlator, section 2.4 and 2.5 in Wahl paper
    (3) normalization, since Eq. (1) in Wahl is not normalized. 
        In itself the correlator returns the number k of coincidences. The average follows from k / N with N the total number of timebins
        Average intensity is simply the total number of photons k1 in the time bins N.
    '''
    lentime=len(times)
    maxtime=times[lentime-1]
    tau, scalingfactor = Taucommensurable(maxtime,coarseningbase)
    result = CorrLog(times,times,tau,scalingfactor)
    normalization=maxtime/(lentime*lentime)
    result=result*normalization
    return tau, result


def CorrLog(times1, times2, tau, scalingfactor):
    '''
        Wahl - logarithmic correlator
        Takes as input tau, shift times and the scaling factors
    ''' 
    weight1 = np.ones(np.shape(times1))
    weight2 = np.ones(np.shape(times2))
    result = np.zeros(np.shape(tau))
    
    
    tauset=np.arange(len(tau))
    print('Running Wahl logarithmic g(2) correlator')
    counter=0
    for lag in tauset:        
        if lag>0:
            tmp = int(np.round(np.log(scalingfactor[lag]/ scalingfactor[lag-1]) /np.log(2)))
            if tmp > 0:
                counter=counter+1
           
            
            for steps in np.arange(tmp):
                # merge the adjacent times,  account for the weights
                (times1,weight1)=shrinktimes(times1,weight1)
                (times2,weight2)=shrinktimes(times2,weight2)
                      
        
        # do a single step of the Wahl section 2.4, 2.5 correlator
        TauEff = math.floor(tau[lag]/scalingfactor[lag])
        result[lag]=SingleCorr(times1,times2,TauEff,weight1,weight2)
            
    result =result/scalingfactor
    return result


def shrinktimes(times,weight):
      '''
        Aux routine -  an input timetag list with weights is doubled in time step size
      '''
      times = np.right_shift(times, 1) # equivalent to dividing by 2, but faster
#     indeq = None
      indeq = np.diff(times)==0
      idx1=np.insert(indeq,len(indeq),False)
      idx2=np.insert(indeq,0,False)
      weight[idx1] = weight[idx1] + weight[idx2]
      
      idx1=np.logical_not(idx2)
      times = times[idx1]
      weight = weight[idx1]
      
      return times, weight


@jit(nopython=True)
def SingleCorr(times1, times2, tau, weight1, weight2):  # single step of Wahl for use with weights
## for use with a single tau.    
        
    result=0
    len1 = len(times1)
    len2 = len(times2)
 
 #   for lag in tauset:
    swi = 1
    times2prime=times2+tau
    p1=0
    p2=0
    
    while (p1<len1 and p2<len2): 
       if swi == 1:
            if times1[p1] >= (times2prime[p2]):
                if times1[p1] == (times2prime[p2]):
                    result = result + weight1[p1] * weight2[p2]
                    p1 = p1+1
                swi = -swi
            else:
                p1 = p1 + 1
                
       else:
            if (times2prime[p2]) >= times1[p1]:
                if times1[p1] == (times2prime[p2]): # shouldn't the indices 1 and 2 be reversed here?
                    result = result + weight1[p1] * weight2[p2]
                    p2 = p2+1
                swi = -swi
            else:
                p2 = p2 + 1
    return result



def Taucommensurable(recordlength, base):
    '''generator of logarithmic time steps and weights / Wahl section 2.5 '''
    tau = [0]
    scalingfactor = [1]
    count = 0
    b = 0
    
    while tau[-1] < recordlength:
        tau = tau + [tau[count]+scalingfactor[count]]
        scalingfactor = scalingfactor + [scalingfactor[count]]
        count = count + 1
        b = b + 1
        if (tau[count]%(scalingfactor[count]*2)==0  and  b>=base):
            b = 0
            scalingfactor[count] = scalingfactor[count]*2
    scalingfactor = [1] + scalingfactor[:-1]
    return (np.array(tau), np.array(scalingfactor))



