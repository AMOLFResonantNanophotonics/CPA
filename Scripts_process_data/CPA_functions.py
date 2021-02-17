'''
Routine implementing the changepoint analysis of 

Watkins & Yang, J. Phys. Chem. B, vol 109, no 1, 2005

outlined and benchmarked in Palstra et al.

 

''' 


import numpy as np

def IdentifyOneJump(timestamps, skepticism):
    (B, log10SumB, scale) = BayesianlistTimeStamps(timestamps)

    if len(B)>5: # if there are less than 5 photons in the segment, we assume there are no changepoints, or at least they are not resolvable
        tjump = np.argmax(B[2:-2])+2 # cut off the edges of the curve, so you can't select the first or last two elements of the interval as changepoints
        # a more correct handling of this would be to use the method in fig 3, Watkins and Yang 2005.  
        tjump = tjump*(log10SumB>skepticism) # returns tjump or 0
        # note that tjump is a click nr. If you want the 'real' time, you need to take timestamps[tjump]
    else:
        tjump=0

    return tjump


def BayesianlistTimeStamps(timestamps):
    N = len(timestamps)
    Tmax = np.amax(timestamps)
    k = np.arange(2, N)
    Tk = timestamps[k-1]

    Vk = Tk/Tmax

    #B = (k*np.log(k/Vk) + (N-k)*np.log((N-k)/(1-Vk)) - N*np.log(N))/np.log(10)
    B1 = k*np.log(k/Vk)
    B2 = (N-k)*np.log((N-k)/(1-Vk))
    B3 = N*np.log(N)
    B4 = np.log(10)
    B = (B1+B2-B3)/B4

    scale = np.ceil(max(B))

    B_scaled = 10.**(B-scale)
    log10SumB = np.log10(np.sum(B_scaled))+scale

    return (B_scaled, log10SumB, scale)


def FindMultipleJumps(timestamps, skepticism, printsignal=True):
    if printsignal:
        print('\nFinding changepoints')
    N = len(timestamps)
    tjump = IdentifyOneJump(timestamps, skepticism)
    
    # find your first jump
    # If no jump is found at all, stop the routine
    if tjump==0:
        tslst=[]
        if printsignal:
            print("\nNO JUMPS WERE FOUND\n")
        return ([-1],[-1])
    
    intervals = 2
    startlst = [0,tjump-1] # these are indices, not times
    endlst = [tjump,N] # these are indices, not times
    trial = np.array([1,1])
    recursion = 0
    
    while sum(trial)>0:
        newtslst = []
        newtrial = []
        
        for q in range(intervals):
            insertnewt=[]
            inserttrial=[]
            if endlst[q]-startlst[q]<7:
                trial[q] = 0
            if trial[q]==1:
                tjump = IdentifyOneJump(timestamps[startlst[q]:endlst[q]]-timestamps[startlst[q]], skepticism)
                
                if tjump>0:
                    insertnewt = [tjump+startlst[q]]
                    inserttrial = [1]

            newtslst = newtslst + insertnewt + [endlst[q]]
            newtrial = newtrial + inserttrial + [trial[q]*tjump>0]
 
        recursion += 1
        tslst = newtslst[:-1] 
        trial = newtrial
        intervals = len(trial)
        
        startlst = [0] + tslst
        endlst = tslst + [N-1]

    retrievedjumpindices = np.array(tslst)
    retrievedjumptimes_bin = timestamps[retrievedjumpindices]

    return (retrievedjumpindices, retrievedjumptimes_bin)
