'''
This file contains the functions for the groupig/clustering algorithm from Watkins & Yang, J. Phys. Chem. B, vol 109, no 1, 2005
that do the actual work. THis includes 
(1) initial hierarchical clustering - InitialClustering
(2) with that clustering as input,  Bayesian Inf. Crit. based clustering into m-levels

'''

# import approxpoissondpf
import numpy as np
import numpy.matlib 
import matplotlib.pyplot as plt
from scipy.stats import poisson
import pandas as pd



def Grouping(segmentlengths, segmentcounts, mingroups, maxgroups):
    '''
    input
    segmentlengths  : the lengths of the CPA segments
    segmentcounts   : the nr of counts in each CPA segment
    mingroups       : the minimum nr of groups to test for
    maxgroups       : the maximum nr of groups to test for

    output
    groupingdata :   the nr of states and their likelihoods with 2 different methods from Watkins & Young
    mlikely_tot  :   most likely trajectory though the states
    psummarry    :   the likelihood of drawing a certain state, given m states. This is in TIME occupancy, not NUMBER occupancy (see definition of pm)
    '''

    initialcluster = InitialClustering(segmentlengths, segmentcounts, plot=False)
     
    Schw = np.zeros([maxgroups, 2])
    psummarry = np.zeros([maxgroups, maxgroups])
    mlst = np.arange(mingroups,maxgroups+1) # nr of levels under test

    mlikely_tot = np.zeros((maxgroups, len(segmentlengths)), dtype=int)
    for mgroups in mlst:
    
        (pm, Im, Tm, mlikely, schwartz, pmj) = ExpectationMaximizationClustering(segmentlengths, segmentcounts, mgroups, initialcluster[mgroups-1])
        Schw[mgroups-mingroups] = schwartz
        mlikely_tot[mgroups-1]  = mlikely
        psummarry[mgroups-1, 0:mgroups] = pm

    groupingdata = {'mlst':mlst, 'Schw0':Schw[:,0], 'Schw1':Schw[:,1]}
    groupingdata = pd.DataFrame(groupingdata)
    cols =['mlst', 'Schw0', 'Schw1']
    groupingdata = groupingdata[cols]
    return groupingdata, mlikely_tot, psummarry


def ApproxPoissPDF(k, mu):
    logpdf=k*np.log(mu)-mu-(k*np.log(k)-k+0.5*np.log(2*np.pi*k))
    
    idx = np.where(k==0)
    for i in idx:
        logpdf[i] = np.log(poisson.pmf(0,mu[i]))
    
    pdf = np.exp(logpdf)
    return (pdf, logpdf)


def ClusterOnce(TG, NG, assignment):
    '''
    This function finds the two most similar segments in an array and clusters them into one state
    TG: jumplengths
    NG: rNlevels
    '''
    
    # note that my meshgrids start at zero. Need to add 1 to get other indexing
    [Tm, Tj] = np.meshgrid(TG, TG) 
    [Nm, Nj] = np.meshgrid(NG, NG)
    [m, j]   = np.meshgrid(np.arange(len(NG)), np.arange(len(NG)))
    m = m+1
    j = j+1
    
    #should this be Mmj or Mjm?
    Mmj = (Nm+Nj)*np.log((Nm+Nj)/(Tm+Tj)) - (Nm)*np.log((Nm/Tm))-(Nj)*np.log((Nj/Tj));  # Eq 11
    
    idx = np.where(np.ndarray.flatten(np.triu(m,1),'F')>0)[0] # not the same as Femius' idx, but Python indexes start at 0
    
    winner_idx = np.argmax(np.ndarray.flatten(Mmj)[idx]) # take all nonzero elements of Mmj above the diagonal (which is a vector), and find the index of the largest one
    hitpair = np.unravel_index(idx[winner_idx], np.shape(Mmj)) # put this index back into 2d shape
    hitpair = sorted(hitpair)
    
    # combine most similar entries by adding the lengths and the counts
    TG[hitpair[0]] = TG[hitpair[0]]+TG[hitpair[1]]  # Eq 12
    NG[hitpair[0]] = NG[hitpair[0]]+NG[hitpair[1]]  # Eq 13
    
    # delete the entry that was squeezed into mlist
    TG = np.delete(TG, hitpair[1])
    NG = np.delete(NG, hitpair[1])
    
    idx = np.where(assignment==hitpair[1])
    assignment[idx] = hitpair[0]
    
    idx = np.where(assignment>hitpair[1])
    assignment[idx] = assignment[idx]-1
    
    return (hitpair, TG, NG, assignment)


def InitialClustering(TG, NG, plot=False):
    '''
    This function takes a list of segment times and counts, and clusters them together sequentially until there is only one state
    The output is used as an initial guess for the clustering defiend by Watkins & Yang (J. Phys. Chem. B 2005, 109, 617-628)
    This function is a lot faster than InitialClustering()
    
    TG    : lengths of the jump segments
    NG    : number of counts in the jump segments
    '''
    
    sublen = 300
    intermediatestop = 30
    
    tg2 = [] 
    ng2 = []
    hierarchyassignment2 =[]
    
    N_sublen = int(np.ceil(len(TG)/sublen))
    
    hierarchyassignment = np.zeros([len(NG)+1, len(NG)], dtype=int) # columns are nr of iterations, rows are segment nr
    hierarchyassignment[0] = np.arange(len(NG))
     
    if len(TG) > sublen:
        # =============================================================================
        # Divide the time series into sub-lengths and partially cluster
        # =============================================================================
        
        for i in range(N_sublen):
            # divide into sub-lengths
            tg2.append(TG[i*sublen:i*sublen+sublen])
            ng2.append(NG[i*sublen:i*sublen+sublen])
            hierarchyassignment2.append(hierarchyassignment[:,i*sublen:(i+1)*sublen])
        
            # set the assignments to start staggered
            # not strictly necessary, but makes for pretty pictures
            hierarchyassignment2[i] = hierarchyassignment2[i]%sublen + i*intermediatestop
            
            # do the clustering until you have intermediatestop levels
            p = 0
            while len(ng2[i]) > intermediatestop:  
                (hitpair, tg2[i], ng2[i], assignment) = ClusterOnce(tg2[i], ng2[i], hierarchyassignment2[i][p].copy())
                hierarchyassignment2[i][p+1] = assignment
                p = p+1
        
        # roll last array so the most recently clustered rows line up
        last_elemwidth = np.shape(hierarchyassignment2[-1])[1]
        if last_elemwidth < sublen:
            if last_elemwidth <= intermediatestop:
                rollby = sublen - intermediatestop
            if last_elemwidth > intermediatestop:
                rollby = sublen-last_elemwidth
            hierarchyassignment2[-1] = np.roll(hierarchyassignment2[-1], rollby, axis=0)
    
        # recombine the sub-lengths
        hierarchyassignment = np.concatenate(hierarchyassignment2, axis=1)
        TG = np.concatenate(tg2)
        NG = np.concatenate(ng2)
        p = sublen - intermediatestop
        
    else:
        p=0
    
    # =============================================================================
    # Cluster all the way down
    # =============================================================================
    while len(NG) > 1:  
        (hitpair, TG, NG, assignment) = ClusterOnce(TG, NG, hierarchyassignment[p].copy())
        hierarchyassignment[p+1] = assignment
        p = p+1
    
    lastiteration = np.where(np.sum(hierarchyassignment, axis=1)==0)[0][0]
    hierarchyassignment = hierarchyassignment[:lastiteration+1]
    
    if plot:
        fig, ax = plt.subplots()
        ax.imshow(hierarchyassignment)
        ax.set_xlabel('changepoint segment')
        ax.set_ylabel('clustering iteration')
        ax.set_title('Initial clustering - Fast')
    
    hierarchyassignment = np.flipud(hierarchyassignment)
    return hierarchyassignment

def ExpectationMaximizationClustering(TG, NG, mgroups, initialassignment, printupdates=True):
    '''
    This function takes a list of segment times and counts and their inital assignment/clustering (the latter is created by InitialClusteringFast() or InitialClustering())
    It follows the clustering as described by Watkins & Yang (J. Phys. Chem. B 2005, 109, 617-628)
    
    TG    : lengths of the jump segments
    NG    : number of counts in the jump segments
    initialassignment: the initial clustering
    '''
    Ttotal=np.sum(TG);
    ncp=len(TG);
    ntotal=np.sum(NG);
    
    pmj = np.zeros([mgroups, len(NG)])
    # initialize - eqn 14
    
    for m in range(mgroups):
        pmj[m] = 1.0*(initialassignment==m)
    
    success=0
    qq=0
    qqmax=200 # max nr of iterations
    logLKlast=1
    tolerance = 1E-8
    
    '''start of Figure 9, Watkins & Yang'''
    while success==0:
        qq = qq+1
        
        # M-step
        Tm = np.matmul(pmj, TG) # list of total durations of segments
        Mn = np.matmul(pmj, NG) # list of total counts in segments
        Im = Mn/Tm # list of intensities in segments
        pm = Tm/Ttotal # probability of drawing an Ij number from the mth intensity level
        pmjnew = np.zeros([mgroups, len(NG)])
        
    
        # E-step
        for m in range(mgroups):
            pmjnew[m] = pm[m]*ApproxPoissPDF(NG, TG*Im[m])[0]
    
        # check for tricky ones so you don't divide by zero
        idx = np.where(np.sum(pmjnew,axis=0)==0.) # find columns that add to zero
        
        # divide out the nontricky ones, dummy operation for the tricky ones 
        denom=np.matlib.repmat(np.sum(pmjnew, axis=0),mgroups, 1);
        denom[:,idx]=denom[:,idx]+1.0
        pmjnew=pmjnew/denom
 
#        original version, causes and catches warning:
#        pmjnew=pmjnew/np.matlib.repmat(np.sum(pmjnew, axis=0),mgroups, 1)
         # now we could have a bunch of NaN ones. We set those to the previous iteration's value
#        pmjnew[:,idx] = pmj[:,idx]
        
        '''stop criterion
        from J Phys Chem B 123, 689: If Wastkins Eqn 15 of Watkins hardly changes, you have converged'''
        logLK = np.zeros(np.shape(pmj))
        for m in range(mgroups):
            [g, logg] = ApproxPoissPDF(NG, TG*Im[m])
            logLK[m] = np.log(pm[m]) + logg # not pmjnew?
        logLK = np.sum(pmj*logLK)
        
        if qq>2:
            if abs(logLK/logLKlast-1)<tolerance:
                success = 1
        
        if qq==qqmax:
            success = 1
            if printupdates:
                print('terminated before full convergence at m='+str(mgroups))
        
        logLKlast = logLK
        pmj = pmjnew
        '''End of Figure 9'''
        
    # now sort all the levels in order of intensity
    if mgroups > 1:      
        inds = (-pm).argsort()
        pm = pm[inds]
        Tm = Tm[inds]
        Im = Im[inds]
        pmj = pmj[inds]
        
        # since you have m groups, the most likely trajectory is
        # pmjlikely = np.amax(pmj, axis=0)
        mlikely = np.argmax(pmj, axis=0)
    
    # if you only have one level
    else:
        mlikely = np.zeros(len(NG), dtype=int)
        # pmjlikely = pmj
    
    # the log likelyhood of the data given this trace is
    # Eq 15 of Watkins & Yang
    [g, logg] = ApproxPoissPDF(NG, TG*Im[mlikely]) 
    logLKv1 = np.log(pm[mlikely])+logg
    
    #alternatively take Eq 15 to the letter
    logLK = np.zeros([mgroups,len(NG)])
    for m in range(mgroups):
        [g, logg] = ApproxPoissPDF(NG, TG*Im[m])
        logLK[m]  = np.log(pm[m]) + logg
    logLKv2 = np.sum(pmj*logLK)
    
    
    # Now use this as input for eqn 16
    # to accomodate the weird sentence 8 lines below Eq 16, calculate the reduced number of CPs
    ncp = 1+len(np.where(np.diff(mlikely)!=0)[0])
    schwarz1 = np.sum(logLKv1)*2 - (2*mgroups-1)*np.log(ncp) - ncp*np.log(ntotal)
    schwarz2 = np.sum(logLKv2)*2 - (2*mgroups-1)*np.log(ncp) - ncp*np.log(ntotal)
    
    schwartz = [schwarz1, schwarz2]
    
    return (pm, Im, Tm, mlikely, schwartz, pmj)


def InitialClusteringSlow(TG, NG, plot=False):
    '''
    This function takes a list of segment times and counts, and clusters them together sequentially until there is only one state
    The output is used as an initial guess for the clustering defiend by Watkins & Yang (J. Phys. Chem. B 2005, 109, 617-628)
    
    TG    : lengths of the jump segments
    NG    : number of counts in the jump segments
    '''
    TG = np.random.random(30)
    NG = np.random.random(30)
    
    hierarchyassignment = np.zeros([len(NG), len(NG)], dtype=int)
    hierarchyassignment[0] = np.arange(len(NG))
    
    
    p = 0
    while len(NG) > 1:  # do the full clustering
        (hitpair, TG, NG, assignment) = ClusterOnce(TG, NG, hierarchyassignment[p].copy())
        hierarchyassignment[p+1] = assignment
        
        p = p+1
  
    if plot:
        fig, ax = plt.subplots()
        ax.imshow(hierarchyassignment)
        ax.set_xlabel('changepoint segment')
        ax.set_ylabel('clustering iteration')
        ax.set_title('Initial clustering - Slow')
        
    hierarchyassignment = np.flipud(hierarchyassignment)
    return hierarchyassignment
