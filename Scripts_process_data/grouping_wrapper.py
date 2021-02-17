'''
Wrapper routine for grouping.

The routine takes as input CPA segmented data, meaning a list of 
- Switching moments T1,T2,T3,...
- Counts in each  C1, C2, C3,  which is equiv. to intensities I_n=C_n/T_n

The result is a list of groupings in m levels, with for each m
- The most likely assignment of intervals n to groups m
- The concomitant intensities for each of the levels
-The Bayesian Information criterion

The routine processes these to return plots reporting on 
- The BIC versus m
- The occupancy of m levels versus m
    
The actual work is actually outssources to grouping_functions.py
'''

def grouping_wrapper(MakeGrouping, PlotBIC, PlotOccupation, Show_intermediateplots, CPA_insteadof_binned, Simulated_insteadof_MeasuredData):
    
    import numpy as np
    import os
    # import math
    # from scipy.optimize import minimize
    import pandas as pd
    # import scipy.optimize
    import pandas
    
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
     
    
    import acc_functions as acc
    import grouping_functions as grouping
    import loaddata_functions as load
    
       
    # =============================================================================
    # import preamble
    # =============================================================================
    
    if Simulated_insteadof_MeasuredData:
        import preamble_simulated as pre
    else:
        import preamble_measured as pre

    # =============================================================================
    # set outputfolder
    # =============================================================================
    
    if CPA_insteadof_binned:
        outputfolder = pre.outputfolder_1 + pre.outputfolder_2_CPA
    else:
        outputfolder = pre.outputfolder_1 + pre.outputfolder_2_binned
    
    
    # =============================================================================
    # get to work on the specified dots
    # =============================================================================
    
    
    
    Dotlist = [i for i in os.listdir(pre.timetags_filepath) if i.startswith('Dot_') and pre.sig in i]
    
    # =============================================================================
    # start processing the data
    # =============================================================================
    
   
    print('\n\nRunning grouping/level clustering routine')     
    for dot_file in Dotlist:
        dot_idx = int(dot_file[4:6])
        print('##################################################')
        print('Starting Dot', dot_idx)
    
    
        # =============================================================================
        # create the folder to save the data
        # =============================================================================
        savepath = outputfolder + 'Dot_%02d/' %dot_idx
        if not os.path.exists(savepath):
            os.makedirs(savepath)
     
    
        # =============================================================================
        # Segment bounds (changepoints or bins)
        # =============================================================================
        '''
        seg_idx_extra   : array of the indices of the photons at the found changepoints, Includes 0, and the total photon count
        seg_times_bin_extra : array of the time bin of the found changepoint. Unit is counter card resolution. Includes 0 and the end of the measurement
        segmentlengths_s    : array of the length of the segments, in seconds
        segmentcounts       : array of the number of counts in every segmens
        '''
        
        filepath_segbounds = savepath + pre.segments_filename
        seg_idx_extra, seg_times_bin_extra, segmentlengths,segmentcounts=load.LoadSegments(filepath_segbounds)    
        segmentlengths_s=segmentlengths*pre.dtau_ns*1e-9       # in seconds  [durations]
    
    
        # =============================================================================
        # Bayesian Information Criterion (BIC) and grouping
        # =============================================================================
        '''
        Calculates the likelihood of a a certain number of states to be the best description of the data
    
        input
        mingroups           : smallest nr of groups to test for. Best to set to 1
        maxgroups           : largest nr of groups to test for.
        dtau_ns             : timing resolution of the counter card
        seg_idx_extra       : array of the indices of the photons at the found changepoints, Includes 0, and the total photon count
        seg_times_bin_extra : array of the time bin of the found changepoint. Unit is counter card resolution. Includes 0 and the end of the measurement
    
        output
        BICdata     : dictionary of nr of levels tested for, and their corresponding likelihoods
        mlikely_tot : most likely trajectory though the states
        psummarry   : the likelihood of drawing a certain state, given m states. This is in TIME occupancy, not NUMBER occupancy (see definition of pm)
        '''
        
        filepath_grouping = savepath + pre.grouping_filename
    
        if MakeGrouping:    ## do actual grouping work
            mingroups = 1
            maxgroups = 20
    
            BICdata, mlikely_tot, psummarry = grouping.Grouping(segmentlengths_s, segmentcounts, mingroups, maxgroups)
            
            print('Saving outcome of grouping to file:')
            print(filepath_grouping)
            with open(filepath_grouping, 'w') as f:
                np.savetxt(f, BICdata, header=' '.join(list(BICdata)))
    
                print('\ngiven m states, what is the likelihood of drawing a certain state?', file=f)
                print('Here shown for '+ str(mingroups) +' through '+ str(maxgroups) +' states', file=f)
                print('psummarry', file=f)
                np.savetxt(f, psummarry)
    
                print('\nGiven m states, what is the most likely trajectory of the emitter through these states.', file=f)
                print('Here shown for '+ str(mingroups) +' through '+ str(maxgroups) +' states', file=f)
                print('mlikely_tot', file=f)
                np.savetxt(f, mlikely_tot)
    
        if not MakeGrouping:
            breaks = []
            print('Attempting to load a previously obtained grouping from file:')
            print(filepath_grouping)
    
            with open(filepath_grouping, 'r') as f:
                for idx,line in enumerate(f):
                    if idx==0:
                        BICheaders = line[2:-1].split(' ')
                    if line == '\n':
                        breaks += [idx]
            breaks += [idx] # add length of file
    
            BICdata = np.genfromtxt(filepath_grouping, skip_header=1, max_rows=breaks[0]-2)
            BICdata = pd.DataFrame(data=BICdata, index=np.arange(len(BICdata)), columns=BICheaders)
            psummarry = np.genfromtxt(filepath_grouping, skip_header=breaks[0]+4, max_rows=breaks[1]-breaks[0]-4)
            mlikely_tot = np.genfromtxt(filepath_grouping, skip_header=breaks[1]+4)
  
            mingroups = min(BICdata['mlst'])
            maxgroups = max(BICdata['mlst'])
    
        if PlotBIC:
            # Plot the BIC
            fig_BIC, ax_BIC = plt.subplots(figsize=(8/2.51, 8/2.51))
            ax_BIC.scatter(BICdata['mlst'], BICdata['Schw0'], s=10)
            ax_BIC.plot(BICdata['mlst'], BICdata['Schw0'], '--', lw=1)
            ax_BIC.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
            
            ax_BIC_ins = inset_axes(ax_BIC, width='50%', height='50%', loc=4)
            ax_BIC_ins.scatter(BICdata['mlst'][5:], BICdata['Schw0'][5:], s=10)
            ax_BIC_ins.plot(BICdata['mlst'][5:], BICdata['Schw0'][5:], '--', lw=1)
            ax_BIC_ins.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
            ax_BIC.set_xlabel('$n_G$')
            ax_BIC.set_ylabel('BIC')
            ax_BIC.set_title('m-state likelihood')
            fig_BIC.savefig(savepath+'BIC.png')
            # fig_BIC.savefig(savepath_BIC+dotID+'.png')
    
            
    
        if PlotOccupation:
            # weight in each state
            fig_occupation, ax_occupancy = plt.subplots(figsize=(8/2.51,8/2.51))
            ax_occupancy.imshow(np.flipud(psummarry), extent=(mingroups-1,maxgroups,mingroups-1,maxgroups), cmap=pre.mycmap)
            myticks = np.round(np.linspace(1,maxgroups,4))
            ax_occupancy.set_xticks(myticks-0.5)
            ax_occupancy.set_xticklabels(myticks.astype(int))
            ax_occupancy.set_yticks(myticks-0.5)
            ax_occupancy.set_yticklabels(myticks.astype(int))
            ax_occupancy.set_xlabel('index of states')
            ax_occupancy.set_ylabel('nr of states available')
            ax_occupancy.set_title('weight in each state')
            fig_occupation.savefig(savepath+'grouping_occupancy.png',dpi=300)
            # fig_occupation.savefig(savepath_occupancy+dotID+'.png',dpi=300)
    
        if pre.sig=='_' and not Show_intermediateplots:
            plt.close('all')    

