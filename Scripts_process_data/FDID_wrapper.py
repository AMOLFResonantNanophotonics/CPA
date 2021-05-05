''' 
Function that reads in segmented data from CPA, and plots FDID diagrams
!Take note not to confuse FDID and FLID diagrams

Premise:  CPA has converted time stamps into segmented data, 
consisting of lists with
(A)  the switching instances T1, T2, T3, T4
(B)  the counts in each interval In starting at Tn. This converts to intensity
(C)  the decay rates fitted by a POisson statistics based MLE method

The output is converted into 2D  diagrams histogramming 
fluorescence decay rates and intensities with the following premises:
    
(A)  each entry (In, gamma_n) provides not a unique discrete histogram click, 
but a Gauss with width given by the estimated error in I_n and gamma_n

(B)  the weight of each entry [area under gauss] allows three choices,  i.e.,  
by duration T_n , by total photon contribution I_n*T_n, 
or simply unit weight per segment.  Note that weighting by duration is most 
like conventional FDIDs
    
'''
def FDID_wrapper(MakeFDIDs, PlotFDIDs, MakeFLIDs, PlotFLIDs, Show_intermediateplots, CPA_insteadof_binned, Simulated_insteadof_MeasuredData):
    
    import numpy as np
    import os 
    
    
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib import gridspec as gsp
     
    import FDID_functions as FDID 
     
    
    import acc_functions as acc
    
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
    # start processing the data
    # =============================================================================
    
    Dotlist = [i for i in os.listdir(pre.timetags_filepath) if i.startswith('Dot_') and pre.sig in i]
    print('\n\nRunning routine to plot FDIDs')
    
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
        # Read in fitted lifetimes
        # =============================================================================
        '''
        fitdata contains 4 parameters of the fitted single-exponential decay
        fit_A         : the scaling factor
        fit_gamma     : the decay rate, in ns^-1
        fit_gamma_err : the error on the fitted decay rate
        fit_bg        : the background noise
        '''
     
        filepath_fitdata = savepath + pre.segment_fit_filename
        ccmeta, fitdata = load.LoadFitData(filepath_fitdata)

        # =============================================================================
        # Turn the FDID data into FLID data
        # =============================================================================
        
        fitdata['fit_tau'] = 1/fitdata['fit_gamma'] # unit is ns
        fitdata['fit_tau_err'] = fitdata['fit_gamma_err']/fitdata['fit_gamma']/fitdata['fit_gamma']
        
        
        # =============================================================================
        # Now make the FDID diagrams
        # Based on Decay rates
        # =============================================================================
        '''
        Make the FDID diagrams, a 2D histogram of the CPA segments
    
        input
        fitdata           : dataframe  of the single-exponential decay rates of the CPA segments, their errors, the lifetimes & errors, background, and scaling parameters
        segmentlengths_s  : 1D array of the time urations of the CPA segments
        segmentcounts     : 1D array, number of photon events in each CPA segment
    
        output
        FDID_prio_none : 2D array, the FDID where all CPA segments have the same weight
        FDID_prio_cts  : 2D array, the FDID where the CPA segments are weighted by the amount of counts
        FDID_prio_dur  : 2D array, the FDID where the CPA segments are weighted by their duration
        FDIDextent     : tuple, the (x1, x2, y1, y2) extent of the FDID. Used for plotting
        '''
    
        filepath_FDID = savepath + pre.FDID_filename
        
        
        if MakeFDIDs:
            # Equal Priority
            FDID_prio_none, FDIDextent = FDID.MakeFDID(fitdata['fit_gamma'], fitdata['fit_gamma_err'], segmentcounts/segmentlengths_s, np.sqrt(segmentcounts)/segmentlengths_s, xmax = pre.maxdecayrate)
            # Count Priority
            FDID_prio_cts, FDIDextent  = FDID.MakeFDID(fitdata['fit_gamma'], fitdata['fit_gamma_err'], segmentcounts/segmentlengths_s, np.sqrt(segmentcounts)/segmentlengths_s, countPriority=True, xmax = pre.maxdecayrate)
            # Duration Priority
            FDID_prio_dur, FDIDextent  = FDID.MakeFDID(fitdata['fit_gamma'], fitdata['fit_gamma_err'], segmentcounts/segmentlengths_s, np.sqrt(segmentcounts)/segmentlengths_s, CPAlengthpriority=True, difftimelst=segmentlengths_s, xmax = pre.maxdecayrate)
    
            print('Writing FDIDs to file as 2D matrices in ')
            print(filepath_FDID)
            
            with open(filepath_FDID, 'w') as f:
                print('Extent (gamma_min, gamma_max, I_min, I_max): '+ str(FDIDextent), file=f)
    
                print('\nNo Priority', file=f)
                np.savetxt(f, FDID_prio_none)
    
                print('\nCount Priority', file=f)
                np.savetxt(f, FDID_prio_cts)
    
                print('\nDuration Priority', file=f)
                np.savetxt(f, FDID_prio_dur)
    
        if not MakeFDIDs:
            print('Attempting to read FDIDs from file for replotting: ')
            print(filepath_FDID)
            with open(filepath_FDID, 'r') as f:
                FDIDdata = f.read().split('\n\n')
                FDIDextent = FDIDdata[0].split(': ')[1]
                FDIDextent = tuple(map(float, FDIDextent[1:-1].split(', ')))
    
            breaks = []
            with open(filepath_FDID, 'r') as f:
                for idx,line in enumerate(f):
                    if line == '\n':
                        breaks += [idx]
            breaks += [idx] # add length of file
            FDID_prio_none = np.genfromtxt(filepath_FDID, skip_header=breaks[0]+2, max_rows=breaks[1]-breaks[0]-2)
            FDID_prio_cts = np.genfromtxt(filepath_FDID, skip_header=breaks[1]+2, max_rows=breaks[2]-breaks[1]-2)
            FDID_prio_dur = np.genfromtxt(filepath_FDID, skip_header=breaks[2]+2)


        if PlotFDIDs:
            figFDIDall = plt.figure(figsize=(11/2.51,4/2.51))
            gs0 = gsp.GridSpec(1, 4, width_ratios=[1,1,1,0.1], wspace=0.1)
            axFDID1 = figFDIDall.add_subplot(gs0[0])
            axFDID2 = figFDIDall.add_subplot(gs0[1])
            axFDID3 = figFDIDall.add_subplot(gs0[2])
            ax_cbar = figFDIDall.add_subplot(gs0[3])
    
            im1 = axFDID1.imshow(np.flipud(FDID_prio_none/np.amax(FDID_prio_none)), extent=FDIDextent, aspect='auto', cmap=pre.mycmapFDID)
            im2 = axFDID2.imshow(np.flipud(FDID_prio_cts/np.amax(FDID_prio_cts)), extent=FDIDextent, aspect='auto', cmap=pre.mycmapFDID)
            im3 = axFDID3.imshow(np.flipud(FDID_prio_dur/np.amax(FDID_prio_dur)), extent=FDIDextent, aspect='auto', cmap=pre.mycmapFDID)
            axFDID1.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
            axFDID2.set_yticklabels([])
            axFDID3.set_yticklabels([])
            axFDID1.set_xlabel(r'$\gamma$ (ns$^{-1}$)')
            axFDID2.set_xlabel(r'$\gamma$ (ns$^{-1}$)')
            axFDID3.set_xlabel(r'$\gamma$ (ns$^{-1}$)')
    
            cbar =plt.colorbar(im3, cax=ax_cbar)
            axFDID1.set_ylabel('I (counts/s)')
            axFDID1.set_title('Equal weight')
            axFDID2.set_title('Intensity weight')
            axFDID3.set_title('Duration weight')
    
            axFDID1.set_xlim([0,pre.maxdecayrate])
            axFDID2.set_xlim([0,pre.maxdecayrate])
            axFDID3.set_xlim([0,pre.maxdecayrate])
    
            figFDIDall.savefig(savepath+'FDID_all.png')
    
            if not Show_intermediateplots and pre.sig == '_':
                plt.close(figFDIDall)
    


                
        # =============================================================================
        # Now make the FLID diagrams
        # Based on Lifetimes!
        # =============================================================================
        '''
        Make the FLID diagrams, a 2D histogram of the CPA segments
    
        input
        fitdata           : dataframe  of the single-exponential decay rates of the CPA segments, their errors, the lifetimes & errors, background, and scaling parameters
        segmentlengths_s  : 1D array of the time urations of the CPA segments
        segmentcounts     : 1D array, number of photon events in each CPA segment
    
        output
        FLID_prio_none : 2D array, the FLID where all CPA segments have the same weight
        FLID_prio_cts  : 2D array, the FLID where the CPA segments are weighted by the amount of counts
        FLID_prio_dur  : 2D array, the FLID where the CPA segments are weighted by their duration
        FLIDextent     : tuple, the (x1, x2, y1, y2) extent of the FLID. Used for plotting
        '''           

    
        filepath_FLID = savepath + pre.FLID_filename
        if MakeFLIDs:   
            # Equal Priority
            FLID_prio_none, FLIDextent = FDID.MakeFDID(fitdata['fit_tau'], fitdata['fit_tau_err'], segmentcounts/segmentlengths_s, np.sqrt(segmentcounts)/segmentlengths_s, xmax = pre.maxlifteime)
            # Count Priority
            FLID_prio_cts, FLIDextent  = FDID.MakeFDID(fitdata['fit_tau'], fitdata['fit_tau_err'], segmentcounts/segmentlengths_s, np.sqrt(segmentcounts)/segmentlengths_s, countPriority=True, xmax = pre.maxlifteime)
            # Duration Priority
            FLID_prio_dur, FLIDextent  = FDID.MakeFDID(fitdata['fit_tau'], fitdata['fit_tau_err'], segmentcounts/segmentlengths_s, np.sqrt(segmentcounts)/segmentlengths_s, CPAlengthpriority=True, difftimelst=segmentlengths_s, xmax = pre.maxlifteime)
    
            print('Writing FLIDs to file as 2D matrices in ')
            print(filepath_FLID)
            
            with open(filepath_FLID, 'w') as f:
                print('Extent (gamma_min, gamma_max, I_min, I_max): '+ str(FLIDextent), file=f)
    
                print('\nNo Priority', file=f)
                np.savetxt(f, FLID_prio_none)
    
                print('\nCount Priority', file=f)
                np.savetxt(f, FLID_prio_cts)
    
                print('\nDuration Priority', file=f)
                np.savetxt(f, FLID_prio_dur)

    
        if not MakeFLIDs:
            print('Attempting to read FLIDs from file for replotting: ')
            print(filepath_FLID)
            with open(filepath_FLID, 'r') as f:
                FLIDdata = f.read().split('\n\n')
                FLIDextent = FLIDdata[0].split(': ')[1]
                FLIDextent = tuple(map(float, FLIDextent[1:-1].split(', ')))
    
            breaks = []
            with open(filepath_FLID, 'r') as f:
                for idx,line in enumerate(f):
                    if line == '\n':
                        breaks += [idx]
            breaks += [idx] # add length of file
            FLID_prio_none = np.genfromtxt(filepath_FLID, skip_header=breaks[0]+2, max_rows=breaks[1]-breaks[0]-2)
            FLID_prio_cts = np.genfromtxt(filepath_FLID, skip_header=breaks[1]+2, max_rows=breaks[2]-breaks[1]-2)
            FLID_prio_dur = np.genfromtxt(filepath_FLID, skip_header=breaks[2]+2)
    
        if PlotFLIDs:
            figFLIDall = plt.figure(figsize=(11/2.51,4/2.51))
            gs0 = gsp.GridSpec(1, 4, width_ratios=[1,1,1,0.1], wspace=0.1)
            axFLID1 = figFLIDall.add_subplot(gs0[0])
            axFLID2 = figFLIDall.add_subplot(gs0[1])
            axFLID3 = figFLIDall.add_subplot(gs0[2])
            ax_cbar = figFLIDall.add_subplot(gs0[3])
    
            im1 = axFLID1.imshow(np.flipud(FLID_prio_none/np.amax(FLID_prio_none)), extent=FLIDextent, aspect='auto', cmap=pre.mycmapFDID)
            im2 = axFLID2.imshow(np.flipud(FLID_prio_cts/np.amax(FLID_prio_cts)), extent=FLIDextent, aspect='auto', cmap=pre.mycmapFDID)
            im3 = axFLID3.imshow(np.flipud(FLID_prio_dur/np.amax(FLID_prio_dur)), extent=FLIDextent, aspect='auto', cmap=pre.mycmapFDID)
            axFLID1.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
            axFLID2.set_yticklabels([])
            axFLID3.set_yticklabels([])
            axFLID1.set_xlabel(r'$\tau$ (ns)')
            axFLID2.set_xlabel(r'$\tau$ (ns)')
            axFLID3.set_xlabel(r'$\tau$ (ns)')
    
            cbar =plt.colorbar(im3, cax=ax_cbar)
            axFLID1.set_ylabel('I (counts/s)')
            axFLID1.set_title('Equal weight')
            axFLID2.set_title('Intensity weight')
            axFLID3.set_title('Duration weight')
    
            axFLID1.set_xlim([0,pre.maxlifteime])
            axFLID2.set_xlim([0,pre.maxlifteime])
            axFLID3.set_xlim([0,pre.maxlifteime])
    
            figFLIDall.savefig(savepath+'FLID_all.png')
    
            if not Show_intermediateplots and pre.sig == '_':
                plt.close(figFLIDall)
    
