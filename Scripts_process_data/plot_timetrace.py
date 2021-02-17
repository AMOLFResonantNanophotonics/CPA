'''
Postprocessing function that plots CPA segmented data.

CPA segmentation has resulted in lists of:
    Jumptimes T1, T2, T3, ...
    Counts per interval C1, C2, C3, equiv to intensities I_n = C_n/T_n
    Decay rates  gamma_1, gamma_2, ...

To plot time traces, the parquet raw timestamps are read in and histogrammed
to obtain conventional histogrammed instantaneous intensities. This is 
overplotted with the CPA-segmented values.

The plot routine also adds occurrence histograms.
    


'''

def plot_timetrace(PlotCombinedTT, Show_intermediateplots, CPA_insteadof_binned, Simulated_insteadof_MeasuredData):
    
    import numpy as np
    import os
    import math
    
    import matplotlib.pyplot as plt
    from matplotlib import gridspec as gsp
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
    
    
    
    Dotlist = [i for i in os.listdir(pre.timetags_filepath) if i.startswith('Dot_') and pre.sig in i]
    # =============================================================================
    # start processing the data
    # =============================================================================
    print('\n\nRunning routine to plot time trace (binned) alongside CPA Segmentation')  
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
        # Load the timestamps,  purely for plotting traditionally binned data. CPA results are retrieved separately
        # =============================================================================
        '''
        timestamps_chX_bin : all events in channel X E (A, B, R)
        timestamps_bin     : all events in channels A and B, chronologically
        '''
    
        timestamps_chA_bin, timestamps_chB_bin, timestamps_chR_bin = load.LoadTimeStamps(pre.timetags_filepath+'Dot_%02d/' %dot_idx, pre.timetags_filenames, pre.timetags_headers)
        timestamps_bin = np.sort(np.concatenate((timestamps_chA_bin, timestamps_chB_bin)))
    
        
    
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
        # Prepare plottable time traces  
        # =============================================================================
        '''
        Makes a visualization of the I, gamma behavior of the emitter in time
    
        input
        seg_times_bin_extra : array of the time bin of the found changepoint. Unit is counter card resolution. Includes 0 and the end of the measurement
        gamma_lst         : 1D array of the single-exponential decay rates of the CPA segments
        ccmeta            : the parameters of the correlation calculation and the metadata
        BinLength_s       : some choice for the bin length. This is for visualization ONLY. Unit is s
        timestamps_bin    : list of photon events measured, in units of the counter card resolution
        segmentcounts     : 1D array, number of photon events in each CPA segment
        segmentlengths_s  : 1D array of the time urations of the CPA segments
        maxdecayrate      : free to choose. The maximum decay rate that you consider
        pre.dtau_ns           : the timng resolution of the counter card, in ns
    
        output
        retrievedx_TT  : duplicated CPA jump times, in unit of counter card resolution
        retrievedy_TTI : duplicated intensites, segmented by CPA. They are normalized to Binlength_s
        retrievedy_TTg : duplicated decay rates, segmented by CPA
        TTI_hist_x     : bins for the histogram of the intensities
        TTI_Binhist    : histogram of the intensities, data binned in bins as wide as BinLength_s
        TTI_CPAhist    : histogram of the intensities, found by CPA
        TTg_hist_x     : the values of the bins used for the lifetime histogram
        TTg_CPAhist     : histogram of the decay rates as found by CPA
        '''
    
    
        BinLength_s = 2e-3 # in seconds
        maxT = timestamps_bin[-1]*pre.dtau_ns*1e-9
        Nbins = math.floor(maxT/BinLength_s)
        
            
        retrievedx_TT = np.repeat(seg_times_bin_extra,2)[1:-1]
        retrievedy_TTg = np.repeat(fitdata['fit_gamma'],2)

        retrieved_norm_intensity = segmentcounts/segmentlengths_s*BinLength_s
        retrievedy_TTI = np.repeat(retrieved_norm_intensity,2)

        # plot the timetrace
        if PlotCombinedTT:
            ## automatically generate plot bounds
            ## if you are not happy with them:  hardcode your own override
            
            max_I = 5*np.mean(retrievedy_TTI)
            maxdecayrate=np.mean(4*retrievedy_TTg)    
            
            ## some arbitrary plot range 
     
            ## e.g., zoom in on 1 sec or so.
            
            tmin=16.1
            tmax=16.8
            
#            tmin=23
#            tmax=24


            #tmin=40.1
            #tmax=41.1


            print('\ntime axis of time trace arbitrarily zoomed in, for esthetic plot')
            print('histogram is over full time span however\n ')
            
            figTT = plt.figure(figsize=(10/2.51, 7/2.51))
            gs0    = gsp.GridSpec(2, 2, height_ratios=[1,1], width_ratios=[1,0.3], wspace=0, hspace=0.1)
    
            axTT_I = figTT.add_subplot(gs0[0,0])
            axTT_I_h = figTT.add_subplot(gs0[0,1], sharey=axTT_I)
            axTT_L = figTT.add_subplot(gs0[1,0], sharex=axTT_I)
            axTT_L_h = figTT.add_subplot(gs0[1,1], sharey=axTT_L)
    
            #plot a binned timetrace of the photon events
            # Remember that this binning has _no_ impact on the rest of the calculations
            axTT_I.plot(np.linspace(0, timestamps_bin[-1]*(pre.dtau_ns), num=Nbins)*1e-9, np.histogram(timestamps_bin,Nbins)[0], label='data', linewidth=1.2, color='C0')
            
            # Intensities
            # time trace
            axTT_I.plot(retrievedx_TT*pre.dtau_ns*1e-9, retrievedy_TTI, label='retrieved levels', color='C2')
            axTT_I.set_ylabel('cts/'+str(BinLength_s*1e3)+' ms')
            axTT_I.tick_params(axis='both', which='major', labelbottom=False)
            axTT_I.set_ylim((-1,max_I))
            axTT_I.set_xlim((tmin,tmax))
            axTT_I.legend(loc=1)
    
            # Intensity histogram
            numbin=20
            axTT_I_h.hist(retrieved_norm_intensity, bins=np.linspace(0,max_I, numbin), color='C2', orientation='horizontal')
            axTT_I_h.tick_params(axis='both', which='major', labelleft=False)
            uplim=len(retrieved_norm_intensity)/numbin*5
            lowlim=-uplim/30
            axTT_I_h.set_xlim(lowlim,uplim)
            axTT_I_h.set_xticks([])
    
            # Lifetimes
            # time trace
            
           
    
            axTT_L.plot(retrievedx_TT*pre.dtau_ns*1e-9, retrievedy_TTg, color='C2')
            axTT_L.set_xlabel('measurement time (s)')
            axTT_L.set_ylabel('$\gamma$ (ns$^{-1}$)')
            axTT_L.set_ylim((0,1.1*maxdecayrate))
    
    
            # Lifetime histogram        
            numbin=20
            axTT_L_h.hist(fitdata['fit_gamma'], bins=np.linspace(0,1.1*maxdecayrate, numbin), color='C2', orientation='horizontal')
            axTT_L_h.tick_params(axis='both', which='major', labelleft=False)
            uplim=len(fitdata['fit_gamma'])/numbin*5
            lowlim=-uplim/30
            axTT_I_h.set_xlim(lowlim,uplim)
            axTT_L_h.set_xticks([])
    
            figTT.savefig(savepath +'TCSPC_timetrace.png')
            # figTT.savefig(savepath_TT +dotID+'.png')
    
            if not Show_intermediateplots and pre.sig=="_'":
                plt.close(figTT)
    
