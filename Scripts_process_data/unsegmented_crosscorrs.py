''' 
Reads in timestamps from parquet files
Constructs fluorescence decay trace of entire data set
[Uses correlate-jit]

Fits double exponential decay, assuming the Poisson maximum likelyhood estimator
[uses fit_MLE]
'''


def unsegmented_crosscorrs(MakeTotalDecayCurve, PlotTotalDecayCurve, Show_intermediateplots, CPA_insteadof_binned, Simulated_insteadof_MeasuredData):
    
    import numpy as np
    import os
    # import pyarrow.parquet as pq
    
    
    # import matplotlib as mpl
    import matplotlib.pyplot as plt
    
    import fitexp_MLE as fitexp
    import correlate_jit
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
        outputfolder = pre.outputfolder_1 + pre. outputfolder_2_CPA
    else:
        outputfolder = pre.outputfolder_1 + pre.outputfolder_2_binned
     
    
    
    # =============================================================================
    # start processing the data
    # =============================================================================
    
    
    Dotlist = [i for i in os.listdir(pre.timetags_filepath) if i.startswith('Dot_') and pre.sig in i]
    print('\n\nRunning routine to construct, plot  and fit decay trace from full data set')
    for dot_file in Dotlist:
        dot_idx = int(dot_file[4:6])
        print('###############################################')
        print('Starting Dot', dot_idx)
    
    
        # =============================================================================
        # create the folder to save the data
        # =============================================================================
        savepath = outputfolder + 'Dot_%02d/' %dot_idx
        if not os.path.exists(savepath):
            os.makedirs(savepath)
    
    
        # =============================================================================
        # Load the timestamps
        # =============================================================================
        '''
        timestamps_chX_bin : all events in channel X E (A, B, R)
        timestamps_bin     : all events in channels A and B, chronologically
        '''
    
        timestamps_chA_bin, timestamps_chB_bin, timestamps_chR_bin = load.LoadTimeStamps(pre.timetags_filepath+'Dot_%02d/' %dot_idx, pre.timetags_filenames, pre.timetags_headers)
        # Tau0A_bin, Tau0B_bin= fitexp.Quickdetectorzerotime(timestamps_chA_bin,timestamps_chB_bin, timestamps_chR_bin,int(pre.minitimebin_ns/pre.dtau_ns))
        
        
        # =============================================================================
        # Load metadata
        # =============================================================================
        file_metadata = load.ReadDict(savepath+pre.metadata_filename, 20)
    
        # =============================================================================
        # The lifetime of the system, non-segmented
        # =============================================================================
        '''
        Calculate the lifetime curve of the entire measurement, with no segmentation
        Note that this can also be done by summing the segmened curves from either channel
    
        input
        crosscorrs    : the calculated correlations, segmented by CPA
        cctaus_ns     : the shift times used for the correlation, in ns
        pre.mygap_bin         : a time interval to exclude from the fitting procedure, to remove the effects of electronic ringing
        rise_time_bin : the number of bins at the end of the trace to exclude from the fitting procedure. This accounts for the rise time of the measurement apparatus
    
        output
        crosscorrsAB_flat  : the 1D decay trace of the entire mresurement
        fit : the fit parameters of an n-expoential decay to the correlation.
        '''
    
        filepath_totdecay = savepath + pre.totaldecay_filename
    
        if MakeTotalDecayCurve:
            cctaus_bin =  np.arange(0,int(pre.minitimebin_ns/pre.dtau_ns)).astype(np.int64)
            cctaus_ns  = cctaus_bin* pre.dtau_ns
    
            # make the correlation
            print('\nRun Wahl correlator to obtain decay histogram')
            crosscorrsA = correlate_jit.Corr(timestamps_chR_bin, timestamps_chA_bin, cctaus_bin)
            crosscorrsB = correlate_jit.Corr(timestamps_chR_bin, timestamps_chB_bin, cctaus_bin)
    
            # prepare the data for fitting
            crosscorrsA_prep, cctaus_ns_prep = fitexp.PrepCrossCorrs([crosscorrsA], cctaus_ns, -file_metadata['Tau0A_bin'], gap=pre.mygap_bin, rise_time_bin=int(5/pre.dtau_ns))
            crosscorrsB_prep, cctaus_ns_prep = fitexp.PrepCrossCorrs([crosscorrsB], cctaus_ns, -file_metadata['Tau0B_bin'], gap=pre.mygap_bin, rise_time_bin=int(5/pre.dtau_ns))
            
            crosscorrsA_prep = crosscorrsA_prep[0]
            crosscorrsB_prep = crosscorrsB_prep[0]
            
            crosscorrsAB_prep = crosscorrsA_prep + crosscorrsB_prep
            
            # fit the decay curve
            fit = fitexp.minimize_func(cctaus_ns_prep, crosscorrsAB_prep, decay='double', decayguess=0.5, decayguess2=0.2)
    
            load.SaveUnsegmentedCorr(filepath_totdecay, fit, cctaus_ns_prep, crosscorrsAB_prep)
    
        if not MakeTotalDecayCurve:
            
            fit, cctaus_ns_prep, crosscorrsAB_prep = load.LoadUnsegmentedCorr(filepath_totdecay)
            
            
    
        if PlotTotalDecayCurve:
            figlt, axlt = plt.subplots(figsize=(10/2.51, 10/2.51))
            axlt.scatter(cctaus_ns_prep, crosscorrsAB_prep, c='C0')
            axlt.plot(cctaus_ns, fitexp.exp_2(fit['fit_result'], cctaus_ns), c='C2', linewidth=3)
            axlt.set_xlabel('\u03C4 (ns)')
            axlt.set_ylabel('counts')
            axlt.set_title('total emitter lifetime')
            axlt.set_yscale('log')
            axlt.set_ylim( (crosscorrsAB_prep[-1]*0.9, max(crosscorrsAB_prep)/0.9) )
            figlt.savefig(savepath+'unsegmented_correlation.png')
    
            if not Show_intermediateplots and pre.sig == '_':
                plt.close(figlt)
    
    
    
    
    












