'''
Wrapper routine that reads in the segment bounds

    -  construction of fluorescence decay traces per segment, and fit

'''


def Segmented_crosscorrs(MakeCrossCorrs, MakeLifetimelst, Simulated_insteadof_MeasuredData, CPA_insteadof_binned):
    
    #import matplotlib.pyplot as plt
    import numpy as np
    import os
    #import pyarrow.parquet as pq
    #import pandas as pd
    
    import loaddata_functions as load
    import CPA_functions as findjumps
    import fitexp_MLE as fitexp
    import correlate_jit
    #import acc_functions as acc
    
    

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
    
    
    print('\nRunning CPA routine, start by reading in data\n')
    
    for dot_file in Dotlist:
        dot_idx = int(dot_file[4:6])
        print('\n##################################################')
        print('Starting Dot', dot_idx)
    
    
        # =============================================================================
        # create the folder to save the data
        # =============================================================================
        savepath = outputfolder + 'Dot_%02d/' %dot_idx
        if not os.path.exists(savepath):
            os.makedirs(savepath)
        
        
        # =============================================================================
        # load the timestamps
        # =============================================================================
    
        timestamps_chA_bin, timestamps_chB_bin, timestamps_chR_bin = load.LoadTimeStamps(pre.timetags_filepath+'Dot_%02d/' %dot_idx, pre.timetags_filenames, pre.timetags_headers)
        timestamps_bin = np.sort(np.concatenate((timestamps_chA_bin, timestamps_chB_bin)))
    
        # =============================================================================
        # Load metadata
        # =============================================================================
        file_metadata = load.ReadDict(savepath+pre.metadata_filename, 20)
    
    
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
        # Segmented cross-correlation for segment-decay traces
        # =============================================================================
        '''
        correlating the photon arrival times with the laser reference times
        this data will become the decay curves
    
        input
        fname_tcspc                  : the file name of the TCSPC file
        ccparams                     : the parameters for the cross-correlation calculation
        seg_times_bin_extra : array of the time bin of the found changepoint. Unit is counter card resolution. Includes 0 and the end of the measurement
        pre.dtau_ns                      : the tiing resolution of the counter card, in ns
        ccmeta                       : the parameters of the correlation calculation and the metadata
    
        output
        crosscorrs : 2D aray. The calculated correlation, per channel
        cctaus_ns  : the shift times used for the correlation, in ns
        '''
    
        filepath_cc = savepath + pre.crosscorrs_segmented_filename
        if MakeCrossCorrs:
            print('\nCorrelating detector and laser in segments, for decay traces')
            cctaus_bin =  np.arange(0,int(pre.minitimebin_ns/pre.dtau_ns)+1).astype(np.int64)
            cctaus_ns  = cctaus_bin* pre.dtau_ns
            crosscorrs = np.zeros((2,len(seg_times_bin_extra)-1, int(pre.minitimebin_ns/pre.dtau_ns)+1))
    
            for a in range(len(seg_times_bin_extra)-1):
                # take only the timestamps inside the relevant segment
                timestamps_chA_bin_part = timestamps_chA_bin[(seg_times_bin_extra[a]<=timestamps_chA_bin)*(timestamps_chA_bin<seg_times_bin_extra[a+1])]
                timestamps_chB_bin_part = timestamps_chB_bin[(seg_times_bin_extra[a]<=timestamps_chB_bin)*(timestamps_chB_bin<seg_times_bin_extra[a+1])]
                timestamps_chR_bin_part = timestamps_chR_bin[(timestamps_chR_bin>=seg_times_bin_extra[a])*(timestamps_chR_bin<seg_times_bin_extra[a+1])]
    
                crosscorrs[0,a] = correlate_jit.Corr(timestamps_chR_bin_part, timestamps_chA_bin_part, cctaus_bin)
                crosscorrs[1,a] = correlate_jit.Corr(timestamps_chR_bin_part, timestamps_chB_bin_part, cctaus_bin)
    
            # save the data
            load.SaveSegmentedCC(filepath_cc, crosscorrs, cctaus_ns)
    
        if not MakeCrossCorrs:
            # load the data
            crosscorrs, cctaus_ns = load.LoadSegmentedCC(filepath_cc)
    
    
        # =============================================================================
        # Segmented intensity lists, and fitted lifetimes/decay rates
        # =============================================================================
        '''
        Create a list of the lifetimes and intensities of the segmented data
    
        input
        crosscorrs                   : the calculated correlations, segmented by CPA
        cctaus_ns                    : the shift times used for the correlation, in ns
        mygap                        : a time interval to exclude from the fitting procedure, to remove the effects of electronic ringing
        rise_time_bin                : the number of bins at the end of the trace to exclude from the fitting procedure. This accounts for the rise time of the measurement apparatus
    
    
        output
        crosscorrs_merged : 1D array, correlation of both channels with the reference pulses. Corrected for possible differences in cable length between detectors
        gamma_lst         : 1D array of the single-exponential decay rates of the CPA segments
        gamma_err_lst     : 1D array of the errors of the single-exponential decay rates of the CPA segments
        segmentlengths_s  : 1D array of the time urations of the CPA segments
        segmentcounts     : 1D array, number of photon events in each CPA segment
        '''
    
       
        filepath_seg = savepath + pre.segment_fit_filename
        if MakeLifetimelst:
    
            crosscorrsA_prep, cctaus_ns_prep = fitexp.PrepCrossCorrs(crosscorrs[0], cctaus_ns, Tau0_bin=-file_metadata['Tau0A_bin'], gap=pre.mygap_bin, rise_time_bin=int(5/pre.dtau_ns))
            crosscorrsB_prep, cctaus_ns_prep = fitexp.PrepCrossCorrs(crosscorrs[1], cctaus_ns, Tau0_bin=-file_metadata['Tau0B_bin'], gap=pre.mygap_bin, rise_time_bin=int(5/pre.dtau_ns))
            crosscorrsAB_prep = crosscorrsA_prep + crosscorrsB_prep
            
            #crosscorrsAB_prep=crosscorrs[0]+crosscorrs[1]
            #cctaus_ns_prep=100-cctaus_ns
            
            fitdata = fitexp.Make_seg_fit_single(crosscorrsAB_prep, cctaus_ns_prep)
            
            load.SaveFitData(filepath_seg, pre.ccmeta, fitdata)
    
        # if not MakeLifetimelst:
            # load the data
            # obsolete here, but left in as a reference of how to load the data
            # ccmeta, fitdata = load.LoadFitData(filepath_seg)
            
