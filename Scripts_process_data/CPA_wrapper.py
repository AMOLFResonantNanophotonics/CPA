'''
Wrapper routine that reads in timestamps from parquet data files and 
next runs segementation by CPA. The actual CPA is implemented in findjumps.py

This wrapper executes
    
    -  segmentation by changepoint analysis, executed in findjumps
    
It provides as output the segment list with switching times and counts per segment
    
'''

def CPA_wrapper(MakeChangepoints, Simulated_insteadof_MeasuredData, CPA_insteadof_binned):
    
    #import matplotlib.pyplot as plt
    import numpy as np
    import os
    #import pyarrow.parquet as pq
    #import pandas as pd
    
    import loaddata_functions as load
    import CPA_functions as findjumps
    import fitexp_MLE as fitexp
    # import correlate_jit
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
        # make a metadata file with basic summary info
        # =============================================================================
        file_metadata = {'cts_chA':len(timestamps_chA_bin), 'cts_chB':len(timestamps_chB_bin), 'cts_chR':len(timestamps_chR_bin),
                          'totalPhotonCount':len(timestamps_chA_bin)+len(timestamps_chB_bin),
                          'maxEventTime_bin':timestamps_bin[-1]}
        file_metadata['maxEventTime_ns'] = file_metadata['maxEventTime_bin']*pre.dtau_ns
        file_metadata['dtau_ns'] = pre.dtau_ns
        
        cctaus_bin =  np.arange(0,int(pre.minitimebin_ns/pre.dtau_ns)).astype(np.int64)
        
        # for experimental data in particular,  the zero time, i.e., laser pulse arrival,  might not be in an extremal channel.
        # run routine to find the zero times using (a fraction) of the data
        
        Tau0A_bin, Tau0B_bin= fitexp.Quickdetectorzerotime(timestamps_chA_bin,timestamps_chB_bin, timestamps_chR_bin,int(pre.minitimebin_ns/pre.dtau_ns))
        file_metadata['Tau0A_bin'] = Tau0A_bin
        file_metadata['Tau0B_bin'] = Tau0B_bin
        
        # save the metadata
        print('Saving basic summary data extracted from data in metadatafile:')
        print(savepath+pre.metadata_filename)
        load.WriteDict(savepath+pre.metadata_filename, file_metadata)
    
     
    
        # =============================================================================
        # Segment bounds (i.e., find changepoints. Or, just for plot purposes, fixed width bins)
        # =============================================================================
    
        '''
        output
        seg_idx         : array of the indices of the photons at the found changepoints
        seg_idx_extra   : array of the indices of the photons at the found changepoints, Includes 0, and the total photon count
        seg_times_bin       : array of the time bin of the found changepoint. Unit is counter card resolution
        seg_times_bin_extra : array of the time bin of the found changepoint. Unit is counter card resolution. Includes 0 and the end of the measurement
        '''
    
        filepath_segbounds = savepath + pre.segments_filename
        if CPA_insteadof_binned:
            
            if MakeChangepoints:
                (seg_idx, seg_times_bin) = findjumps.FindMultipleJumps(timestamps_bin, pre.skepticism)
        
                seg_idx_extra   = np.append(np.insert(seg_idx,0,0), len(timestamps_bin))
                seg_times_bin_extra = np.append(np.insert(seg_times_bin,0,0), timestamps_bin[-1])
                
                # save the data
                load.SaveSegments(filepath_segbounds,seg_idx_extra,seg_times_bin_extra)
    
            # else:
                # load the data
                # obsolete here, but left in as a reference of how to load the data
                # seg_idx_extra, seg_times_bin_extra = np.loadtxt(filepath_segbounds).transpose().astype(np.int64)
        
        if not CPA_insteadof_binned:
            seg_times_ns_extra = np.arange(0, file_metadata['maxEventTime_ns'] + pre.analysis_binwidth_ns, pre.analysis_binwidth_ns).astype(np.int64)
            seg_times_bin_extra = (seg_times_ns_extra/pre.dtau_ns).astype(np.int64)
            seg_idx_extra = np.searchsorted(timestamps_bin, seg_times_bin_extra)
            
            load.SaveSegments(filepath_segbounds,seg_idx_extra,seg_times_bin_extra)

