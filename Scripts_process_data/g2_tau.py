'''
 Wrapper routine that creates g(2) from timestamped photon data. 
 
 Reads in parquet time stamped data
 Correlates them using Wahl algorithm for which it calls  correlate_jit
 Plots and saves the result

'''
def g2_tau(MakeG2, PlotG2, Show_intermediateplots, CPA_insteadof_binned, Simulated_insteadof_MeasuredData):
    
    import numpy as np
    import os
    import matplotlib.pyplot as plt
    #from numba import jit,int64
    
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
        outputfolder = pre.outputfolder_1 + pre.outputfolder_2_CPA
    else:
        outputfolder = pre.outputfolder_1 + pre.outputfolder_2_binned
    
    
    # =============================================================================
    # get to work on the specified dots
    # =============================================================================
     
    Dotlist = [i for i in os.listdir(pre.timetags_filepath) if i.startswith('Dot_') and pre.sig in i]
    print('\n\nRunning routine to construct and plot g(2) from entire time trace')
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
        # Load the timestamps
        # =============================================================================
        '''
        timestamps_chX_bin : all events in channel X E (A, B, R)
        timestamps_bin     : all events in channels A and B, chronologically
        '''
    
        timestamps_chA_bin, timestamps_chB_bin, timestamps_chR_bin = load.LoadTimeStamps(pre.timetags_filepath+'Dot_%02d/' %dot_idx, pre.timetags_filenames, pre.timetags_headers)
        
    
        
    
    
        # =============================================================================
        # The g2(tau)
        # =============================================================================
        '''
        input
        timestamps_chX_bin : all events in channel X E (A, B, R)
        g2params           : the parameters for the g2 calculation
        g2meta             : the parameters of the g2 calculation and the metadata
    
        output
        g2          : the caclulated g2
        g2taus_ns   : the time shifts used for calculating the g2, in ns
        '''
        
        filepath_g2 = savepath + pre.g2_filename
        
        if MakeG2:
    
            g2taus_bin =  np.arange(pre.g2meta['tauStart_bin'],pre.g2meta['tauEnd_bin']).astype(np.int64)
            g2taus_ns = g2taus_bin* pre.dtau_ns
            
            print("\nStart Wahl-based correlator for g(2)")
            g2 = correlate_jit.Corr(timestamps_chB_bin, timestamps_chA_bin, g2taus_bin)
            
            print("Save the outcome to: ")
            print(filepath_g2)
            with open(filepath_g2,"w") as f:
                for key in pre.g2meta.keys():
                    print("%s: %s" % (key, pre.g2meta[key]), file=f)
                print("\ntau, g2(tau)", file=f)
                np.savetxt(f, np.array([g2taus_ns,g2], dtype=float).transpose() )
            
            
        if not MakeG2:
            print("Attempting to load a previously calculated g(2) from: ")
            print("filepath_g2")
            breaks = load.Findlinebreaks(filepath_g2)
            g2taus_ns, g2 = np.loadtxt(filepath_g2, skiprows=breaks[0]+2, delimiter=' ').transpose()
    
        if PlotG2:
            figg2, axg2 = plt.subplots(figsize=(6/2.51, 6/2.51))
            axg2.scatter(g2taus_ns, g2)
            axg2.set_xlabel('\u03C4 (ns)')
            axg2.set_ylabel('$g^2$(\u03C4)')
            axg2.set_title('autocorrelation signal')
            figg2.savefig(savepath+'g2.png')
            
            if not Show_intermediateplots and pre.sig == '_':
                plt.close(figg2)
