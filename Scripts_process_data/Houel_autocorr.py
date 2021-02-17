'''
Implementation of long-time intensity autocorrelation analysis according to
Houel et al. ACS Nano 2015, 9, 1, 886–893


Fitting Eq. 3 therein to long-time-scale (> milliseconds)  autocorrelation
which for simple two-level dots gives a measure related to the power law exponent of switching

Autocorrelations are obtained using Wahl algorithm with logarithmic coarsening
'''


def Houelautocorrelationanalysis(MakeLongCorrs, PlotLongCorrs, Show_intermediateplots, CPA_insteadof_binned, Simulated_insteadof_MeasuredData):
    
    
    import numpy as np
    import os
    from scipy.optimize import minimize
    import matplotlib.pyplot as plt
    
    import acc_functions as acc
    import loaddata_functions as load
    import correlate_jit as corr
    
    
    #  fit range for the Houel et al analysis
    shortcutoff = 1e-3 # seconds
    longcutoff  = 1 # seconds
    coarsening  = 5 # to calculate long time g(2) you exponentially coarsen the time scale
                    # every n points you double the time step. Good values are from 2 to 10. 
                    # Bigger is slower
                    
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
    print('\n\nRunning routine to perform autocorrelation analysis [Houel et al. ]')
    
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
        timestamps_bin = np.sort(np.concatenate((timestamps_chA_bin, timestamps_chB_bin)))
    
                
        filepath_g2long = savepath + pre.g2long_filename
    
        if MakeLongCorrs:
            
        # =============================================================================
        # do the long time scale correlations
        # =============================================================================
      
            print("\nRun logarithmic coarsening Wahl g(2) autocorrelator") 
            tau, g2 = corr.NormalizedAutoCorrLog(timestamps_bin, coarsening)   #normalized autocorrelation Wahl/Enderlein algorithm
            
        # =============================================================================
        # convert time axis to seconds, clip to cut out
        # =============================================================================
      
        
            tau=tau*pre.dtau_ns*1e-9
            idx=(tau>=shortcutoff)*(tau<longcutoff)
            longcorr_xdata=tau[idx]
            longcorr_ydata=g2[idx]-1
    
        # =============================================================================
        # Run a fit algorithm 
        # =============================================================================
            
            print("\nFit Houel et al. ACS Nano 2015, 9, 1, 886–893 Eq 3") 
            longcorr_x0    = np.array([0.1, 0.003, 0.04]) # Initial guess for the fit
            my_args = (longcorr_xdata, longcorr_ydata, acc.PowerLaw_exp, longcorr_x0)
            res_Longcorr = minimize(acc.lstsq_PowerLaw_exp, longcorr_x0, args=my_args)
            longcorr_fitresults = res_Longcorr['x']
            
            
       # save the data
            print("\nSave autocorrelation and fit to file:") 
            print(filepath_g2long)
            with open(filepath_g2long, 'w') as f:
                print('binned autocorrelation in the style of Houel et al, ACS Nano 9 (2015)\n', file=f)
                print('params = [A, B, C]\nf(t) = (A*(x)**-C) * np.exp(-B*x)', file=f)
                print(str(longcorr_fitresults) + '\n', file=f)
                print('time(s) \t\t\t\t occurence(AU)', file=f)
                np.savetxt(f, np.array([longcorr_xdata,longcorr_ydata]).transpose())
    
        if not MakeLongCorrs:
            print("\n Attempting to load autocorrelation and fit from file:") 
            print(filepath_g2long)
            with open(filepath_g2long, 'r') as f:
                a = f.read().split('\n')
                longcorr_fitresults = np.array(a[4][1:-1].split()).astype(float)
            [longcorr_xdata,longcorr_ydata] = np.loadtxt(filepath_g2long, skiprows=7).transpose()
    
        if PlotLongCorrs: 
    
            
            fig_longcorr, ax_longcorr = plt.subplots(figsize=(4,4))
            yfit = np.array(acc.PowerLaw_exp(longcorr_fitresults, longcorr_xdata))
            ax_longcorr.scatter(longcorr_xdata, longcorr_ydata, label='data')
            ax_longcorr.plot(longcorr_xdata, yfit, 'C1', label='fit')
            
            ax_longcorr.set_xscale('log')
            ax_longcorr.set_yscale('log')
            ax_longcorr.set_xlabel('time (s)')
            ax_longcorr.set_ylabel('occurence (AU)')
            ax_longcorr.legend()
            ax_longcorr.text(0.1,0.15, 'B = %1.2f' %longcorr_fitresults[1], transform=ax_longcorr.transAxes)
            ax_longcorr.text(0.1,0.1, 'C = %1.2f' %longcorr_fitresults[2], transform=ax_longcorr.transAxes)
            ax_longcorr.set_title('Intensity autocorrelation (log-log)')
            fig_longcorr.savefig(savepath+'Houel_autocorr.png')
           
    
            if pre.sig=='_' and not Show_intermediateplots:
                plt.close(fig_longcorr)
    
