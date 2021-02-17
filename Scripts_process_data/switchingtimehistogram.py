'''
Routine that converts the switching times in CPA-segments 
constructed by CPA_step1,  into histograms of dwell times, fitting a power law

'''
def switchingtimehistogram(MakeSwitchingHist, PlotSwitchingHist, Show_intermediateplots, CPA_insteadof_binned, Simulated_insteadof_MeasuredData):
    
    
    import numpy as np
    import os
    import scipy.optimize
    import matplotlib.pyplot as plt
    
    import loaddata_functions as load
    import acc_functions as acc
    
    
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
    
    # =============================================================================
    # start processing the data
    # =============================================================================
    print('\n\nRunning routine to plot residence time histograms')
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
        # Switching histogram
        # =============================================================================
        '''
        Calculate the histogram of the times between retrieved switching events
    
        input
        retrievedjumptimes_bin_extra : array of the time bin of the found changepoint. Unit is counter card resolution. Includes 0 and the end of the measurement
        pre.dtau_ns                      : timing resolution of the counter card
        rollofftime                  : can be chosen. The shortest time where the historgam can stll be
    
        output
        powerlawfit_SH     : expontnt of the powerlaw fitted tot the switching histogram
        powerlawfit_SH_var : error on the exponent
        SH_xdata           : the x-values of the histogram bins in which the switching times were binned
        SH_ydata           : the histogram of the switching times
        '''
        
        filepath_switching = savepath + pre.switchinghist_filename
        
        if MakeSwitchingHist:
            histedges = 10**(np.linspace(-3,np.log10(np.amax(segmentlengths_s)), num=50))# we're making bins that grow exponentially wider)
    
            Nt  = np.histogram(segmentlengths_s,bins=histedges)[0] # the nr of occurences in the variable bins
            xd = 10**((np.log10(histedges[1::])+np.log10(histedges[0:-1]))/2) # the position of th bins
            dt = np.diff(histedges) # the size of the variable-sized bins
    
    
            # take the last few elements (around where the power-law ends)
            # if there are elements with value 0, take them out. They mess with the plot
            SH_ydata = (Nt/dt)
            select = np.nonzero(SH_ydata)
            SH_ydata = SH_ydata[select]
            SH_ydata = np.log10(SH_ydata)
            SH_xdata = np.log10(xd[select])
            x0    = np.array([np.amax(SH_ydata), 2.5])
            sigma = np.ones(len(SH_ydata))
    
    
            # now fit a power law. Typically the CPA itself causes a short time roll off
            # as a rule of thumb, exclude segments < 10 ms from the fit
            
            rollofftime = 1e-2
            rolloffbin = (np.abs(xd - rollofftime)).argmin()
    
    
            powerlawfit_SH, powerlawfit_SH_var = scipy.optimize.curve_fit(acc.linear, SH_xdata[rolloffbin:], SH_ydata[rolloffbin:], x0, sigma[rolloffbin:])
            powerlaw_height = powerlawfit_SH[0]
            powerlaw_slope = powerlawfit_SH[1]
            powerlaw_error = powerlawfit_SH_var[1,1]
    
            # saving the powerlaw data
            print('Saving the switching time histogram to file')
            print(filepath_switching)
            with open(filepath_switching, 'w') as f:
                print('The powerlaw exponent of the long tail of this dot: '+str(powerlaw_slope), file=f)
                print('The height of the powerlaw: '+str(powerlaw_height), file=f)
                print('\nThe error on the fit: \n'+str(powerlaw_error), file=f)
                print('\nThe bin nr after which the slope is powerlaw: '+str(rolloffbin), file=f)
                print('\nThe fitted data\nlog(xdata), log(ydata)', file=f)
                np.savetxt(f, np.array([SH_xdata, SH_ydata]).transpose())
    
        if not MakeSwitchingHist:
            #loading the powerlaw data
            print('Attempting to load previously saved switching time histogram from file')
            print(filepath_switching)
    
            with open(filepath_switching) as f:
                a = f.read().split('\n')
                powerlaw_slope  = float(a[0].split(': ')[1])
                powerlaw_height = float(a[1].split(': ')[1])
                rolloffbin = int(a[6].split(': ')[1])
            SH_xdata, SH_ydata = np.loadtxt(filepath_switching, skiprows=10).transpose()
    
        if PlotSwitchingHist:
            figSH, axswitching = plt.subplots(figsize=(8/2.51, 8/2.51))
            axswitching.scatter((10**SH_xdata), 10**SH_ydata)
            axswitching.plot((10**SH_xdata[rolloffbin:]), 10**acc.linear(SH_xdata[rolloffbin:], powerlaw_height, powerlaw_slope), color='C1')
            axswitching.set_xscale('log')
            axswitching.set_yscale('log')
            axswitching.set_xlim((5e-4,1.1e1))
            axswitching.set_ylim((1e0,1.5e7))
            axswitching.set_xlabel('$\Delta t$ (s)')
            axswitching.set_ylabel('occurrence')
            axswitching.set_title('switching times')
            # axswitching.plot(10**SH_xdata, 10**acc.linear(SH_xdata, x0[0], x0[1]), color='r')
            axswitching.text(0.1, 0.1, r'$\alpha$ = %1.1f' %powerlaw_slope, transform=axswitching.transAxes, fontsize=8)
            figSH.savefig(savepath+'switchingHist.png') 
    
     
    
