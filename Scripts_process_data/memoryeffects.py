def memoryeffects(Show_intermediateplots, CPA_insteadof_binned, Simulated_insteadof_MeasuredData):
    
    import numpy as np
    import os 
    import matplotlib.pyplot as plt
    from matplotlib import gridspec as gsp
    
    import matplotlib as mpl
    
    
    import loaddata_functions as load
    import acc_functions as acc
    
    

    # =============================================================================
    # specify range of intensities you want to see w.r.t. the mean
    # =============================================================================
    aa = 2.5
    
    
    
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
    print('\n\nRunning routine to plot overview of memory effects in observables from CPA Segmentation')
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
        segmentcounts       : array of the number of counts in every segment
        '''
        
        filepath_segbounds = savepath + pre.segments_filename
        seg_idx_extra, seg_times_bin_extra, segmentlengths,segmentcounts=load.LoadSegments(filepath_segbounds)
        
           
        segmentlengths_s=segmentlengths*pre.dtau_ns*1e-9       # in seconds  [durations]
        seg_times_s = seg_times_bin_extra[:-1]*pre.dtau_ns*1e-9 # in seconds   [cumulative time]
        
        # =============================================================================
        # Read in previously fitted lifetimes and intensity lists
        # =============================================================================
        '''
        fitdata contains 4 parameters of the fitted single-exponential decay
        fit_A         : the scaling factor
        fit_gamma     : the decay rate, in ns^-1
        fit_gamma_err : the error on the fitted decay rate
        fit_bg        : the background noise
        '''
    
        filepath_fitdata = savepath + pre.segment_fit_filename
        ccmeta,fitdata=load.LoadFitData(filepath_fitdata)
    
        # =============================================================================
        # Make Figures
        # =============================================================================
        totalcount=sum(segmentcounts)
        totalduration=max(seg_times_s)
        meanI   = totalcount/totalduration
        scale_I = int(np.log10(aa*meanI))
        meanI   = meanI/(10**scale_I)
        meanlogT = np.mean(np.log10(segmentlengths_s))
        
        meangamma       = np.mean(fitdata['fit_gamma'])
        
        
        
        seg_times_s = seg_times_bin_extra[:-1]*pre.dtau_ns*1e-9
        intensity = segmentcounts/np.diff(seg_times_bin_extra)/pre.dtau_ns*1e9/(10**scale_I)
        gamma = fitdata['fit_gamma']
        
        fig = plt.figure(figsize=(18/2.51,6/2.51))
        gs0 = gsp.GridSpec(1, 3, width_ratios=[1.3,1,1], wspace=0.4)
        gsL = gsp.GridSpecFromSubplotSpec(3, 1, gs0[0], height_ratios=[1,1,1], hspace=0.3) 
        
        axL1 = fig.add_subplot(gsL[0])
        axL2 = fig.add_subplot(gsL[1], sharex=axL1)
        axL3 = fig.add_subplot(gsL[2], sharex=axL1)
        
        ax2 = fig.add_subplot(gs0[1])
        ax3 = fig.add_subplot(gs0[2])
        
        out = np.histogram2d(seg_times_s, intensity, bins=100, range=[[0,totalduration], [0,3*meanI]])
        axL1.imshow(np.transpose(out[0]), extent=[out[1][0], out[1][-1], out[2][0], out[2][-1]], aspect='auto', origin='lower',cmap=pre.mycmap)
        axL1.set_ylabel("I ($10^{%1d}$/s)" %scale_I)
        plt.setp(axL1.get_xticklabels(), visible=False)
        
        out = np.histogram2d(seg_times_s, gamma, bins=100, range=[[0,seg_times_s[-1]], [0,pre.maxdecayrate]])
        axL2.imshow(np.transpose(out[0]), extent=[out[1][0], out[1][-1], out[2][0], out[2][-1]], aspect='auto', origin='lower',cmap=pre.mycmap)
        axL2.set_ylabel('$\gamma$ (ns$^{-1}$)')
        plt.setp(axL2.get_xticklabels(), visible=False)
        axL2.set_yticks([0, pre.maxdecayrate/2-0.05, pre.maxdecayrate-0.1])
        axL2.set_ylim([0,pre.maxdecayrate])
        
        out = np.histogram2d(seg_times_s, np.log10(segmentlengths_s), range=[[0,totalduration], [-4.5,0.5]], bins=100)
        axL3.imshow(np.transpose(out[0]), extent=[out[1][0], out[1][-1], out[2][0], out[2][-1]], aspect='auto', origin='lower',cmap=pre.mycmap)
        axL3.set_xlabel('time (s)')
        axL3.set_ylabel('$\log_{10}(T_n)$')
        axL3.set_ylim([-4.5,0.5])
        axL3.set_yticks([-4,-2,0])
        
        out = np.histogram2d(np.log10(segmentlengths_s), intensity, bins=100, range=[[-5,1],[0,aa*meanI]])
        ax2.imshow(np.transpose(out[0]), extent=[out[1][0], out[1][-1], out[2][0], out[2][-1]], aspect='auto', origin='lower',cmap=pre.mycmap)
        ax2.set_xlabel('$\log_{10}(T_n)$')
        ax2.set_ylabel("I ($10^{%1d}$ cts/s)" %scale_I)
        
        out = np.histogram2d(np.log10(segmentlengths_s), gamma, bins=100, range=[[meanlogT-3,meanlogT+3],[0,pre.maxdecayrate]])
        ax3.imshow(np.transpose(out[0]), extent=[out[1][0], out[1][-1], out[2][0], out[2][-1]], aspect='auto', origin='lower',cmap=pre.mycmap)
        ax3.set_xlabel('$\log_{10}(T_n)$')
        ax3.set_ylabel('$\gamma$ (ns$^{-1}$)')
        
        fig.align_ylabels([axL1, axL2, axL3])
        
        fig.savefig(savepath+'memoryeffects_aging_Tn.pdf',bbox_inches='tight')
        
        ################
        fig = plt.figure(figsize=(18/2.51, 15.5/2.51))
        gs0 = gsp.GridSpec(2, 1, height_ratios=[2,0.7], hspace=0.3)
        gs1 = gsp.GridSpecFromSubplotSpec(2, 3, gs0[0], width_ratios=[1,1,1], wspace=0.4, height_ratios=[1,1], hspace=0.1) 
        gs2 = gsp.GridSpecFromSubplotSpec(1, 3, gs0[1], width_ratios=[1,1,1], wspace=0.4) 
        ax00 = fig.add_subplot(gs1[0,0])
        ax01 = fig.add_subplot(gs1[0,1])
        ax02 = fig.add_subplot(gs1[0,2])
        ax10 = fig.add_subplot(gs1[1,0])
        ax11 = fig.add_subplot(gs1[1,1])
        ax12 = fig.add_subplot(gs1[1,2])
        ax20 = fig.add_subplot(gs2[0])
        ax21 = fig.add_subplot(gs2[1])
        ax22 = fig.add_subplot(gs2[2])
        
        shift = 1
        out = np.histogram2d(intensity[:len(intensity)-shift], intensity[shift:], bins=50, range=[[0,2.5*meanI],[0,2.5*meanI]])
        out2 = out[0] /(np.sum(out[0],axis=0)[:,None]+1E-10) ## the sum usually should be integer. If it is 0, add infinitesimal offset, so that out2 is 0
        ax00.imshow(np.transpose(out2), extent=[out[1][0], out[1][-1], out[2][0], out[2][-1]], origin='lower',cmap=pre.mycmap)
        # ax00.set_xlabel('$I_n$ ($10^{%1d}$ cts/s)' %scale_I)
        ax00.set_ylabel('$I_{n+1}$ ($10^{%1d}$ cts/s)' %scale_I)
        ax00.set_xticklabels([])
        
        out = np.histogram2d(gamma[:len(intensity)-shift], gamma[shift:], bins=50, range=[[0,pre.maxdecayrate],[0,pre.maxdecayrate]])
        out2 = out[0] /(np.sum(out[0],axis=0)[:,None]+1E-10) ## the sum usually should be integer. If it is 0, add infinitesimal offset, so that out2 is 0
        ax01.imshow(np.transpose(out2), extent=[out[1][0], out[1][-1], out[2][0], out[2][-1]], origin='lower', vmax=0.1,cmap=pre.mycmap)
        # ax01.set_xlabel('$\gamma_n$ (ns$^{-1}$)')
        ax01.set_ylabel('$\gamma_{n+1}$ (ns$^{-1}$)')
        ax01.set_xticklabels([])
        ax01.set_yticks([0,0.1,0.2,0.3,0.4])
        ax01.set_yticklabels([0,0.1,0.2,0.3,0.4])
        
        out = np.histogram2d(np.log10(segmentlengths_s)[:len(intensity)-shift], np.log10(segmentlengths_s)[shift:], bins=50, range=[[-4,-0.75],[-4,-0.75]])
        norm = np.sum(out[0],axis=0)
        norm[np.where(norm==0)] = 1
        out2 = out[0] / norm[:,None]
        ax02.imshow(np.transpose(out2), extent=[out[1][0], out[1][-1], out[2][0], out[2][-1]], origin='lower', vmax=0.2,cmap=pre.mycmap)
        # ax02.set_xlabel('$\log_{10}(T_{n} \,(\mathrm{s}))$')
        ax02.set_ylabel('$\log_{10}(T_{n+1} \,(\mathrm{s}))$')
        ax02.set_xticklabels([])
        ax02.set_yticks([-4, -3, -2, -1])
        ax02.set_yticklabels([-4, -3, -2, -1])
        
        shift = 2
        out = np.histogram2d(intensity[:len(intensity)-shift], intensity[shift:], bins=50, range=[[0,2.5*meanI],[0,2.5*meanI]])
        out2 = out[0] /(np.sum(out[0],axis=0)[:,None]+1E-10) ## the sum usually should be integer. If it is 0, add infinitesimal offset, so that out2 is 0
        ax10.imshow(np.transpose(out2), extent=[out[1][0], out[1][-1], out[2][0], out[2][-1]], origin='lower',cmap=pre.mycmap)
        ax10.set_xlabel('$I_n$ ($10^{%1d}$ cts/s)' %scale_I)
        ax10.set_ylabel('$I_{n+2}$ ($10^{%1d}$ cts/s)' %scale_I)
        
        out = np.histogram2d(gamma[:len(intensity)-shift], gamma[shift:], bins=50, range=[[0,pre.maxdecayrate],[0,pre.maxdecayrate]])
        out2 = out[0] /(np.sum(out[0],axis=0)[:,None]+1E-10) ## the sum usually should be integer. If it is 0, add infinitesimal offset, so that out2 is 0
        ax11.imshow(np.transpose(out2), extent=[out[1][0], out[1][-1], out[2][0], out[2][-1]], origin='lower', vmax=0.1,cmap=pre.mycmap)
        ax11.set_xlabel('$\gamma_n$ (ns$^{-1}$)')
        ax11.set_ylabel('$\gamma_{n+2}$ (ns$^{-1}$)')
        ax11.set_yticks([0,0.1,0.2,0.3,0.4])
        ax11.set_yticklabels([0,0.1,0.2,0.3,0.4])
        
        out = np.histogram2d(np.log10(segmentlengths_s)[:len(intensity)-shift], np.log10(segmentlengths_s)[shift:], bins=50, range=[[-4,-0.75],[-4,-0.75]])
        norm = np.sum(out[0],axis=0)
        norm[np.where(norm==0)] = 1
        out2 = out[0] / norm[:,None]
        ax12.imshow(np.transpose(out2), extent=[out[1][0], out[1][-1], out[2][0], out[2][-1]], origin='lower', vmax=0.2,cmap=pre.mycmap)
        ax12.set_xlabel('$\log_{10}(T_{n} \,(\mathrm{s}))$')
        ax12.set_ylabel('$\log_{10}(T_{n+2} \,(\mathrm{s}))$')
        ax12.set_yticks([-4, -3, -2, -1])
        ax12.set_yticklabels([-4, -3, -2, -1])
        
        out = np.correlate(intensity,intensity, 'full')/(len(intensity)*meanI**2)
        ax20.plot(np.arange(-len(intensity)+1,len(intensity)), out)
        ax20.set_xlim((-20, 20))
        # ax20.set_ylim((out[], 1.3))
        ax20.set_xlabel('relative segment nr')
        ax20.set_ylabel('autocorr($I$)')
        
        ax21.plot(np.arange(-len(intensity)+1,len(intensity)), np.correlate(gamma,gamma, 'full')/(len(gamma)*meangamma**2))
        ax21.set_xlim((-20, 20))
        ax21.set_xlabel('relative segment nr')
        ax21.set_ylabel('autocorr($\gamma$)')
        
        ax22.plot(np.arange(-len(segmentlengths_s)+1,len(segmentlengths_s)), np.correlate(segmentlengths_s, segmentlengths_s, 'full')/(len(intensity)*np.mean(segmentlengths_s)**2))
        ax22.set_xlim((-20, 20))
        ax22.set_xlabel('relative segment nr')
        ax22.set_ylabel('autocorr(T)')
        
        fig.savefig(savepath+'memoryeffects_correlations.pdf',bbox_inches='tight')
        
        if pre.sig=='_' and not Show_intermediateplots:
            plt.close('all')
