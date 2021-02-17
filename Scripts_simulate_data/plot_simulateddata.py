
def PlotSimulatedTrace():
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import gridspec as gsp
    
    import preamble_for_simulation as pre
    import loaddata_functions as load
    
    
    binsize_ns    = 1e6     # binning size, in s. For visualization purposes ONLY
    Nbins = int(pre.NN*pre.minibinsize_ns/binsize_ns)
    
    
    filepath_timestamps   = pre.simulated_parquetpath + pre.dot_foldername
    filepath_segmentation = pre.simulated_parquetpath + pre.dot_foldername + pre.seg_filename
    
    # =============================================================================
    # Load the jumptimes
    # =============================================================================
    [jumptimes, trace_I, trace_g] = np.loadtxt(filepath_segmentation, skiprows=1, delimiter=',').transpose()
    jumptimes = np.insert(jumptimes, 0, 0)
    
    
    # =============================================================================
    # load the timestamps
    # =============================================================================

    timestamps_chA_bin, timestamps_chB_bin, timestamps_chR_bin = load.LoadTimeStamps(filepath_timestamps, pre.timetags_filenames, pre.timetags_headers)
    timestamps_bin = np.sort(np.concatenate((timestamps_chA_bin, timestamps_chB_bin)))
    
    
    # =============================================================================
    # Make the plot
    # =============================================================================
    
    photons_x = np.arange(0,pre.NN*pre.minibinsize_ns,binsize_ns)*1E-9 # the photon trace
    photons_y = np.histogram(timestamps_bin,Nbins)[0]/binsize_ns*1E9 # nr of photons/s
    
    Itrace_plot_x = np.repeat(jumptimes*pre.minibinsize_ns,2)[1:-1]*1E-9
    Itrace_plot_y = np.repeat(trace_I,2)
    
    gtrace_plot_x = np.repeat(trace_g,2)
    
    
    fig = plt.figure(figsize=(20/2.51, 6/2.51))
    gs0 = gsp.GridSpec(2, 1, height_ratios=[1,1], hspace=0.05)
    ax_I = fig.add_subplot(gs0[0])
    ax_g = fig.add_subplot(gs0[1], sharex=ax_I)
    ax_I.plot(photons_x, photons_y)
    ax_I.plot(Itrace_plot_x, Itrace_plot_y)

    
    ax_g.plot(Itrace_plot_x, gtrace_plot_x, c='C1')
    
    I_ticks = ax_I.get_yticks()[ax_I.get_yticks() > 0]
    I_tickslabels = [str(int(x/10**int(np.log10(x)))) + '$\cdot 10^{%d}$' % int(np.log10(x)) for x in I_ticks]
    
    I_ticks=np.append([0.0], I_ticks)
    ax_I.set_yticks(I_ticks)
    I_ticks = ax_I.get_yticks()
   
    I_tickslabels = ['0'] + I_tickslabels
    ax_I.set_yticklabels(I_tickslabels)
    ax_I.set_ylabel('I (cts/)')
    ax_g.set_xlabel('Time (s)')
    ax_g.set_ylabel('$\gamma$ (ns$^{-1}$)')

    fig.savefig(pre.simulated_parquetpath + pre.dot_foldername + 'trace.pdf')
