'''
Handler for reading and writing files
'''

def Findlinebreaks(loadFilepath):
    breaks = []
    print(loadFilepath)
    with open(loadFilepath, 'r') as f:
        for idx,line in enumerate(f):
            if line == '\n' or len(line)<3:
                breaks += [idx]
    return breaks

def WriteDict(filename, mydict):
    with open(filename,"w") as f:
        for key in mydict.keys():
            print("%s: %s" % (key, mydict[key]), file=f)

def ReadDict(loadfile, mybreak):
    mydict = {}
    with open(loadfile, 'r') as f:
        for idx,line in enumerate(f):
            if idx<mybreak:
                (key, val) = line.split(': ')
                if val.startswith('[') and val.endswith(']\n'):
                    mydict[str(key)] = eval(val)
                    continue
                if '.' in val:
                    mydict[str(key)] = float(val)
                    continue
                else:
                    mydict[str(key)] = int(val)
                    continue
    return mydict

# =============================================================================
# Timestamps
# =============================================================================

def LoadTimeStamps(datafolder, filenames, headers):
    import numpy as np
    import pyarrow.parquet as pq
    
    print('Reading parquet files')
    print('Folder: '+ datafolder)
    print('Files: ')
    print(filenames)
    
    timestamps_chA_bin = np.array(pq.read_pandas((datafolder + filenames[0])).to_pandas()[headers[0]])
    timestamps_chB_bin = np.array(pq.read_pandas((datafolder + filenames[1])).to_pandas()[headers[1]])# + int(channelshift_bin)
    timestamps_chR_bin = np.array(pq.read_pandas((datafolder + filenames[2])).to_pandas()[headers[2]])
    
    return timestamps_chA_bin, timestamps_chB_bin, timestamps_chR_bin




# =============================================================================
# CPA segments and counts
# =============================================================================

def SaveSegments(filepath_segbounds,seg_idx_extra,seg_times_bin_extra):
    import numpy as np
    print('Saving CPA segment bounds to')
    print(filepath_segbounds)
    np.savetxt(filepath_segbounds, np.array([seg_idx_extra, seg_times_bin_extra]).transpose(), fmt='%d', header='segment photon idx, segment time(bin)')


def LoadSegments(filepath_segbounds):
    import numpy as np
    
    
    print('Reading CPA segmentation from:')
    print(filepath_segbounds)
    seg_idx_extra, seg_times_bin_extra = np.loadtxt(filepath_segbounds).transpose().astype(np.int64)
    
    segmentlengths = np.diff(seg_times_bin_extra)
    segmentcounts  = (np.diff(seg_idx_extra)).astype(float)
    
    return     seg_idx_extra, seg_times_bin_extra, segmentlengths,segmentcounts




# =============================================================================
# Segmented cross-correlations
# =============================================================================

def SaveSegmentedCC(filepath_cc, crosscorrs, cctaus_ns):
    import numpy as np
    print('Leaving decay histograms for each segment in')
    print(filepath_cc)
    
    with open(filepath_cc,"w") as f:    
        fmtfloat = "{0:f}, {1}"
        print("\nrows->tau (relative time in correlation), columns-> CPA segment or bin number", file=f)
    
        print("Lifetime A", file=f)
        for x in range(np.shape(crosscorrs)[2] ):
            print(fmtfloat.format(cctaus_ns[x], ' '.join(str(i) for i in crosscorrs[0,:,x])), file=f)
    
        print("\nLifetime B", file=f)
        for x in range(np.shape(crosscorrs)[2] ):
            print(fmtfloat.format(cctaus_ns[x], ' '.join(str(i) for i in crosscorrs[1,:,x])), file=f)

def LoadSegmentedCC(filepath_cc):
    import numpy as np
    
    print('Trying to open segmented decay traces ')
 
    
    breaks = Findlinebreaks(filepath_cc)
    crosscorrs = np.array([
            np.genfromtxt(filepath_cc, delimiter=' ', dtype=int, skip_header=breaks[0]+3, max_rows=breaks[1]-breaks[0]-3).transpose()[1:],
            np.genfromtxt(filepath_cc, delimiter=' ', dtype=int, skip_header=breaks[1]+2).transpose()[1:]
            ])
    cctaus_ns = np.genfromtxt(filepath_cc, delimiter=', ', skip_header=breaks[0]+3, max_rows=breaks[1]-breaks[0]-3, usecols=0)    
    
    return crosscorrs, cctaus_ns

# =============================================================================
# Fits to the segmented cross-correlations
# =============================================================================

def SaveFitData(filepath_seg, ccmeta, fitdata):
    import numpy as np
    
    print('Leaving segment info incl. fitted rates in')
    print(filepath_seg)
        
    with open(filepath_seg,"w") as f:
        for key in ccmeta.keys():
            print("%s: %s" % (key, ccmeta[key]), file=f)
        print(' ', file=f)
        print((', ').join(fitdata.columns.values), file=f)
        np.savetxt(f, fitdata.values, delimiter=', ')


def LoadFitData(filepath_seg):
    import pandas as pd

    print('Trying to open segment info incl. fitted rates [output of an earlier CPA wrapper run]: ')
    
    
    breaks = Findlinebreaks(filepath_seg)
    ccmeta = ReadDict(filepath_seg, breaks[0])

    fitdata = pd.read_csv(filepath_seg, skiprows=breaks[0]+1, header=0, sep='\s*,\s*', engine='python')
    
    return ccmeta, fitdata

# =============================================================================
# Unsegmented Correlations
# =============================================================================

def SaveUnsegmentedCorr(filepath_totdecay, fit, cctaus_ns, crosscorrsAB_flat):
    import numpy as np
    
    
    print('Saving total time trace decay histogram, and fit function: ')
    print(filepath_totdecay)
    
    
    # save the correlation and the fit
    with open(filepath_totdecay, 'w') as f:
        for i in fit:
            if i=='func':
                f.write(str(i) + ': ' + fit[i])
            else:
                f.write(str(i) + ': ' + str([j for j in fit[i]])+ '\n')
        f.write('\n')
        f.write('cctaus (ns), correlation\n')
        filepath_totdecay
        np.savetxt(f,np.vstack((cctaus_ns,crosscorrsAB_flat)).T)


def LoadUnsegmentedCorr(filepath_totdecay):
    import numpy as np

     
    print('Trying to open total time trace decay histogram, and fit function [output of an earlier CPA wrapper run]: ')
    
    breaks = Findlinebreaks(filepath_totdecay)
    
    fit = ReadDict(filepath_totdecay, breaks[0])

    [cctaus_ns, crosscorrsAB_flat] = np.loadtxt(filepath_totdecay, skiprows=breaks[0]+2).transpose()
    
    return fit, cctaus_ns, crosscorrsAB_flat











