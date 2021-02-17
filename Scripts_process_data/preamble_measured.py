'''
Definition file specifying the whereabouts and metadata
of:
    
- Example measured data [CsPbBr quantum dot, Ref. ....]
      Parquet format,  int64 timetags for two detector channels and laser pulses
      
- Processed results made by the processing scripts, and used by subsequent scripts.
    - CPA, leads to segmentation into intervals delineated by switching times, with counts, and fitted rates
    
    - g(2), on basis of timetags
    - Totaldecaytrace, on basis of time tags
    - Houel et al.style ACS Nano 2015, 9, 1, 886–893 long time intensity autocorrelation
    
     Postprocessing of segments, i.e.
    - Grouping, on basis of segmented data
    - FDID diagrams
    - Switching time histograms
    - Trace plots
'''

# =============================================================================
# Figure layout
# =============================================================================

import matplotlib as mpl
import acc_functions as acc
exec(open("../figure_layout.py").read())
mycmap = mpl.colors.ListedColormap(acc.Isacolmap(), name='myColorMap', N=acc.Isacolmap().shape[0])  #color map for color plots except FDID
mycmapFDID = mpl.colors.ListedColormap(acc.IsacolmapFDID(), name='myColorMap', N=acc.IsacolmapFDID().shape[0]) # color map for FDID


# =============================================================================
# Specifying where the parquet data lives
# =============================================================================

timetags_filepath = '../Data_measurements_parquet/'
timetags_filenames = ['timestamps_chA_bin.parquet', 'timestamps_chB_bin.parquet', 'timestamps_chR_bin.parquet']
timetags_headers  = ['chA', 'chB', 'chR']


# =============================================================================
# specify a dot you want to study
# =============================================================================
#sig = '_'
sig = '_08'
# If you do not want to set a filter, set sig='_'


# =============================================================================
# Specify where the data should be saved
# =============================================================================
# outputfolder_1 will be combined with either of outputfolder_2 to form a file name,
# depending on whether the data will be binned or segmented with CPA
outputfolder_1 = '../Output_measurements'
outputfolder_2_CPA = '_CPA/'
outputfolder_2_binned = '_binned/'

metadata_filename = 'metadata_dataset.txt'
segments_filename = 'segmentbounds.txt'
crosscorrs_segmented_filename = 'corr_seg.txt'
segment_fit_filename = 'segment_fit_data.csv'

totaldecay_filename = 'corr_unseg.txt'
g2_filename = 'g2.txt'
g2long_filename = 'Houel_autocorr.txt'
grouping_filename = 'grouping.txt'
switchinghist_filename = 'switchinghist.txt'
FDID_filename = 'FDID.txt'
FDID_master_filename = 'FDID_master.txt'


# =============================================================================
# Properties of the measurement/measurement card
# =============================================================================

dtau_ns = 0.16461 # timing resolution of the counter card
minitimebin_ns = 100 # the repetition rate of the pump laser


# =============================================================================
# Electroninc ringing range
# =============================================================================

mygap_bin =  (162,200)  # range in which our data has electronic 
                        # ringing of the DPC230 correlator card. Cut out.

# =============================================================================
# set the skepticism for the CPA algorithm
# =============================================================================
skepticism = 8 # decision threshold for accepting/rejection changepoints. Can be tuned


# =============================================================================
# If you are binning the data, set an analyis binwidth
# =============================================================================
analysis_binwidth_ns = 1e7


# =============================================================================
# Cross-correlation specs, time range
# =============================================================================
cc_tauStart_ns = 0
cc_tauStop_ns = 100
cc_g2range_ns = 210


# metadata for the cross-correlations
ccmeta = {'tauStart_ns':cc_tauStart_ns, 'tauEnd_ns':cc_tauStop_ns,
            'tauStart_bin':int(cc_tauStart_ns/dtau_ns), 'tauEnd_bin':int(cc_tauStop_ns/dtau_ns)
            }

# metadata for the autocorrelations
g2meta = {'tauStart_ns':-cc_g2range_ns, 'tauEnd_ns':cc_g2range_ns,
            'tauStart_bin':-int(cc_g2range_ns/dtau_ns), 'tauEnd_bin':int(cc_g2range_ns/dtau_ns)
            }


# =============================================================================
# memory effects
# =============================================================================
maxdecayrate = 0.5 # in ns^-1

    

