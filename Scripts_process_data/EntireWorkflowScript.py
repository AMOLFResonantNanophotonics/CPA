import CPA_wrapper
import segmented_crosscorrs
import grouping_wrapper
import FDID_wrapper
import plot_timetrace
import switchingtimehistogram
import memoryeffects
import unsegmented_crosscorrs
import g2_tau
import Houel_autocorr

from importlib import reload as reload

reload(CPA_wrapper)
reload(segmented_crosscorrs)
reload(grouping_wrapper)
reload(FDID_wrapper)
reload(plot_timetrace)
reload(switchingtimehistogram)
reload(memoryeffects)
reload(unsegmented_crosscorrs)
reload(g2_tau)
reload(Houel_autocorr)
  


# =============================================================================
# Simulated or measured data
# =============================================================================
Simulated_insteadof_MeasuredData = False

# =============================================================================
# CPA or binned analysis
# =============================================================================
CPA_insteadof_binned = True

# =============================================================================
# Do you want to see all the plots?
# =============================================================================
Show_intermediateplots = False

###############################################################
#   
#   Run scipts for CPA or binning
#   Make segmentation in flowchart
#
###############################################################

MakeChangepoints = True

CPA_wrapper.CPA_wrapper(MakeChangepoints, Simulated_insteadof_MeasuredData, CPA_insteadof_binned)


###############################################################
#   
#   Run scipts for segmented correlations
#   Segmented correlations in flowchart
#
###############################################################

MakeCrossCorrs = True
MakeLifetimelst  = True
segmented_crosscorrs.Segmented_crosscorrs(MakeCrossCorrs, MakeLifetimelst, Simulated_insteadof_MeasuredData, CPA_insteadof_binned)

###############################################################
#   
#   Run scipts for cluster analysis and grouping
#   Grouping in flowchart
#
###############################################################
MakeGrouping = True
PlotBIC = True
PlotOccupation = True

grouping_wrapper.grouping_wrapper(MakeGrouping, PlotBIC, PlotOccupation, Show_intermediateplots, CPA_insteadof_binned, Simulated_insteadof_MeasuredData)


###############################################################
#   
#   Run scipts that plot the results of CPA, using segmented data
#   In flowchart, this is called Segmentation related functions
#
###############################################################



# Plot results from segmentation - time trace, blinking
PlotCombinedTT = True

plot_timetrace.plot_timetrace(PlotCombinedTT, Show_intermediateplots, CPA_insteadof_binned, Simulated_insteadof_MeasuredData)

# Plot results from segmentation - switchinghistogram
# This script won't work if you are doing the binned analysis
MakeSwitchingHist = True
PlotSwitchingHist = True
Show_intermediateplots = True
if CPA_insteadof_binned:
    switchingtimehistogram.switchingtimehistogram(MakeSwitchingHist, PlotSwitchingHist, Show_intermediateplots, CPA_insteadof_binned, Simulated_insteadof_MeasuredData)


# Plot results from segmentation - FDID
MakeFDIDs = True
PlotFDIDs = True
Show_intermediateplots = True
FDID_wrapper.FDID_wrapper(MakeFDIDs, PlotFDIDs, Show_intermediateplots, CPA_insteadof_binned, Simulated_insteadof_MeasuredData)



# Plot results from segmentation - screening for memory effects
memoryeffects.memoryeffects(Show_intermediateplots, CPA_insteadof_binned, Simulated_insteadof_MeasuredData)



###############################################################
#   
#   Run scipts that pertain to the full, unsegmented data
#   In flowchart, this is called Timestamp related functions
#
###############################################################

## fluorescence decay histogram of entire data set

MakeTotalDecayCurve = True
PlotTotalDecayCurve = True

unsegmented_crosscorrs.unsegmented_crosscorrs(MakeTotalDecayCurve, PlotTotalDecayCurve, Show_intermediateplots, CPA_insteadof_binned, Simulated_insteadof_MeasuredData)



## g2 of entire data set
MakeG2 = True
PlotG2 = True

g2_tau.g2_tau(MakeG2, PlotG2, Show_intermediateplots, CPA_insteadof_binned, Simulated_insteadof_MeasuredData)


## Houel et al.   ACS Nano 2015, 9, 1, 886–893  autocorrelation analysis 
MakeLongCorrs = True
PlotLongCorrs = True
    
Houel_autocorr.Houelautocorrelationanalysis(MakeLongCorrs, PlotLongCorrs, Show_intermediateplots, CPA_insteadof_binned, Simulated_insteadof_MeasuredData)
