A toolbox that provides functionality to perform changepoint analysis (CPA) on data containing time-tagged events. 
Initially developed for time-correlated single-photon counting (TCSPC) measurements, but the principle holds generally for all types of events. 
The full description of the method will be uploaded to arXiv shortly. 

An arXiv submission containing benchmarks and an explanation of the method can be found here: https://arxiv.org/abs/2102.11044

An arXiv submission showing the method applied to TCSPC data from lead halide perovskites can be found here: http://arxiv.org/abs/2102.09333


Changelog 20210505: 
(1)added comments to preamble files

(2) made a new routine that calculates FLIDs as well as FDIDs. 
It is purely optional, and is inside the existing FDID_wrapper script. 
Mostly identical to the routine that makes FDIDs

(3) Fixed a mistake in the file preamble_for_simulation
NN set to 1e9 and the comment corrected

(4) fixed up grouping_wrapper.py to prevent errors when there are too few changepoints

(5) made mingroups and maxgroups (for the grouping algorithm) variables that can be set in the preamble
