

# =============================================================================
# Figure layout
# =============================================================================
exec(open("../figure_layout.py").read())


# =============================================================================
# file locations
# =============================================================================
simulated_parquetpath = '../Data_simulated_parquet/'
dot_foldername = 'Dot_01/' 

seg_filename = 'nominalsegmentation.csv'

timetags_filenames = ['timestamps_chA_bin.parquet', 'timestamps_chB_bin.parquet', 'timestamps_chR_bin.parquet']
timetags_headers  = ['chA', 'chB', 'chR']

# =============================================================================
# Initialization parameters - measurement time scales
# =============================================================================
minibinsize_ns = 100  #  time between laser pulses [inverse of repetition rate, ns]. Size of minitime bin. This is equivalent to the laser repetition rate
NN = int(6E8)           # total number of minitimebins of size minibinsize (100 ns in this example)
# measurement time in ns = NN* minibinsize_ns
dtau_ns = 0.16461 # timing resolution of the counter card  (taken in this example identical to Becker and Hickl DPC230 card)



# =============================================================================
# Emitter params
# =============================================================================

nlevels = 4 # number of intensity levels
alpha_lst = [1.5, 1.5, 1.5, 1.5] # powerlaws of the levels, dimensionless
I_lst = [5e3, 2E4, 3E4, 5E4] # intensities of the levels in cts/s
g_lst = [0.4,0.25,0.15, 0.05] # decay rates of the levels, in ns^-1
I_noise = 1e2 # background noise in cts/s


# nlevels = 2
# alpha_lst = [1.5, 1.5] # dimensionless
# I_lst = [2E4, 8E4] # in cts/s
# g_lst = [5, 0.1] # in ns^-1
# I_noise = 1e2 # in cts/s



# Jumptime statistics cut off 
Tshortest_ns = 10E6 # shortest switching time in ns
Tlongest_ns =  5E9 # longest switching time, in ns

