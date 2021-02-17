'''  
Routines to fit fluorescence lifetime using the Poisson statistics maximum likelihood estimator
Bazjer  Eur Biophys J 20, 247

''' 

def Quickdetectorzerotime(timeA_bin,timeB_bin,timeR_bin,nbins):
        # for experimental data in particular,  the zero time, i.e., laser pulse arrival,  might not be in an extremal channel.
        # also different detectors might have different cable lengths, i.em zero times.
        #
        # a working method to find tau0 is to correlate a small part of the dataset.
        # here, the first "target" photons are used
        # 
        # input:  detectorA, B and reference time tag list as int64
        # number of bins matching the laser period
        
        import numpy as np
        import correlate_jit
        target=int(1E4)
        cutoffA=min((len(timeA_bin)),target)
        cutoffB=min((len(timeB_bin),target))
        cutoffR=min((len(timeR_bin),3*target))
        
        cctaus_bin =  np.arange(0,int(nbins)).astype(np.int64)
        Tau0A_bin = np.argmax(correlate_jit.Corr(timeR_bin[0:cutoffR], timeA_bin[0:cutoffA], cctaus_bin))
        Tau0B_bin = np.argmax(correlate_jit.Corr(timeR_bin[0:cutoffR], timeB_bin[0:cutoffB], cctaus_bin))
        
        return Tau0A_bin,Tau0B_bin

def PrepCrossCorrs(crosscorrs, taus, Tau0_bin=0, gap=(0,0), rise_time_bin=0):
    '''
    This function cleans fluorescence decay traces in the sense that:
        * shifts data to the zero-point in time 
        * flips time axis, since all data is specified as reverse start stop. (Undoes the reverse)
        * clips out detector rise time
        * clips out interval "gap" for data with artefacts (electronic ringing)
        
    For synthetic data this routine does in essence nothing at all.
    '''


    import numpy as np

    # flip around the data and set the correct Tau0
    crosscorrs = np.roll(crosscorrs,Tau0_bin,axis=1)[:,::-1]

    # get rid of the electronic ringing
    if not gap == (0,0):
        taus = np.delete(taus, np.arange(gap[0],gap[1]))
        crosscorrs = np.delete(crosscorrs, np.arange(gap[0],gap[1]), axis=1)


    # get rid of the rise of the peak. Cut a couple elements off the END of the seqence for the rise of the curve
    if not rise_time_bin == 0:
        taus = taus[:-rise_time_bin]
        crosscorrs = crosscorrs[:,:-rise_time_bin]
    return crosscorrs, taus


def Make_seg_fit_single(crosscorrs, taus, intermediateplots=False):
    '''
    This function fits a fluorescence decay trace corresponding to a single data segment
    It is currently hardcoded to fit single exponentials

    '''
    import numpy as np
    import pandas as pd
    import warnings
    # suppressing these warnings. The fit works, even in cases where the warning is thrown
    warnings.filterwarnings("ignore", message="divide by zero encountered in log")
    warnings.filterwarnings("ignore", message="invalid value encountered in multiply")
    warnings.filterwarnings("ignore", message="invalid value encountered in subtract")

    A_lst         = np.zeros(len(crosscorrs))
    gamma_lst     = np.zeros(len(crosscorrs))
    gamma_err_lst = np.zeros(len(crosscorrs))
    bg_lst        = np.zeros(len(crosscorrs))

    #type casting so as not to confuse the fit function
    taus=taus.astype(float)
    crosscorrs=crosscorrs.astype(float)

    for i in range(len(crosscorrs)):
        fitresults = minimize_func(taus, crosscorrs[i], decay='single')
        A_lst[i]         = fitresults['fit_result'][0]
        gamma_lst[i]     = fitresults['fit_result'][1]
        gamma_err_lst[i] = fitresults['se'][1]
        bg_lst[i]        = fitresults['fit_result'][2]

        if intermediateplots:
            show_fit_res(taus, crosscorrs[i], fitresults, decay='single')

    allresults = pd.DataFrame(data={'fit_A':A_lst, 'fit_gamma':gamma_lst, 'fit_gamma_err':gamma_err_lst, 'fit_bg':bg_lst})

    return allresults



def minimize_func(taus, ydata, decay='single', decayguess=0.1, decayguess2=0.1):
    import numpy as np
    from scipy.optimize import minimize
    from scipy.optimize import LbfgsInvHessProduct as invHessProd

    # This function takes the data and minimizes the difference between the fit and the data
    # Effectively, it minimizes errfunc_rel()
    # func is a string so that it can be a variable. Will be evaluated with an eval function

    if decay == 'single': # whether it is a single or double exponential
        func = 'exp_1'
        fit_guess = np.array([np.mean(ydata[:5]), decayguess, np.mean(ydata[-30:])])
        bound_guess = ((0.001, None), (0,None), (0, None)) # if the amplitude is a hard 0, it will throw a warning

    if decay == 'double':
        func = 'exp_2'
        fit_guess = np.array([np.mean(ydata[:5]), decayguess, np.mean(ydata[:10])*0.3, decayguess2, 10])
        bound_guess = ((0.001, None), (0, None), (0.001, None), (0,None), (0, None))
        
    if decay == 'triple':
        func = 'exp_3'
        fit_guess = [ydata[0], 0.65, np.mean(ydata[20:40]), 0.11, np.mean(ydata[100:150]), 0.017, 0.96]
        bound_guess = ((0.001, None), (0, None), (0.001, None), (0,None), (0.001, None), (0,None), (0, None))
    
    if decay == 'stretched':
        func = 'exp_stretched'
        fit_guess = np.array([np.mean(ydata[:5]), 0.3, 1.05, 10])
        bound_guess = ((0.001, None), (0, None), (0.001, None), (0.001, None))

    my_args = (taus, ydata, func)


    # The minimize function uses fit_guess as the first argument of errfunc_rel.
    # args are the rest of the arguments of efffunc_rel
    fit = minimize(errfunc_rel, fit_guess, args=my_args, bounds = bound_guess)

    fit_result =  fit['x']
    se = np.sqrt(np.diag(invHessProd.todense(fit['hess_inv'])))#*fit_guess      #se = standard error

    res = {'fit_result':fit_result, 'fit_guess':fit_guess, 'se':se}

    return res


def errfunc_rel(fit_guess, xdata, ydata, func_name):
    import numpy as np

    # Merit function to be minimized given that
    # we are using maximum likelihood procedure for Poisson distributed data
    # Bazjer  Eur Biophys J 20, 247

    yfit = np.array(eval(func_name)(fit_guess, xdata))
    err = np.sum(- ydata*np.log(yfit) + yfit)
    return err


def show_fit_res(taus, lifetimedata, fitresults, decay='single'):
    # solely for plotting purposes, to assess fit results
    import numpy as np
    import matplotlib.pyplot as plt

    if decay == 'single': # whether it is a single or double exponential
        func = 'exp_1'

    if decay == 'double':
        func = 'exp_2'
        
    if decay == 'triple':
        func = 'exp_3'
    
    if decay == 'stretched':
        func = 'exp_stretched'
   
    yfit= np.array(eval(func)(fitresults['fit_result'], taus))
   
    plt.figure()
    plt.plot(taus, yfit, 'C1', label='fit')

    plt.scatter(taus, lifetimedata, label='data', s=5)
    plt.legend(prop={'size': 18})
    plt.yscale('log')
    plt.ylim(bottom=0.5)
    plt.xlabel('t (ns)', fontsize = 22)
    plt.ylabel('counts', fontsize = 22)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    #plt.title(myfile[:-8])
    #plt.show()

    times = []
    if decay == 'single':
        times = [1./fitresults['fit_result'][1]]
    if decay == 'double':
        times = [1./fitresults['fit_result'][1], 1./fitresults['fit_result'][3]]
    return(times)


def exp_1(L, x):
    # a single exponential, to be used as fit function
    import numpy as np

    [A, gamma, bg] = L
    return A*np.exp(-gamma * x) + bg


def exp_2(L, x):
    # a double exponential, to be used as fit function
    import numpy as np

    [A1, gamma1, A2, gamma2, bg] = L
    return A1*np.exp(-gamma1 * x) + A2*np.exp(-gamma2 * x) + bg

def exp_3(L, x):
    # a double exponential, to be used as fit function
    import numpy as np

    [A1, gamma1, A2, gamma2, A3, gamma3, bg] = L
    return A1*np.exp(-gamma1 * x) + A2*np.exp(-gamma2 * x) + A3*np.exp(-gamma3 * x) + bg

def exp_stretched(L, x):
    import numpy as np
    [A1, gamma, beta, bg] = L
    return A1*np.power(np.exp(-gamma*x),beta) + bg
