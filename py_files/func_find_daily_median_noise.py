import numpy as np

from func_find_mean_HS_noise import func_find_mean_HS_noise

# This routine finds the median noise of all spectra during this day.

#Inputs:
# spc_input             - input spectra, 3D: (mprof, nhts, ppts)
# Npts                  - number of points in spectrum

# Outputs:
# median_noise_power_dB

# updated: 23-November-2022
# ========================================================================

def func_find_daily_median_noise(spc_input, Nspc_input):

    #%% Cell: Determine size of spc_input

    m_profiles  = np.size(spc_input,0)
    nhts    = np.size(spc_input,1)
    ppts    = np.size(spc_input,2)

    spc_nos_mean    = np.ones(m_profiles * nhts) * np.nan
    initial_seed_length  = ppts // 8

    #%% Cell: Process each spectrum
    
    # Process each profile

    for r in range(m_profiles):
        
        # process each range gate of this profile
        for c in range(nhts):
      
            # get the spectrum at this (profile, range)
            # Get spectra for this range gate
            spc_lin = np.squeeze(spc_input[r, c, :])
            
            # if not full of nans, find the noise level for this spectrum
            f   = ~np.isnan(spc_lin)
            
            if(np.sum(f) > 1):
                
                # Set really small values to 1/8'tile value
                sort_spc    = np.sort(spc_lin)
                sort_spc[0:initial_seed_length]  = np.ones(initial_seed_length) * sort_spc[initial_seed_length]
   
                # Perform the HS noise search
                [nos_mean, nos_max] = func_find_mean_HS_noise(
                    sort_spc, Nspc_input, initial_seed_length)
                
                spc_nos_mean[(r)*nhts + c] = nos_mean

    #%% Cell: Find the median of all noise estimates
                
    f  = ~np.isnan(spc_nos_mean)

    if(np.sum(f) >2):
        #sort_nos_mean = np.sort(spc_nos_mean[f])
        #median_nos_mean = np.median(sort_nos_mean)    
        
        median_nos_mean = np.median(spc_nos_mean[f])    
        
    #%% Cell: Convert to power (noise_level*Npts) & express in dB
        
    median_noise_power_dB    = 10.0*np.log10(median_nos_mean * ppts)

    #%% Cell: Define the variable to return to calling routine
    
    return [median_noise_power_dB]
