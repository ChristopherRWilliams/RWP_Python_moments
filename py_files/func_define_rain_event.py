#%% Cell: Define the imports and needed functions

import numpy as np
#import matplotlib.pyplot as plt

#from func_fill_time_gaps_with_Nan_profiles import func_fill_time_gaps_with_Nan_profiles

#%% Cell: Start the function

def func_define_rain_event(
            snr_min, range_min, range_max, uaz_lo, moment_lo):


    #%% Cell: Define the initial height range to look for rain

    # Get the snr for all valid range gates
    index_range_min = np.argmin(np.abs(uaz_lo.range - range_min))
    index_range_max = np.argmin(np.abs(uaz_lo.range - range_max))

    meso_ref_snr = moment_lo.snr[:, index_range_min:index_range_max+1]

    # How many profiles and range gates are in this limited set?
    [m, n] = np.shape(meso_ref_snr)

    # Initialize the profile indices to opposite values
    start_prof = m-1
    end_prof = 1

    #%% Count the number of profiles with snr > threshold
    
    # found valid number of pixels in each profile
    f = meso_ref_snr >= snr_min
    sum_valid_obs = np.sum(f, 1)  # count profiles

    # set the threshold to 1/3 of the number of gates
    threshold = np.ceil(n/4)

    #%% Find the begining and end times of the valid profiles
    
    for i in range(m):
        if (sum_valid_obs[i] > threshold):

            # find first index
            if(i < start_prof):
                start_prof = i
            # end if(i < start_prof):

            # find last index
            if(i > end_prof):
                end_prof = i
            # end if(i > end_prof):
        # end if (sum_valid_obs[i] > threshold):
    # end for i in range(m) loop
    
    #%% define the output values
    
    event_start_hour    = uaz_lo.timestamp[start_prof, 4]
    event_end_hour      = uaz_lo.timestamp[end_prof, 4] + 1

    if(event_start_hour < event_end_hour):
        event_rain_flag = 1
    else:
        event_rain_flag = 0
    # end if(event_start_hour < event_end_hour):


    #%% Return variables to calling program
    
    return event_rain_flag, event_start_hour, event_end_hour

# end of func_define_rain_event