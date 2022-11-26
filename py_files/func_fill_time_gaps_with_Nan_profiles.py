# % This routine places a NaN profile when the time between profiles is too
# % large.

# % if there are time-gaps in profiles, the plot will have
# % streaks. Put a profile of NaNs between any profiles that are
# % too far apart

# % Inputs
# % input_time      - time vector for plotting
# % input_data_1    - 2D matrix of data to be plotted
# % input_data_2    - 2D matrix of data to be plotted
# % input_data_3    - 2D matrix of data to be plotted

# % Outputs
# % output_time      - updated time vector for plotting
# % output_data_1    - updated 2D matrix of data to be plotted
# % output_data_2    - updated 2D matrix of data to be plotted
# % output_data_3    - updated 2D matrix of data to be plotted

# % updated: 20-June-2021
# % ========================================================================

from matplotlib.pyplot import axis
import numpy as np
from numpy.lib import median


def func_fill_time_gaps_with_Nan_profiles(input_time, input_data1, input_data2, input_data3):

    # Define output variables
    output_time = input_time
    output_data1 = input_data1
    output_data2 = input_data2
    output_data3 = input_data3

    # Nan insert
    nan_insert1 = np.ones([np.size(output_data1, axis=1)]) * np.nan
    nan_insert2 = np.ones([np.size(output_data2, axis=1)]) * np.nan
    nan_insert3 = np.ones([np.size(output_data3, axis=1)]) * np.nan

    # Determine if any NaN profiles are needed
    n = len(input_time)
    dif_time = input_time[1:] - input_time[0:-1]
    median_dif_time = np.median(dif_time)
    time_thres = 3*np.median(dif_time)

    # Are there any time steps greater than 3*expected time gap?
    f = dif_time > time_thres
    if(np.sum(f) > 0):
        put_NaN_profile_flag = 1
    else:
        put_NaN_profile_flag = 0

    # Pad extra profiles if needed
    if(put_NaN_profile_flag):
        # f points to profiles with gaps > time_thres
        # get these indices
        gap_profile_index = np.argwhere(f)

        # Start from the end and go backwards-
        # put a profile of NaN into the datasets
        num_NaN_profile = np.sum(f) - 1

        for i in range(num_NaN_profile, -1, -1):
            # Get the index of last good profile
            gap_index = gap_profile_index[i]

            # Insert NaN profile at each gap index
            output_time = np.insert(
                output_time, gap_index+1, median_dif_time + output_time[gap_index])
            output_data1 = np.insert(
                output_data1, gap_index+1, nan_insert1, axis=0)
            output_data2 = np.insert(
                output_data2, gap_index+1, nan_insert2, axis=0)
            output_data3 = np.insert(
                output_data3, gap_index+1, nan_insert3, axis=0)

    return output_time, output_data1, output_data2, output_data3
