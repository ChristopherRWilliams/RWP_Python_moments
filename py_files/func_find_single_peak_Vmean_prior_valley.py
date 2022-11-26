import numpy as np


# %function [spk_single, index_limit] = find_single_peak_Vmean_prior(spk, peak_index)

# % Find the largest magnitude peak in spectra.

# % Input Variables
# % spk             - input spectrum, no weighting
# % valley_thres    - depth of valley threshold

# % Output Variables
# % spk_clean     the same as input spk except outside points are NaN'd
# % index_limit   End index values for the found peak
# %               index_limit(1) -first velocity bin used in peak
# %               index_limit(2) -last veocity bin used in peak

# % This routine finds the single peak by finding the maximum value and
# %   marching down both sides of the peak until a NaN is found.
# % The new index is used to filter the spectrum

# % Steps:
# % 1. Check to see that there are valid points in spectra
# % 2. Find the maximum value in the spectra.
# % 3. Move to the left of the max value
# % 4. Move to the right of the max value
# % 5. NaN out all values outside of the indices


def func_find_single_peak_Vmean(spk, valley_thres):

    ### Check to see that there are valid points in spectra ###
    f = ~np.isnan(spk)

    if np.sum(f) < 2:
        spk_single = spk
        index_limit = np.array([np.nan, np.nan])
        return np.nan, index_limit
    ### Find max value in spectra ###

    i_max = np.nanargmax(spk)
    index_limit = np.array([0, 0], dtype=int)
    ### Move to the left of max value ###

    current_index = i_max
    current_value = spk[current_index]
    next_index = current_index
    get_next = 1

    while(get_next):

        # decrement index
        next_index = next_index - 1

        if(next_index < 1):
            next_index = 0
            current_index = 1
            get_next = 0
            break

        # Get the value of next spectral point
        next_value = spk[next_index]

        # Is the next_value below the noise threshold?
        if(np.isnan(next_value)):
            get_next = 0
            # point to this noise value as we want to know the noise
            current_index = next_index+1
        elif(next_value < current_value):  # is the next value decreasing or increasing
            # is the next value < the current value
            current_value = next_value
            current_index = next_index
        else:
            # Is the next value going up the next peak?
            delta_value = 10*np.log10(next_value) - 10*np.log10(current_value)
            if(delta_value >= valley_thres):
                get_next = 0
    # Place the current_index into the output variable
    index_limit[0] = current_index

    ### Move to the right of max value ###

    current_index = i_max
    current_value = spk[current_index]
    next_index = current_index
    get_next = 1

    while(get_next):
        # increment index
        next_index = next_index + 1

        if(next_index >= np.size(spk) - 1):
            next_index = np.size(spk)
            current_index = np.size(spk) - 1
            get_next = 0
            break

        # Get the value of the next spectral point
        next_value = spk[next_index]

        # Is the next value below noise threshold?
        if(np.isnan(next_value)):
            get_next = 0
            current_index = next_index-1
        elif(next_value < current_value):
            # Is the next value less than current value?
            current_value = next_value
            current_index = next_index
        else:
            # Is the next value going up the next peak?
            delta_value = 10*np.log10(next_value) - 10*np.log10(current_value)
            if(delta_value >= valley_thres):
                get_next = 0
        # Place current_index into output variable
    index_limit[1] = current_index

    # NaN out all values outside of the indices
    spk_single = spk
    # NaN out values to the left
    if(index_limit[0] > 1):
        end_index = index_limit[0]-1
        spk_single[0:end_index+1] = np.nan

    # look to right of peak
    if(index_limit[1] < np.size(spk)):
        start_index = index_limit[1]+1
        spk_single[start_index:np.size(spk_single)+1] = np.nan

    return [spk_single, index_limit]
