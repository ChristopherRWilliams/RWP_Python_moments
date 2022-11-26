import numpy as np

# Estimate the noise floor for a single spectrum

# % Input:
# % spk           - Values from a single spectrum.  Linear and positive.
# %
# % Output:
# % mean_noise    - Mean of all values in noise.
# % max_noise     - Maximum value of noise values.
# %
# % From Tony's H-S method.
# % Updated: 26-March-2014
# %
# % Testing variables:
# % spk = data_lin(j,:);
# % nfft = NumSpectraAveraged(i);
# %
# % The spectrum may have negative values.  If negative values exist,
# % then shift spectrum so that min value is zero.


def func_find_mean_HS_noise(spk, nfft, initial_seed_length):

    min_spk = min(spk)
    if (min_spk < 0):
        spk = spk + np.absolute(min_spk)

    # Sort spectrum into accending order
    sort_spk = np.sort(spk)
    rtest = (1+nfft)/nfft
    xsum = 0
    xsum_sq = 0

    # Pre-load the xsum and xsum_sq values
    # to process the first few samples without testing for signal
    n_load = initial_seed_length

    for i in range(n_load):
        xsum = xsum + sort_spk[i]
        xsum_sq = xsum_sq + sort_spk[i]*sort_spk[i]
        #print(f'loop index = {i}')

    # point to the next value
    i = n_load

    noise_index = len(spk) - 1
    get_more = 1

    while(get_more):
        xsum = xsum + sort_spk[i]
        xsum_sq = xsum_sq+sort_spk[i]*sort_spk[i]

        if(xsum_sq*(i+1) > rtest*xsum*xsum):  # i+1 for 0 index offset
            noise_index = i-1
            get_more = 0

            if(noise_index < 1):
                noise_index = 1

        i = i + 1
        if(i > noise_index):
            get_more = 0
    # Save the mean and max noise values for output
    # Remember to shift the noise values beacuse the spectrum shifted
    if(min_spk < 0):
        mean_noise = np.mean(sort_spk[0:noise_index]) - np.absolute(min_spk)
        max_noise = sort_spk[noise_index] - np.absolute(min_spk)
    else:
        mean_noise = np.mean(sort_spk[0:noise_index+1])
        max_noise = sort_spk[noise_index]

    return [mean_noise, max_noise]
