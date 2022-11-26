import numpy as np


# % This routine finds the noise adjusted snr and reflectivity with PRC = 1.
# % Estimate the reflectivity for a gain of PRC = 1;
# % Z = 10*log10(PRC * r^2 * SNR_linear)
# % Z = 10*log10(PRC) + 20*log10(r) + 10*log10(SNR_linear)
# % Z = 10*log10(PRC) + 20*log10(r) + SNR_dB
# %    if PRC = 1, then:
# % Z = 0 + 20*log10(r) + SNR_dB
# % The radar constant is given in dB, so, reflectivity is:
# % Z = radar_const_dB + 20*log10(r) + SNR_dB;

# % The SNR_dB needs to be adjusted by the mean noise observed at other hts:
# % SNR_dB_adj = SNR_dB_old + noise_old(in dB) - noise_new(in dB);

# % updated:05-Feb-2014


def func_find_adjusted_zdb_snr_with_const_noise(snr_input, nos_input, range_input, radar_const_dB, nos_const_dB):

    ## Define output variables ##
    n = np.size(snr_input, axis=0)
    zdb_adj = np.empty([n, 1])*np.nan
    snr_adj = np.empty([n, 1])*np.nan

    # Define the range in log ##
    range_20log10r = 20*np.log10(range_input)

    # Use the constant noise power that was used at input
    #nos_mean_adj = nos_const_dB

    ## Find the new adjusted snr and Z (with PRC = 1) ##
    # proces each range
    for i in range(n):
        # adjust the SNR_dB
        snr_adj[i] = snr_input[i] + nos_input[i] - nos_const_dB

        # find the reflectivity with PRC = 1
        zdb_adj[i] = radar_const_dB + range_20log10r[i] + snr_adj[i]

    return zdb_adj, snr_adj, nos_const_dB
