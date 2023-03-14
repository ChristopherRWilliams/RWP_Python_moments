
import matplotlib.pyplot as plt

import numpy as np
import warnings

from func_find_mean_HS_noise import func_find_mean_HS_noise
from func_find_noise_adjusted_zdb_and_snr_with_const_noise import \
    func_find_adjusted_zdb_snr_with_const_noise
from func_find_single_peak_Vmean_prior_valley import func_find_single_peak_Vmean

#%% Define the output variables

# This routine finds the moments of each spectrum

# % Inputs
# % spc_input               = spectrum (profiles, hts, npts)
# % Vd_input                = Doppler velocity;
# % range_input                = prof_center_ht;
# % Ncoh_input              = number_coherent_integration ;
# % Nspc_input              = number of individual spectra per average spectrum
# % valid_pts_thres         = valid threshold
# % cal_input               = RWP_radar_const_dB
# % nos_cont_dB             = noise

# % Outputs
# % .snr_pop         = snr from stright POP processing
# % .snr             = snr after adjusting for mean noise profile (best value)
# % .nos_power       = noise power after adjusting, should be equal to nos_cont_dB
# % .zdb             = reflectivity
# % .nos_pop         = noise power from POP processing: nos_mean * Npts
# % .Vmean           = mean Velocity
# % .Vsig            = Velocity standard deviation
# % .index_lt        = index of Vd for left valid point
# % .index_rt        = index of Vd for right valid point
# % .Vleft           = velocity of left valid point in spectrum
# % .Vright          = velocity of right valid point in spectrum
# % .nos_mean        = mean noise from HS method
# % .nos_max         = max noise value from HS method
# % .nos_std         = STD of noise samples at and below spc_nos_max

# Updated: 18 June 2021
# % Note about dealeasing of radial velocities...
# % Usually, the radar operates so that it observes noise samples above
# % the highest expected range with signal. At these higher heights, the
# % hydrometeor motion typically approaches zero. But Since the short pulse
# % mode (and the wind mode) does not observe above convective clouds,
# % there is no guarentee that the hydrometeor motion will be zero near the
# % top of range gates. So, for each profile, the processing starts at the
# % bottom and progresses up the profile.

# % This code implements the first of two aliasing algorithms. The first
# % algorithm is applied to a profile and tries to follow a continuous track
# % up (or down) the profile. The second aliasing algorithm will be a
# % temporal routine that is applied after all of the profiles for the day
# % have been processed.

# % The profile aliasing routine uses the variable Vmean_prior to represent
# % the prior Vmean estimate and is used to set boundaries on where to look
# % for Vmean in the current range gate.

# % Since the original spectrum is replicated to make a synthetic spectrum
# % 2-spectra wide, the estimated spectrum SNR and spectrum width will be
# % calculated correctly. The Vmean and the intigration limits may be shifted
# % by (+/-) 2*V_Nyquist.


def func_find_2spc_tda_moments(spc_input, Vd_input, range_input, Ncoh_input,
                               Nspc_input, valid_pts_thres, cal_input_dB, nos_const_dB):

    # define moment class
    class moment:
        def __init__(self, snr_pop, snr, nos_power, zdB, nos_pop, Vmean, Vsig, Vskew, Vkurt, index_lt, index_rt, Vleft, Vright, nos_mean, nos_max, nos_std):
            self.snr_pop = snr_pop
            # snr from stright POP processing
            self.snr = snr
            # snr after adjusting for mean noise profile
            self.nos_power = nos_power
            # noise after adjusting
            self.zdB = zdB
            # reflectivity
            self.nos_pop = nos_pop
            # noise from POP processing
            self.Vmean = Vmean
            # mean velocity
            self.Vsig = Vsig
            # spectrum standard deviation
            self.Vskew = Vskew
            # spectrum skewness
            self.Vkurt = Vkurt
            # spectrum kurtosis
            self.index_lt = index_lt
            # index of doppler velocity left valid point
            self.index_rt = index_rt
            # index of doppler velocity right valid point
            self.Vleft = Vleft
            # doppler velocity left valid point
            self.Vright = Vright
            # doppler velocity right valid point
            self.nos_mean = nos_mean
            # mean noise
            self.nos_max = nos_max
            # max noise
            self.nos_std = nos_std
            # noise standard deviation samples at or below nos_max

    #%% Define outputs
    
    m_profiles = np.size(spc_input, 0)
    nhts = np.size(spc_input, 1)

    spc_snr_pop     = np.ones([m_profiles, nhts]) * np.nan
    spc_snr         = np.ones([m_profiles, nhts]) * np.nan
    spc_nos_power   = np.ones([m_profiles, nhts]) * np.nan
    spc_zdb         = np.ones([m_profiles, nhts]) * np.nan
    spc_nos         = np.ones([m_profiles, nhts]) * np.nan
    spc_Vmean       = np.ones([m_profiles, nhts]) * np.nan
    spc_Vsig        = np.ones([m_profiles, nhts]) * np.nan
    spc_Vskew       = np.ones([m_profiles, nhts]) * np.nan
    spc_Vkurt       = np.ones([m_profiles, nhts]) * np.nan
    spc_index_lt    = np.ones([m_profiles, nhts]) * np.nan
    spc_index_rt    = np.ones([m_profiles, nhts]) * np.nan
    spc_Vleft       = np.ones([m_profiles, nhts]) * np.nan
    spc_Vright      = np.ones([m_profiles, nhts]) * np.nan

    spc_nos_mean    = np.ones([m_profiles, nhts]) * np.nan
    spc_nos_max     = np.ones([m_profiles, nhts]) * np.nan
    spc_nos_std     = np.ones([m_profiles, nhts]) * np.nan

    #%% Define Vd for the expanded spectrum

    Vd_orig = Vd_input
    Vd = np.transpose(Vd_orig)  # column form
    dVd = Vd[6] - Vd[5]

    # Expand velocity bins to have length 2*Vd
    Npts_expand = np.size(Vd,0) * 2
    Vd = np.zeros(Npts_expand)
    
    # The DC point is in location (Npts_expand//2)

    # Relabel the velocity bins. Start from zero and move in neg direction
    for c in range(Npts_expand//2-1, -1, -1):
        Vd[c] = Vd[c+1] - dVd

    # Relabel the velocity bins. Start from zero and move in pos direction
    for c in range(Npts_expand//2+1, Npts_expand):
        Vd[c] = Vd[c-1] + dVd

    # Define the left and right indices for the original spectra
    # Original spectrum is in the middle
    Npts_orig = np.size(Vd_orig)
    # [first set orignal index, second set orignal index]
    #index_orig = np.array([Npts_orig, 2*Npts_orig])
    #index_dc = Npts_expand/2

    # index to original DC value (zero velocity)
    orig_spc_index_dc = Npts_orig//2

    # Define the indices to produce the expanded spectrum.
    # The expanded spectrum has the size 2*length(orig_spc):
    # The expansion takes the right half and appends it to the left side, and
    # appends the left have to rhte right side.
    # spc_expand = [right_half, full spectrum, left_half]

    # seg1_indices = DC to Npts
    # seg2_indices = 1  to Npts
    # seg3_indices = 1 to (DC-1)

    # These indices are referenced to the original spectrum.
    orig_spc_seg1_indices   = np.zeros(Npts_orig//2,int)
    for c in range(0, Npts_orig//2):
        index = Npts_orig//2 + c 
        orig_spc_seg1_indices[c] = index

    orig_spc_seg2_indices   = np.zeros(Npts_orig,int)
    for c in range(0, Npts_orig):
        index = c 
        orig_spc_seg2_indices[c] = index

    orig_spc_seg3_indices   = np.zeros(Npts_orig//2,int)
    for c in range(0, Npts_orig//2):
        index = c 
        orig_spc_seg3_indices[c] = index
        

    # Define interval for finding peak in expanded spectrum
    peak_interval_step   = Npts_orig//2 + Npts_orig//4


    #%% Rename input variables to local variable names
    
    Ncoh = Ncoh_input
    Nspc = Nspc_input

    #%%  Calcualte TDA Correction
    
    # The index needs to span all of expanded Vd.
    # The index needs to span from -Npts_orig to +Npts_orig
    # The TDA correction is non-zero if NCOH > 1
    index = np.arange(-1*Npts_orig, Npts_orig)

    # find x term
    x = (np.pi * index / (Ncoh*Npts_orig))

    # find y term
    y = (np.pi * index/Npts_orig)

    # There is a discontinutity at y(1), such that sin(y(1)) = 0.
    # Then, the tda correction will have a division by zero.
    # set y[0]) to the next index value:
    y[0] = y[1];
    # Also, y(DC) = 0, which will cause a division by zero.
    # Set y(DC) to the next value (which will cause the TDA correction to be 
    # approximately equal to 1.
    y[Npts_orig] = y[Npts_orig + 1];

    # calculate numerator and denominator for TDA correction
    numerator = (Ncoh**2)*(np.sin(x)*np.sin(x))
    denominator = (np.sin(y)*np.sin(y))

    # calculate TDA correction
    f = (numerator > 0)
    tda_2spc_correction = numerator * np.nan
    tda_2spc_correction[f] = numerator[f]/denominator[f]
    if(np.sum(~f) > 0):
        tda_2spc_correction[~f] = np.sum(~f)

    # start_time = time.time()  # Time per profile check
    # time_list = []

    #%% Process each Profile
    
    for prof_num in range(m_profiles):
        if(prof_num == 0):
            print(
                f'Processing profile {prof_num+1} out of {m_profiles} profiles')
        else:
            if(np.remainder(prof_num, 250) == 0):
                print(
                    f'Processing profile {prof_num} out of {m_profiles} profiles')
                # stop_time = time.time()-start_time  # debug excess time
                # time_list.append(stop_time)
                #print(f'Time to process 250 profiles is {stop_time} seconds')
                # start_time = time.time()  # reset timer

        # Start from the bottom and find the moments of each spectrum
        # Use the previous average Vmean estimates for inital guess

        # Set Vmean_prior to 0
        Vmean_prior = 0
        # Save the Vmean_prior for debugging plots
        Vmean_prior_prof = np.ones([nhts, 1]) * np.nan

        # Set the number of gates to average for estimating Vmean_prior
        drange = range_input[6] - range_input[5]
        Vmean_prior_ave_depth = 2*drange
        num_Vmean_prior_gates = round(Vmean_prior_ave_depth/drange)

        # define snr and Vsig thresholds for which samples to estimate Vmean_prior
        snr_Vmean_prior_threshold = -25  # so low all samples are used
        Vsig_Vmean_prior_threshold = 2*dVd

        for c in range(nhts):

            # Get spectra for this range gate
            spc_lin = np.squeeze(spc_input[prof_num, c, :])

            # Replace DC valvue with mean of neighbors
            spc_lin[orig_spc_index_dc] = .5 * \
                (spc_lin[orig_spc_index_dc-1]+spc_lin[orig_spc_index_dc+1])

            # Expand spectra to correct for aliasing
            spc = np.concatenate((spc_lin[orig_spc_seg1_indices], spc_lin[orig_spc_seg2_indices], spc_lin[orig_spc_seg3_indices]), axis=0)

            # Find noise statistics of original spectrum: mean, max, and std ###
            f = ~np.isnan(spc)
            if(np.sum(f) > 0):
                # Get the original central spectrum
                #spc_orig = spc[index_orig[0]:index_orig[1]]

                # Perform the HS noise search
                initial_seed_length = len(spc_lin)//8
                spc_sort    = np.sort(spc_lin)
                spc_sort[0:initial_seed_length] = np.ones(initial_seed_length) * spc_sort[initial_seed_length-1]
                
                [nos_mean, nos_max] = func_find_mean_HS_noise(
                    spc_sort, Nspc, initial_seed_length)

                # find the noise standard deviation
                f = spc_lin < nos_max
                if (np.sum(f) > 3):
                    nos_std = np.std(spc_lin[f], ddof=1)
                else:
                    nos_std = np.nan

                # Set noise threshold
                std_nos_threshold = nos_mean + 3.0*nos_std
                nos_threshold = np.maximum(nos_max, std_nos_threshold)

                spc_nos_mean[prof_num, c] = nos_mean
                spc_nos_max[prof_num, c] = nos_max
                spc_nos_std[prof_num, c] = nos_std
                
            #%% Use Vmean_prior to isolate the valid peak (of two options)

            # find the peak in the original spectrum
            orig_peak_index     = np.argmax(spc_lin)
            
         
            # This peak occurs in two locations in the expanded spectrum
            # It depends on whether orig_peak_index is the right or left of
            # DC, to determine left and right indices.
            if (orig_peak_index < orig_spc_index_dc):
                left_peak_index     = orig_peak_index + Npts_orig//2
                right_peak_index    = left_peak_index + Npts_orig
            else:
                right_peak_index     = orig_peak_index + Npts_orig//2
                left_peak_index      = right_peak_index - Npts_orig
         
            # find the index corresponding to Vmean_prior
         
            # save the Vmean_prior for debugging plots
            Vmean_prior_prof[c] = Vmean_prior

            # Find index for Vmean_prior
            if ~np.isnan(Vmean_prior):
                peak_index = np.argmin(np.absolute(Vd - Vmean_prior))
            else:
                peak_index = np.argmin(np.absolute(Vd - 0))

            # Is right_peak_index or left_peak_index closer to peak_index?
            dif_right   = np.abs(right_peak_index - peak_index)
            dif_left    = np.abs(left_peak_index  - peak_index)         

            # Keep a limited range of spectra centered on Vmean_prior
            # left_index and right_index point to valid values
            # The peak_interval_step is (3/4)*Npts_orig, because a step of
            # Npts_orig will be the Nyquist peak.
            if(dif_right < dif_left):
                # adjust the peak_index to the true peak
                peak_index  = right_peak_index
                left_index  = peak_index - peak_interval_step
                right_index = Npts_expand
            else:
                peak_index  = left_peak_index
                left_index  = 0
                right_index = peak_index + peak_interval_step
         
            # make sure left_index and right_index are valid
            if left_index > 0:
                spc[:left_index][:] = nos_mean

            if right_index < Npts_expand:
                spc[right_index + 1:] = nos_mean

            #%% Find the single peak

            f = spc < nos_threshold
            spc[f] = spc[f]*np.nan

            valley_thres = 100

            [_, index_limit] = func_find_single_peak_Vmean(
                spc, valley_thres)

            #%% Get the spectrum within index_limit and do TDA correction 
            
            # make sure there are enough consecutive points above noise floor
            num_pts = index_limit[1] - index_limit[0] + 1

            if(num_pts >= valid_pts_thres):
                # Construct short arrays

                spk_short = spc[index_limit[0]:index_limit[1]+1]

                Vd_short = Vd[index_limit[0]:index_limit[1]+1].flatten()
                # flatten for element wise multiply below

                tda_short = tda_2spc_correction[index_limit[0]:index_limit[1]+1]

                # Replace the values with the tda values
                spk_tda_short = (spk_short - nos_mean) * tda_short + nos_mean

                # Estimate the moments using the TDA values
                noise_power = nos_mean * Npts_orig

                # ignore negative noise power
                if noise_power < 0:
                    warnings.filterwarnings("ignore", category=RuntimeWarning)
                    # this flag will be tripped anytime a complex log is calculated (line 393 and 394)

                if(np.sum(spk_tda_short-nos_mean) > 0):
                    snr_tda = 10 * \
                        np.log10(np.sum(spk_tda_short-nos_mean)/noise_power)
                    nos_tda = 10*np.log10(noise_power)
                    Vmean_tda = np.sum(spk_tda_short*Vd_short)/np.sum(spk_tda_short)
                    Vsig_tda = np.sqrt(
                        np.sum((Vd_short-Vmean_tda)**2*spk_tda_short)/np.sum(spk_tda_short))

                    # % Calculate the skewness
                    # % skewness is: sum((Vd_short - Vmean_tda).^3 .* spk_tda_short) /sum(spk_tda_short)
                    # %              ----------------------------------------------
                    # %                 Vsig_tda.^3
            
                    skewness_top      = (sum((Vd_short - Vmean_tda)**3 * spk_tda_short))/ sum(spk_tda_short);
                    skewness_bot      = Vsig_tda**3;
                    Vskew_tda         = skewness_top / skewness_bot;
            
                    # % Calculate the kurtosis
                    # % kurtosis is: sum((Vd_short - Vmean_tda).^4 .* spk_tda_short) /sum(spk_tda_short)
                    # %              ----------------------------------------------
                    # %                 Vsig_tda.^4
            
                    kurtosis_top      = (sum((Vd_short - Vmean_tda)**4 * spk_tda_short))/ sum(spk_tda_short);
                    kurtosis_bot      = Vsig_tda**4;
                    Vkurt_tda         = kurtosis_top / kurtosis_bot;
            

                # Save the values if there are any valid moments
                if(np.isfinite(snr_tda)) and (np.isreal(snr_tda)):
                    spc_snr_pop[prof_num, c]    = snr_tda
                    spc_nos[prof_num, c]        = nos_tda
                    spc_Vmean[prof_num, c]      = Vmean_tda
                    spc_Vsig[prof_num, c]       = Vsig_tda
                    spc_Vskew[prof_num, c]      = Vskew_tda
                    spc_Vkurt[prof_num, c]      = Vkurt_tda
                    spc_index_lt[prof_num, c]   = index_limit[0]
                    spc_index_rt[prof_num, c]   = index_limit[1]
                    spc_Vleft[prof_num, c]      = Vd[index_limit[0]]
                    spc_Vright[prof_num, c]     = Vd[index_limit[1]]

                #%% Find new Vmean_prior value
                
                if (~np.isnan(spc_Vmean[prof_num, c])):

                    # get range gate index for previous valid Vmean estimate
                    gate_indices = np.arange(
                        c - (num_Vmean_prior_gates-1), c+1)
                    f = gate_indices >= 0
                    gate_indices = gate_indices[f]

                    # get moments for these range gates
                    array_Vmean = spc_Vmean[prof_num, gate_indices]
                    array_snr = spc_snr_pop[prof_num, gate_indices]
                    array_Vsig = spc_Vsig[prof_num, gate_indices]
                    # apply thresholds to keep only valid moments
                    f = ~np.isnan(array_Vmean)
                    g = array_snr > snr_Vmean_prior_threshold
                    k = array_Vsig > Vsig_Vmean_prior_threshold
                    # Check for enough samples
                    if(np.sum(f & g & k)) >= num_Vmean_prior_gates//2:
                        Vmean_prior = np.median(array_Vmean[f & g & k])

        #%% Estimate the reflectivity

        # Now that you have the profile, do noise adjustment
        # Adjust SNR by the noise and calculate reflectivity

        snr_profile = spc_snr_pop[prof_num, :]
        nos_profile = spc_nos[prof_num, :]
        ht_profile = range_input

        [zdb_adj, snr_adj, nos_mean_adj] = func_find_adjusted_zdb_snr_with_const_noise(
            snr_profile, nos_profile, ht_profile, cal_input_dB, nos_const_dB)

        # save values
        spc_snr[prof_num, :] = np.transpose(snr_adj)
        spc_zdb[prof_num, :] = np.transpose(zdb_adj)
        spc_nos_power[prof_num, :] = nos_mean_adj

        #%% Construct a spectrum profile (for testing only)

        plot_flag = 0

        if(plot_flag):

            plot_ht = np.array(range_input/1000)
            dht = plot_ht[6] - plot_ht[5]

            min_ht = 0
            max_ht = np.around(np.amax(plot_ht))

            plot_spc_lin = np.squeeze(spc_input[prof_num, :, :])

            plot_spk_lin = np.ones([nhts, 3*Npts_orig]) * np.nan
            plot_spk_dB = np.ones([nhts, 3*Npts_orig]) * np.nan

            for c in range(nhts):
                spc = np.concatenate((
                    plot_spc_lin[c, :], plot_spc_lin[c, :], plot_spc_lin[c, :]), axis=0)
                plot_spk_lin[c, :] = spc
                plot_spk_dB[c, :] = 10*np.log10(np.absolute(spc))

            # Get the moments
            plot_Vmean = spc_Vmean[prof_num, :]
            plot_Vmean_prior = Vmean_prior_prof
            plot_Vleft = spc_Vleft[prof_num, :]
            plot_Vright = spc_Vright[prof_num, :]

            # Bottom Row - Spectra

            f, (a0, a1) = plt.subplots(
                1, 2, gridspec_kw={'width_ratios': [3, 1]})
            plt.sca(a0)
            plt.pcolormesh(Vd-dVd/2,
                           plot_ht - dht/2, plot_spk_dB[:-1, :-1], shading='flat')
            a0.axis([-45, 45, min_ht, max_ht])

            plt.plot(plot_Vmean, plot_ht, 'bs', label='Vmean')
            plt.plot(plot_Vmean_prior, plot_ht, 'kx', label='Vprior')

            plt.plot(plot_Vleft, plot_ht, 'r+')
            plt.plot(plot_Vright, plot_ht, 'r+')
            plt.plot([Vd[0], Vd[0]], [0, max_ht], 'k', lw=1)
            plt.plot([Vd_orig[0], Vd_orig[0]], [0, max_ht], 'k', lw=1)
            plt.plot([Vd_orig[Npts_orig-1], Vd_orig[Npts_orig-1]],
                     [0, max_ht], 'k', lw=1)
            plt.plot([Vd[-1], Vd[-1]], [0, max_ht], 'k', lw=1)

            labels = ['', '40', '', '30', '', '20', '', '10', '',
                      '0', '', '10', '', '20', '', '30', '', '40', '']

            plt.xticks(np.arange(-45, 50, 5), labels)
            a0.legend()

            plt.xlabel('(upward) Velocity [m/s] (downward)')
            plt.ylabel('Height (km)')
            plt.title(f'SGP, RWP, Prof Num: {prof_num}')

            plt.sca(a1)
            plt.plot(spc_snr_pop[prof_num, :], plot_ht, 'k*-')
            plt.grid(True)
            plt.xlabel('SNR [dB]')
            plt.title('SNR')

            plt.tight_layout()
            plt.show()

    #%% Create moment object
    
    moment_object = moment(spc_snr_pop, spc_snr, spc_nos_power, spc_zdb, spc_nos, spc_Vmean, spc_Vsig,
                           spc_Vskew, spc_Vkurt, spc_index_lt, spc_index_rt, spc_Vleft, spc_Vright, spc_nos_mean, spc_nos_max, spc_nos_std)

    return moment_object
