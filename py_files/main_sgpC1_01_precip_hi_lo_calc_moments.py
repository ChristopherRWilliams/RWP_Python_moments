# Program = main_sgpC1_01_precip_hi_lo_calc_moments.py
# updated = 25-November-2022

# This routine is specific for site SGP C1, near Lamont, Oklahoma.

# This routine reads raw precip mode RWP spectra, divides the spectra into
# to modes (lo and hi), and calculates the moments for each profile.
# **********************************************************************

#%% Cell: Input Definitions

import datetime
# import glob
# import os

#import matplotlib.pyplot as plt
import numpy as np

from func_daterange import func_daterange
# from func_fill_time_gaps_with_Nan_profiles import \
#   func_fill_time_gaps_with_Nan_profiles
from func_find_2spc_tda_moments_of_spc import func_find_2spc_tda_moments
from func_read_multiple_netCDF_files_rwp_precip_spc_2modes import \
    func_read_multiple_netCDF_files_rwp_precip_spc_2modes
from func_time_dealias_on_profile import func_time_dealias_on_profile
from func_define_rain_event import func_define_rain_event
from func_plot_time_ht_3panel import func_plot_time_ht_3panel
from func_save_moments_in_netcdf import func_save_moments_in_netcdf
from func_find_daily_median_noise import func_find_daily_median_noise

#%% Cell: Define Input / Output directories and root filenames

# Define the base path and site specific variables

# Set the basepath to one level above "\py_files"
#basepath                            = 'E:\Projects\SGP\Python\site_sgp_I10'
basepath                            = 'C:\\Users\\Christopher\\Dropbox\\AMT_SGP_RWP_code'
site_id_3CHAR_str                   = 'SGP'
site_id_3char_str                   = 'sgp'
facility_id_str                     = 'C1'

# Define the radar operating frequency
rwp_freq_MHz                        = 915

# Set calibration constants to 0 dB (same as scalar = 1 in linear)
uaz_lo_cal_dB = 0  # in dB
uaz_hi_cal_dB = 0

# Define the constant noise values (site dependent)
# These two user defined constants are not needed: c.williams 2022-11-23
#lo_nos_const_dB                     = 75
#hi_nos_const_dB                     = 69

# Define location and name of input files
# Data on DataStore4
#input_netCDF_directory              = f'G:/DOE_Archive/{site_id_3char_str}/{facility_id_str}/rwp_precip_spec/'
input_netCDF_directory              = f'{basepath}\\rwp_precip_spec'
input_root_filename                 = f'{site_id_3char_str}{rwp_freq_MHz}rwpprecipspec{facility_id_str}.a0.'
#input_suffix                        = '.*.nc'
input_suffix                        = '.*.cdf'

# Define python formated output file (not used after 13-June-2022) 
data_level_str                      = 'crw0'
npz_output_directory                = f'{basepath}\\npz_precip_{data_level_str}'
#output_version                      = 'v20220218'
#output_suffix                       = '.nc'
npz_output_suffix                   = '.npz'

# Define netCDF formated output file 
nc_output_directory                 = f'{basepath}\\nc_precip_{data_level_str}'
#nc_lo_root_filename                = f'{site_3char_str}{site_id}rwppreciplo.{data_level}.'
#nc_hi_root_filename                = f'{site_3char_str}{site_id}rwppreciphi.{data_level}.'
nc_output_suffix                    = '.nc'

# This file type has both lo and high modes
output_lo_root_filename             = f'{site_id_3char_str}{facility_id_str}rwppreciplo.{data_level_str}.'
output_hi_root_filename             = f'{site_id_3char_str}{facility_id_str}rwppreciphi.{data_level_str}.'

# Define output image directory
image_day_directory                 = f'{basepath}\\images_precip_day/'
image_lo_root_ZVW_filename          = f'fig01_{site_id_3char_str}{facility_id_str}_rwp_precip_lo_ZVmeanVsig_day_'
image_lo_root_ZSK_filename          = f'fig02_{site_id_3char_str}{facility_id_str}_rwp_precip_lo_ZSkewKurt_day_'
image_hi_root_ZVW_filename          = f'fig03_{site_id_3char_str}{facility_id_str}_rwp_precip_hi_ZVmeanVsig_day_'
image_hi_root_ZSK_filename          = f'fig04_{site_id_3char_str}{facility_id_str}_rwp_precip_hi_ZSkewKurt_day_'
image_event_directory               = f'{basepath}\\images_precip_event/'
image_event_lo_root_ZVW_filename    = f'fig05_{site_id_3char_str}{facility_id_str}_rwp_precip_event_lo_ZVmeanVsig_day_'
image_event_lo_root_ZSK_filename    = f'fig06_{site_id_3char_str}{facility_id_str}_rwp_precip_event_lo_ZSkewKurt_day_'
image_event_hi_root_ZVW_filename    = f'fig07_{site_id_3char_str}{facility_id_str}_rwp_precip_event_hi_ZVmeanVsig_day_'
image_event_hi_root_ZSK_filename    = f'fig08_{site_id_3char_str}{facility_id_str}_rwp_precip_event_hi_ZSkewKurt_day_'

#%% Cell: Define the days to process

# Define the year month day to process

start_year  = 2018
end_year    = 2018
start_month =    6
end_month   =    6
start_day   =    7
end_day     =    7

# Verify dates are valid
try:
    start_date = datetime.date(start_year, start_month, start_day)
    end_date = datetime.date(end_year, end_month, end_day)
    delta = end_date - start_date
    num_days = delta.days + 1
except ValueError:
    print('...The end day is out of range for the month/months selected...')

#%% Cell: User should not need to change anything below this cell


#%% Cell: Loop through all days

# loop through all dates within defined date range above
for single_date in func_daterange(start_date, end_date):

    # Define date string
    # Output filename format: YearMonthDay 20190401
    date_str = single_date.strftime('%Y%m%d')
    
    # Plotting format: Year-Month-Day 2019-04-01
    date_plot_str = single_date.strftime('%Y-%m-%d')

    # print message to user
    print(' ')
    print(
        f'...starting to process day: {date_plot_str} by reading raw spectra files... ')
    
    # Read multiple input files for this day...
    uaz_lo, uaz_hi, process_spc_flag = func_read_multiple_netCDF_files_rwp_precip_spc_2modes(
            input_netCDF_directory, input_root_filename, single_date, 
            input_suffix, rwp_freq_MHz)
    
    #%% Cell: If there are any valid profiles, find their moments
    
    if(process_spc_flag):

        #%% Cell: Define the output filenames for this day
        
        # Set python formated output filenames
        npz_output_lo_filename = f'{npz_output_directory}\{output_lo_root_filename}{date_str}{npz_output_suffix}'
        npz_output_hi_filename = f'{npz_output_directory}\{output_hi_root_filename}{date_str}{npz_output_suffix}'

        # Set netCDF formated output filenames
        nc_output_lo_filename = f'{nc_output_directory}\{output_lo_root_filename}{date_str}{nc_output_suffix}'
        nc_output_hi_filename = f'{nc_output_directory}\{output_hi_root_filename}{date_str}{nc_output_suffix}'

        #%% Cell: Find the lo mode median noise value for this day
        
        print(' ')
        print('Estimating lo mode median noise for whole day...')
        
        # Process the uaz lo mode 
        spc_input       = uaz_lo.spc_pop
        Nspc_input      = uaz_lo.Nspc
        
        median_lo_noise_power_dB = func_find_daily_median_noise(spc_input,Nspc_input)
                
        #%% Cell: Process lo mode spectra to estimate moments
        
        # Calculate the moments of the uaz low mode
        # Define common input variables
        spc_input       = uaz_lo.spc_pop
        Vd_input        = uaz_lo.Vd
        range_input     = uaz_lo.range
        Ncoh_input      = uaz_lo.Ncoh
        Nspc_input      = uaz_lo.Nspc
        cal_input_dB    = uaz_lo_cal_dB
        valid_pts_thres = 3
        nos_const_dB    = median_lo_noise_power_dB
    
        # Process Lo moments        
        print(' ')
        print('Processing Moments in lo mode...')
        moment_lo = func_find_2spc_tda_moments(spc_input, Vd_input, range_input, Ncoh_input,
                                               Nspc_input, valid_pts_thres, cal_input_dB, nos_const_dB)

        #%% Cell: Do temporal deliasing on lo mode velocity

        # Process the lo mode
        Vmean_input    = moment_lo.Vmean
        index_lt_input = moment_lo.index_lt
        index_rt_input = moment_lo.index_rt
        Npts           = uaz_lo.Npts
        Ncoh           = uaz_lo.Ncoh
        IPP            = uaz_lo.ipp
               
        num_prior_profiles          = 1               
        wavelength                  = (3 * 10**8) / (915 * 10**6)
        V_Nyquist                   = (wavelength/4) * (1/(Ncoh * IPP))
        num_range_gate_threshold    = 5
        V_dif_threshold             = 1.5 * V_Nyquist

        print(' ')
        print('Removing aliased profiles in lo mode...')               
        Vmean_corrected, index_lt_corrected, index_rt_corrected = func_time_dealias_on_profile(
            Vmean_input, index_lt_input, index_rt_input, 
            V_Nyquist, Npts, num_prior_profiles, V_dif_threshold, num_range_gate_threshold)

        # Replace pre-dealiased values with corrected values
        moment_lo.Vmean      = Vmean_corrected
        moment_lo.index_lt   = index_lt_corrected
        moment_lo.index_rt   = index_rt_corrected

        #%% Cell: Save the lo mode moments - netCDF format
        
        # Global attributes in netCDF output file
        location_description_str    = 'DOE ARM Southern Great Plains (SGP), Lamont, Oklahoma, NW Radar Wind Profiler Intermediate Site (I10)';
        today                       = datetime.date.today()        
        history_str                 = 'created by Christopher Williams on ' + today.strftime("%Y.%m.%d") + ', using software version 2022.04.26'
        source_str                  = f'{input_netCDF_directory}{input_root_filename}YYYYMMDD.HHMMSS.nc'
        site_id_str                 = site_id_3char_str
        datastream_str              = f'{output_lo_root_filename}YYYYMMDD.nc'

        # rename the variables before writing to netCDF file
        #single dimension variables
        bad_flag                            = (-1)*(2**15 - 1)
        r_calib_radar_constant              = uaz_lo_cal_dB
        num_time_domain_integrations        = uaz_lo.Ncoh
        num_frequency_domain_integrations   = uaz_lo.Nspc
        num_fft_points                      = uaz_lo.Npts
        interpulse_period                   = uaz_lo.ipp
        num_code_bits                       = uaz_lo.code_bits
        reference_noise_power               = moment_lo.nos_power[0,0]
        operating_frequency_MHz             = rwp_freq_MHz
        nyquist_velocity                    = uaz_lo.Nyquist_vel
        # rename pulse_length to pulse_width
        pulse_width                         = uaz_lo.pulse_length
        range_resolution                    = uaz_lo.range_resolution
        lat                                 = uaz_lo.lat
        lon                                 = uaz_lo.lon
        alt                                 = uaz_lo.alt
               
        # calculate base time in Epoch
        base_time   = datetime.datetime(int(uaz_lo.timestamp[0,0]), int(uaz_lo.timestamp[0,1]), int(uaz_lo.timestamp[0,2]), 0, 0, tzinfo=datetime.timezone.utc).timestamp()
        
        # seconds into day
        time_sec       = uaz_lo.timestamp[:,4]*(60*60) + \
                  uaz_lo.timestamp[:,5] * (60) + uaz_lo.timestamp[:,6] + \
                  uaz_lo.timestamp[:,7] / 1000
        
        # (time_sec is both 'time_offset' and 'time' in netCDF file)
        time_offset    = time_sec
        time           = time_sec

        # range dimension variables
        range_along_beam            = uaz_lo.range
        height_above_radar          = uaz_lo.range

        # These values are fixed for vertical mode. They will be time dependent in the wind mode...
        beam_azimuth_angle          = 0
        beam_elevation_angle        = 90
               
               
        # (time,range) dimension variables
        signal_to_noise_ratio   = moment_lo.snr
        reflectivity_factor     = moment_lo.zdB
        radial_velocity         = moment_lo.Vmean
        spectrum_width          = (2) * moment_lo.Vsig
        spectrum_skewness       = moment_lo.Vskew
        spectrum_kurtosis       = moment_lo.Vkurt
        spectrum_mean_noise_level   = 10 * np.log10(moment_lo.nos_mean)

        # do the actual netCDF file writting                              
        func_save_moments_in_netcdf(nc_output_lo_filename, 
                  bad_flag, 
                  source_str, 
                  site_id_str, 
                  facility_id_str, 
                  data_level_str, 
                  location_description_str, 
                  datastream_str, 
                  history_str, 
                  base_time, 
                  time_offset, 
                  time, 
                  lat, 
                  lon, 
                  alt, 
                  r_calib_radar_constant, 
                  num_time_domain_integrations, 
                  num_frequency_domain_integrations, 
                  num_fft_points, 
                  interpulse_period, 
                  operating_frequency_MHz, 
                  nyquist_velocity, 
                  pulse_width, 
                  range_resolution, 
                  num_code_bits, 
                  reference_noise_power, 
                  beam_azimuth_angle, 
                  beam_elevation_angle, 
                  range_along_beam, 
                  height_above_radar, 
                  signal_to_noise_ratio, 
                  reflectivity_factor, 
                  radial_velocity, 
                  spectrum_width, 
                  spectrum_skewness, 
                  spectrum_kurtosis, 
                  spectrum_mean_noise_level, 
                  date_str)

        #%% Cell: Find the lo mode median noise value for this day
        
        print(' ')
        print('Estimating hi mode median noise for whole day...')
        
        # Process the uaz lo mode 
        spc_input       = uaz_hi.spc_pop
        Nspc_input      = uaz_hi.Nspc
        
        median_hi_noise_power_dB = func_find_daily_median_noise(spc_input,Nspc_input)

        #%% Cell: Process hi mode spectra to estimate moments

        # Calculate the moments of the uaz high mode
        # Define common input variables
        spc_input = uaz_hi.spc_pop
        Vd_input = uaz_hi.Vd
        ht_input = uaz_hi.range
        Ncoh_input = uaz_hi.Ncoh
        Nspc_input = uaz_hi.Nspc
        cal_input_dB = uaz_hi_cal_dB
        valid_pts_thres = 3
        nos_const_dB = median_hi_noise_power_dB

        # Process Hi moments
        print(' ')
        print('Processing Moments in hi mode...')
        moment_hi = func_find_2spc_tda_moments(spc_input, Vd_input, ht_input, Ncoh_input,
                                               Nspc_input, valid_pts_thres, cal_input_dB, nos_const_dB)

        #%% Cell: Do temporal deliasing on hi mode velocity

        # Process the lo mode
        Vmean_input    = moment_hi.Vmean
        index_lt_input = moment_hi.index_lt
        index_rt_input = moment_hi.index_rt
        Npts           = uaz_hi.Npts
        Ncoh           = uaz_hi.Ncoh
        IPP            = uaz_hi.ipp
               
        num_prior_profiles          = 1               
        wavelength                  = (3 * 10**8) / (915 * 10**6)
        V_Nyquist                   = (wavelength/4) * (1/(Ncoh * IPP))
        num_range_gate_threshold    = 5
        V_dif_threshold             = 1.5 * V_Nyquist

        print(' ')        
        print('Removing aliased profiles in hi mode...')               
        Vmean_corrected, index_lt_corrected, index_rt_corrected = func_time_dealias_on_profile(
            Vmean_input, index_lt_input, index_rt_input, 
            V_Nyquist, Npts, num_prior_profiles, V_dif_threshold, num_range_gate_threshold)

        # Replace pre-dealiased values with corrected values
        moment_hi.Vmean      = Vmean_corrected
        moment_hi.index_lt   = index_lt_corrected
        moment_hi.index_rt   = index_rt_corrected

        #%% Cell: Save the hi mode moments - netCDF format
        
        # Global attributes in netCDF output file
        location_description_str    = 'DOE ARM Southern Great Plains (SGP), Lamont, Oklahoma, NW Radar Wind Profiler Intermediate Site (I10)';
        today                       = datetime.date.today()        
        history_str                 = 'created by Christopher Williams on ' + today.strftime("%Y.%m.%d") + ', using software version 2022.11.23'
        source_str                  = f'{input_netCDF_directory}{input_root_filename}YYYYMMDD.HHMMSS.nc'
        site_id_str                 = site_id_3char_str
        datastream_str              = f'{output_hi_root_filename}YYYYMMDD.nc'

        # rename the variables before writing to netCDF file
        #single dimension variables
        bad_flag                            = (-1)*(2**15 - 1)
        r_calib_radar_constant              = uaz_hi_cal_dB
        num_time_domain_integrations        = uaz_hi.Ncoh
        num_frequency_domain_integrations   = uaz_hi.Nspc
        num_fft_points                      = uaz_hi.Npts
        interpulse_period                   = uaz_hi.ipp
        num_code_bits                       = uaz_hi.code_bits
        reference_noise_power               = moment_hi.nos_power[0,0]
        operating_frequency_MHz             = rwp_freq_MHz
        nyquist_velocity                    = uaz_hi.Nyquist_vel
        # rename pulse_length to pulse_width
        pulse_width                         = uaz_hi.pulse_length
        range_resolution                    = uaz_hi.range_resolution
        lat                                 = uaz_hi.lat
        lon                                 = uaz_hi.lon
        alt                                 = uaz_hi.alt
               
        # calculate base time in Epoch
        base_time   = datetime.datetime(int(uaz_hi.timestamp[0,0]), int(uaz_hi.timestamp[0,1]), int(uaz_hi.timestamp[0,2]), 0, 0, tzinfo=datetime.timezone.utc).timestamp()
        
        # seconds into day
        time_sec       = uaz_hi.timestamp[:,4]*(60*60) + \
                  uaz_hi.timestamp[:,5] * (60) + uaz_hi.timestamp[:,6] + \
                  uaz_hi.timestamp[:,7] / 1000
        
        # (time_sec is both 'time_offset' and 'time' in netCDF file)
        time_offset    = time_sec
        time           = time_sec

        # range dimension variables
        range_along_beam            = uaz_hi.range
        height_above_radar          = uaz_hi.range

        # These values are fixed for vertical mode. They will be time dependent in the wind mode...
        beam_azimuth_angle          = 0
        beam_elevation_angle        = 90
               
               
        # (time,range) dimension variables
        signal_to_noise_ratio   = moment_hi.snr
        reflectivity_factor     = moment_hi.zdB
        radial_velocity         = moment_hi.Vmean
        spectrum_width          = (2) * moment_hi.Vsig
        spectrum_skewness       = moment_hi.Vskew
        spectrum_kurtosis       = moment_hi.Vkurt
        spectrum_mean_noise_level   = 10 * np.log10(moment_hi.nos_mean)

        # do the actual netCDF file writting                              
        func_save_moments_in_netcdf(nc_output_hi_filename, 
                  bad_flag, 
                  source_str, 
                  site_id_str, 
                  facility_id_str, 
                  data_level_str, 
                  location_description_str, 
                  datastream_str, 
                  history_str, 
                  base_time, 
                  time_offset, 
                  time, 
                  lat, 
                  lon, 
                  alt, 
                  r_calib_radar_constant, 
                  num_time_domain_integrations, 
                  num_frequency_domain_integrations, 
                  num_fft_points, 
                  interpulse_period, 
                  operating_frequency_MHz, 
                  nyquist_velocity, 
                  pulse_width, 
                  range_resolution, 
                  num_code_bits, 
                  reference_noise_power, 
                  beam_azimuth_angle, 
                  beam_elevation_angle, 
                  range_along_beam, 
                  height_above_radar, 
                  signal_to_noise_ratio, 
                  reflectivity_factor, 
                  radial_velocity, 
                  spectrum_width, 
                  spectrum_skewness, 
                  spectrum_kurtosis, 
                  spectrum_mean_noise_level, 
                  date_str)

        #%% Cell: Make some early plots of the moments
        
        #%% Cell: Define whether there was any rain during this day
        
        # Determine whether this day has a rain event
        snr_min = 20

        # Set the min and max range
        range_min = 500  # meters
        range_max = 2500

        event_rain_flag, event_start_hour, event_end_hour = func_define_rain_event(
            snr_min, range_min, range_max, uaz_lo, moment_lo)
        
        #%% Cell: Plot 24-hour plot of Lo mode moments: SNR, Vmean, Vsig
        
        plot_range      = uaz_lo.range
        plot_drange     = plot_range[5] - plot_range[4]
        start_hour      = 0
        end_hour        = 24
        ylabels = ['0', '', '2', '', '4', '', '6', '', '8', '',
                   '10', '', '12', '', '14', '', '16', '', '18', '', '20']
        #ylabels = ['0', '', '', '', '4', '', '', '', '8', '',
        #           '', '', '12', '', '', '', '16', '', '', '', '20']
        # set flag if xlabels are [0 to 24]
        plot_xlabel_flag = 1    # flag to determine whether to plot xlabels
        xlabels = ['0', '', '', '3', '', '', '6', '', '', '9', '', '',
                   '12', '', '', '15', '', '', '18', '', '', '21', '', '', '24']
        
        input_top_data  = moment_lo.snr
        input_mid_data  = moment_lo.Vmean
        input_bot_data  = moment_lo.Vsig
        input_time      = uaz_lo.timestamp[:,8]
        
        # define the [min, step, max] values for colorbars
        top_color_values = np.array([-10,  10, 50]) 
        mid_color_values = np.array([ -2,   2, 10]) 
        bot_color_values = np.array([  0, 0.5,  3]) 
        
        top_color_ticks = ['-10', '', '10', '', '30', '', '50']
        mid_color_ticks = ['-2', '', ' 2', '', ' 6', '', '10']
        bot_color_ticks = ['0', '', '1', '', '2', '', '3']
        
        title_top_str   = f'{site_id_3CHAR_str} {facility_id_str}, RWP, Precip lo, $\Delta$r=%im, SNR [dB], {date_plot_str}' %plot_drange
        title_mid_str   = f'{site_id_3CHAR_str} {facility_id_str}, RWP, Precip lo, $\Delta$r=%im, Mean Velocity [m/s] (+ downward)' %plot_drange
        title_bot_str   = f'{site_id_3CHAR_str} {facility_id_str}, RWP, Precip lo, $\Delta$r=%im, Spectrum Breadth (Vsig) [m/s]' %plot_drange
        
        plot_filename = f'{image_day_directory}{image_lo_root_ZVW_filename}{date_str}'
        
        func_plot_time_ht_3panel(plot_range, plot_drange, start_hour, end_hour,
                         ylabels, plot_xlabel_flag, xlabels, input_top_data,
                         input_mid_data, input_bot_data, input_time,
                         top_color_values, mid_color_values, bot_color_values,
                         top_color_ticks, mid_color_ticks, bot_color_ticks,
                         title_top_str, title_mid_str, title_bot_str,
                         plot_filename)

        #%% Cell: Plot 24-hour plot of Lo mode moments: SNR, Skewness, Kurtosis
        
        plot_range      = uaz_lo.range
        plot_drange     = plot_range[5] - plot_range[4]
        start_hour      = 0
        end_hour        = 24
        ylabels = ['0', '', '2', '', '4', '', '6', '', '8', '',
                   '10', '', '12', '', '14', '', '16', '', '18', '', '20']
        #ylabels = ['0', '', '', '', '4', '', '', '', '8', '',
        #           '', '', '12', '', '', '', '16', '', '', '', '20']
        # set flag if xlabels are [0 to 24]
        plot_xlabel_flag = 1    # flag to determine whether to plot xlabels
        xlabels = ['0', '', '', '3', '', '', '6', '', '', '9', '', '',
                   '12', '', '', '15', '', '', '18', '', '', '21', '', '', '24']
        
        input_top_data  = moment_lo.snr
        input_mid_data  = moment_lo.Vskew
        input_bot_data  = moment_lo.Vkurt
        input_time      = uaz_lo.timestamp[:,8]
        
        # define the [min, step, max] values for colorbars
        top_color_values = np.array([-10,  10, 50]) 
        mid_color_values = np.array([ -2, 0.5,  2]) 
        bot_color_values = np.array([  1, 0.5,  5]) 
        
        top_color_ticks = ['-10', '', '10', '', '30', '', '50']
        mid_color_ticks = ['-2', '', ' -1', '', ' 0', '', '1', '', '2']
        bot_color_ticks = ['1', '', '2', '', '3', '', '4', '', '5']
        
        title_top_str   = f'{site_id_3CHAR_str} {facility_id_str}, RWP, Precip lo, $\Delta$r=%im, SNR [dB], {date_plot_str}' %plot_drange
        title_mid_str   = f'{site_id_3CHAR_str} {facility_id_str}, RWP, Precip lo, $\Delta$r=%im, Skewness [$m^3$/$s^3$]' %plot_drange
        title_bot_str   = f'{site_id_3CHAR_str} {facility_id_str}, RWP, Precip lo, $\Delta$r=%im, Kurtosis [$m^4$/$s^4$]' %plot_drange
        
        plot_filename = f'{image_day_directory}{image_lo_root_ZSK_filename}{date_str}'
        
        func_plot_time_ht_3panel(plot_range, plot_drange, start_hour, end_hour,
                         ylabels, plot_xlabel_flag, xlabels, input_top_data,
                         input_mid_data, input_bot_data, input_time,
                         top_color_values, mid_color_values, bot_color_values,
                         top_color_ticks, mid_color_ticks, bot_color_ticks,
                         title_top_str, title_mid_str, title_bot_str,
                         plot_filename)
                    
        #%% Cell: Plot 24-hour plot of Hi mode moments: SNR, Vmean, Vsig

        plot_range      = uaz_hi.range
        plot_drange     = plot_range[5] - plot_range[4]
        start_hour      = 0
        end_hour        = 24
        #ylabels = ['0', '', '2', '', '4', '', '6', '', '8', '',
        #           '10', '', '12', '', '14', '', '16', '', '18', '', '20']
        ylabels = ['0', '', '', '', '4', '', '', '', '8', '',
                   '', '', '12', '', '', '', '16', '', '', '', '20']
        # set flag if xlabels are [0 to 24]
        plot_xlabel_flag = 1    # flag to determine whether to plot xlabels
        xlabels = ['0', '', '', '3', '', '', '6', '', '', '9', '', '',
                   '12', '', '', '15', '', '', '18', '', '', '21', '', '', '24']
        
        input_top_data  = moment_hi.snr
        input_mid_data  = moment_hi.Vmean
        input_bot_data  = moment_hi.Vsig
        input_time      = uaz_hi.timestamp[:,8]
        
        # define the [min, step, max] values for colorbars
        top_color_values = np.array([-10,  10, 50]) 
        mid_color_values = np.array([ -2,   2, 10]) 
        bot_color_values = np.array([  0, 0.5,  3]) 
        
        top_color_ticks = ['-10', '', '10', '', '30', '', '50']
        mid_color_ticks = ['-2', '', ' 2', '', ' 6', '', '10']
        bot_color_ticks = ['0', '', '1', '', '2', '', '3']
        
        title_top_str   = f'{site_id_3CHAR_str} {facility_id_str}, RWP, Precip hi, $\Delta$r=%im, SNR [dB], {date_plot_str}' %plot_drange
        title_mid_str   = f'{site_id_3CHAR_str} {facility_id_str}, RWP, Precip hi, $\Delta$r=%im, Mean Velocity [m/s] (+ downward)' %plot_drange
        title_bot_str   = f'{site_id_3CHAR_str} {facility_id_str}, RWP, Precip hi, $\Delta$r=%im, Spectrum Breadth (Vsig) [m/s]' %plot_drange
        
        plot_filename = f'{image_day_directory}{image_hi_root_ZVW_filename}{date_str}'
        
        func_plot_time_ht_3panel(plot_range, plot_drange, start_hour, end_hour,
                         ylabels, plot_xlabel_flag, xlabels, input_top_data,
                         input_mid_data, input_bot_data, input_time,
                         top_color_values, mid_color_values, bot_color_values,
                         top_color_ticks, mid_color_ticks, bot_color_ticks,
                         title_top_str, title_mid_str, title_bot_str,
                         plot_filename)

        #%% Cell: Plot 24-hour plot of Hi mode moments: SNR, Skewness, Kurtosis
        
        plot_range      = uaz_hi.range
        plot_drange     = plot_range[5] - plot_range[4]
        start_hour      = 0
        end_hour        = 24
        #ylabels = ['0', '', '2', '', '4', '', '6', '', '8', '',
        #           '10', '', '12', '', '14', '', '16', '', '18', '', '20']
        ylabels = ['0', '', '', '', '4', '', '', '', '8', '',
                   '', '', '12', '', '', '', '16', '', '', '', '20']
        # set flag if xlabels are [0 to 24]
        plot_xlabel_flag = 1    # flag to determine whether to plot xlabels
        xlabels = ['0', '', '', '3', '', '', '6', '', '', '9', '', '',
                   '12', '', '', '15', '', '', '18', '', '', '21', '', '', '24']
        
        input_top_data  = moment_hi.snr
        input_mid_data  = moment_hi.Vskew
        input_bot_data  = moment_hi.Vkurt
        input_time      = uaz_hi.timestamp[:,8]
        
        # define the [min, step, max] values for colorbars
        top_color_values = np.array([-10,  10, 50]) 
        mid_color_values = np.array([ -2, 0.5,  2]) 
        bot_color_values = np.array([  1, 0.5,  5]) 
        
        top_color_ticks = ['-10', '', '10', '', '30', '', '50']
        mid_color_ticks = ['-2', '', ' -1', '', ' 0', '', '1', '', '2']
        bot_color_ticks = ['1', '', '2', '', '3', '', '4', '', '5']
        
        title_top_str   = f'{site_id_3CHAR_str} {facility_id_str}, RWP, Precip hi, $\Delta$r=%im, SNR [dB], {date_plot_str}' %plot_drange
        title_mid_str   = f'{site_id_3CHAR_str} {facility_id_str}, RWP, Precip hi, $\Delta$r=%im, Skewness [$m^3$/$s^3$]' %plot_drange
        title_bot_str   = f'{site_id_3CHAR_str} {facility_id_str}, RWP, Precip hi, $\Delta$r=%im, Kurtosis [$m^4$/$s^4$]' %plot_drange
        
        plot_filename = f'{image_day_directory}{image_hi_root_ZSK_filename}{date_str}'
        
        func_plot_time_ht_3panel(plot_range, plot_drange, start_hour, end_hour,
                         ylabels, plot_xlabel_flag, xlabels, input_top_data,
                         input_mid_data, input_bot_data, input_time,
                         top_color_values, mid_color_values, bot_color_values,
                         top_color_ticks, mid_color_ticks, bot_color_ticks,
                         title_top_str, title_mid_str, title_bot_str,
                         plot_filename)

        #%% Cell: Plot the Rain Event (if there was rain)

        ### Plot a few hours of data: lo mode: Z, Vmean, Vsig ###

        if(event_rain_flag):

            #%% Cell: Plot 24-hour plot of Lo mode moments: SNR, Vmean, Vsig
            
            plot_range      = uaz_lo.range
            plot_drange     = plot_range[5] - plot_range[4]
            start_hour      = event_start_hour
            end_hour        = event_end_hour
            ylabels = ['0', '', '2', '', '4', '', '6', '', '8', '',
                       '10', '', '12', '', '14', '', '16', '', '18', '', '20']
            #ylabels = ['0', '', '', '', '4', '', '', '', '8', '',
            #           '', '', '12', '', '', '', '16', '', '', '', '20']
            # set flag if xlabels are [0 to 24]
            plot_xlabel_flag = 0    # flag to determine whether to plot xlabels
            xlabels = ['0', '', '', '3', '', '', '6', '', '', '9', '', '',
                       '12', '', '', '15', '', '', '18', '', '', '21', '', '', '24']
            
            input_top_data  = moment_lo.snr
            input_mid_data  = moment_lo.Vmean
            input_bot_data  = moment_lo.Vsig
            input_time      = uaz_lo.timestamp[:,8]
            
            # define the [min, step, max] values for colorbars
            top_color_values = np.array([-10,  10, 50]) 
            mid_color_values = np.array([ -2,   2, 10]) 
            bot_color_values = np.array([  0, 0.5,  3]) 
            
            top_color_ticks = ['-10', '', '10', '', '30', '', '50']
            mid_color_ticks = ['-2', '', ' 2', '', ' 6', '', '10']
            bot_color_ticks = ['0', '', '1', '', '2', '', '3']
            
            title_top_str   = f'{site_id_3CHAR_str} {facility_id_str}, RWP, Precip lo, $\Delta$r=%im, SNR [dB], {date_plot_str}' %plot_drange
            title_mid_str   = f'{site_id_3CHAR_str} {facility_id_str}, RWP, Precip lo, $\Delta$r=%im, Mean Velocity [m/s] (+ downward)' %plot_drange
            title_bot_str   = f'{site_id_3CHAR_str} {facility_id_str}, RWP, Precip lo, $\Delta$r=%im, Spectrum Breadth (Vsig) [m/s]' %plot_drange
            
            plot_filename = f'{image_event_directory}{image_event_lo_root_ZVW_filename}{date_str}'
            
            func_plot_time_ht_3panel(plot_range, plot_drange, start_hour, end_hour,
                             ylabels, plot_xlabel_flag, xlabels, input_top_data,
                             input_mid_data, input_bot_data, input_time,
                             top_color_values, mid_color_values, bot_color_values,
                             top_color_ticks, mid_color_ticks, bot_color_ticks,
                             title_top_str, title_mid_str, title_bot_str,
                             plot_filename)

            #%% Cell: Plot 24-hour plot of Lo mode moments: SNR, Skewness, Kurtosis
        
            plot_range     = uaz_lo.range
            plot_drange = plot_range[5] - plot_range[4]
            start_hour      = event_start_hour
            end_hour        = event_end_hour
            ylabels = ['0', '', '2', '', '4', '', '6', '', '8', '',
                       '10', '', '12', '', '14', '', '16', '', '18', '', '20']
            #ylabels = ['0', '', '', '', '4', '', '', '', '8', '',
            #           '', '', '12', '', '', '', '16', '', '', '', '20']
            # set flag if xlabels are [0 to 24]
            plot_xlabel_flag = 0    # flag to determine whether to plot xlabels
            xlabels = ['0', '', '', '3', '', '', '6', '', '', '9', '', '',
                       '12', '', '', '15', '', '', '18', '', '', '21', '', '', '24']
            
            input_top_data  = moment_lo.snr
            input_mid_data  = moment_lo.Vskew
            input_bot_data  = moment_lo.Vkurt
            input_time      = uaz_lo.timestamp[:,8]
            
            # define the [min, step, max] values for colorbars
            top_color_values = np.array([-10,  10, 50]) 
            mid_color_values = np.array([ -1, 0.25,  1]) 
            bot_color_values = np.array([  1, 0.5,  5]) 
            
            top_color_ticks = ['-10', '', '10', '', '30', '', '50']
            mid_color_ticks = ['-1.0', '', '0.5', '', '0.0', '', '0.5', '', '1.0']
            bot_color_ticks = ['1', '', '2', '', '3', '', '4', '', '5']
            
            title_top_str   = f'{site_id_3CHAR_str} {facility_id_str}, RWP, Precip lo, $\Delta$r=%im, SNR [dB], {date_plot_str}' %plot_drange
            title_mid_str   = f'{site_id_3CHAR_str} {facility_id_str}, RWP, Precip lo, $\Delta$r=%im, Skewness [$m^3$/$s^3$]' %plot_drange
            title_bot_str   = f'{site_id_3CHAR_str} {facility_id_str}, RWP, Precip lo, $\Delta$r=%im, Kurtosis [$m^4$/$s^4$]' %plot_drange
            
            plot_filename = f'{image_event_directory}{image_event_lo_root_ZSK_filename}{date_str}'
            
            func_plot_time_ht_3panel(plot_range, plot_drange, start_hour, end_hour,
                             ylabels, plot_xlabel_flag, xlabels, input_top_data,
                             input_mid_data, input_bot_data, input_time,
                             top_color_values, mid_color_values, bot_color_values,
                             top_color_ticks, mid_color_ticks, bot_color_ticks,
                             title_top_str, title_mid_str, title_bot_str,
                             plot_filename)
                    
            #%% Cell: Plot 24-hour plot of Hi mode moments: SNR, Vmean, Vsig
    
            plot_range      = uaz_hi.range
            plot_drange     = plot_range[5] - plot_range[4]
            start_hour      = event_start_hour
            end_hour        = event_end_hour
            #ylabels = ['0', '', '2', '', '4', '', '6', '', '8', '',
            #           '10', '', '12', '', '14', '', '16', '', '18', '', '20']
            ylabels = ['0', '', '', '', '4', '', '', '', '8', '',
                       '', '', '12', '', '', '', '16', '', '', '', '20']
            # set flag if xlabels are [0 to 24]
            plot_xlabel_flag = 0    # flag to determine whether to plot xlabels
            xlabels = ['0', '', '', '3', '', '', '6', '', '', '9', '', '',
                       '12', '', '', '15', '', '', '18', '', '', '21', '', '', '24']
            
            input_top_data  = moment_hi.snr
            input_mid_data  = moment_hi.Vmean
            input_bot_data  = moment_hi.Vsig
            input_time      = uaz_hi.timestamp[:,8]
            
            # define the [min, step, max] values for colorbars
            top_color_values = np.array([-10,  10, 50]) 
            mid_color_values = np.array([ -2,   2, 10]) 
            bot_color_values = np.array([  0, 0.5,  3]) 
            
            top_color_ticks = ['-10', '', '10', '', '30', '', '50']
            mid_color_ticks = ['-2', '', ' 2', '', ' 6', '', '10']
            bot_color_ticks = ['0', '', '1', '', '2', '', '3']
            
            title_top_str   = f'{site_id_3CHAR_str} {facility_id_str}, RWP, Precip hi, $\Delta$r=%im, SNR [dB], {date_plot_str}' %plot_drange
            title_mid_str   = f'{site_id_3CHAR_str} {facility_id_str}, RWP, Precip hi, $\Delta$r=%im, Mean Velocity [m/s] (+ downward)'%plot_drange
            title_bot_str   = f'{site_id_3CHAR_str} {facility_id_str}, RWP, Precip hi, $\Delta$r=%im, Spectrum Breadth (Vsig) [m/s]'%plot_drange
            
            plot_filename = f'{image_event_directory}{image_event_hi_root_ZVW_filename}{date_str}'
            
            func_plot_time_ht_3panel(plot_range, plot_drange, start_hour, end_hour,
                             ylabels, plot_xlabel_flag, xlabels, input_top_data,
                             input_mid_data, input_bot_data, input_time,
                             top_color_values, mid_color_values, bot_color_values,
                             top_color_ticks, mid_color_ticks, bot_color_ticks,
                             title_top_str, title_mid_str, title_bot_str,
                             plot_filename)
    
            #%% Cell: Plot 24-hour plot of Hi mode moments: SNR, Skewness, Kurtosis
            
            plot_range      = uaz_hi.range
            plot_drange     = plot_range[5] - plot_range[4]
            start_hour      = event_start_hour
            end_hour        = event_end_hour
            #ylabels = ['0', '', '2', '', '4', '', '6', '', '8', '',
            #           '10', '', '12', '', '14', '', '16', '', '18', '', '20']
            ylabels = ['0', '', '', '', '4', '', '', '', '8', '',
                       '', '', '12', '', '', '', '16', '', '', '', '20']
            # set flag if xlabels are [0 to 24]
            plot_xlabel_flag = 0    # flag to determine whether to plot xlabels
            xlabels = ['0', '', '', '3', '', '', '6', '', '', '9', '', '',
                       '12', '', '', '15', '', '', '18', '', '', '21', '', '', '24']
            
            input_top_data  = moment_hi.snr
            input_mid_data  = moment_hi.Vskew
            input_bot_data  = moment_hi.Vkurt
            input_time      = uaz_hi.timestamp[:,8]
            
            # define the [min, step, max] values for colorbars
            top_color_values = np.array([-10,  10, 50]) 
            mid_color_values = np.array([ -1, 0.25,  1]) 
            bot_color_values = np.array([  1, 0.5,  5]) 
            
            top_color_ticks = ['-10', '', '10', '', '30', '', '50']
            mid_color_ticks = ['-1.0', '', '0.5', '', '0.0', '', '0.5', '', '1.0']
            bot_color_ticks = ['1', '', '2', '', '3', '', '4', '', '5']
            
            title_top_str   = f'{site_id_3CHAR_str} {facility_id_str}, RWP, Precip hi, $\Delta$r=%im, SNR [dB], {date_plot_str}' %plot_drange
            title_mid_str   = f'{site_id_3CHAR_str} {facility_id_str}, RWP, Precip hi, $\Delta$r=%im, Skewness [$m^3$/$s^3$]' %plot_drange
            title_bot_str   = f'{site_id_3CHAR_str} {facility_id_str}, RWP, Precip hi, $\Delta$r=%im, Kurtosis [$m^4$/$s^4$]' %plot_drange
            
            plot_filename = f'{image_event_directory}{image_event_hi_root_ZSK_filename}{date_str}'
            
            func_plot_time_ht_3panel(plot_range, plot_drange, start_hour, end_hour,
                             ylabels, plot_xlabel_flag, xlabels, input_top_data,
                             input_mid_data, input_bot_data, input_time,
                             top_color_values, mid_color_values, bot_color_values,
                             top_color_ticks, mid_color_ticks, bot_color_ticks,
                             title_top_str, title_mid_str, title_bot_str,
                             plot_filename)
            
        # end if(event_rain_flag):
    
    # end if(process_spc_flag):

    #%% Print an 'end-of-day processing' message
    
    print(" ")
    print(
        f'Finished processing moments for day {date_plot_str}')

# end for single_date in func_daterange(start_date, end_date):
    
#%% Cell: End of cells

