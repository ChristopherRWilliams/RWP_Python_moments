#%% Write moments in netCDF format

#%% Define imports

#from netCDF4 import Dataset    # Note: python is case-sensitive!
import netCDF4
import numpy as np

def func_save_moments_in_netcdf(output_netCDF_filename, 
                  input_bad_flag, 
                  input_source_str, 
                  input_site_id_str, 
                  input_facility_id_str, 
                  input_data_level_str, 
                  input_location_description_str, 
                  input_datastream_str, 
                  input_history_str, 
                  input_base_time, 
                  input_time_offset, 
                  input_time, 
                  input_lat, 
                  input_lon, 
                  input_alt, 
                  input_r_calib_radar_constant, 
                  input_num_time_domain_integrations, 
                  input_num_frequency_domain_integrations, 
                  input_num_fft_points, 
                  input_interpulse_period, 
                  input_operating_frequency_MHz, 
                  input_nyquist_velocity, 
                  input_pulse_width, 
                  input_range_resolution, 
                  input_num_code_bits, 
                  input_reference_noise_power, 
                  input_beam_azimuth_angle, 
                  input_beam_elevation_angle, 
                  input_range_along_beam, 
                  input_height_above_radar, 
                  input_signal_to_noise_ratio, 
                  input_reflectivity_factor, 
                  input_radial_velocity, 
                  input_spectrum_width, 
                  input_spectrum_skewness, 
                  input_spectrum_kurtosis, 
                  input_spectrum_mean_noise_level, 
                  input_date_str):

    #%% Cell: Set all NaNs to bad_flag value
    
    # Determine size of moment matrices
    m_profiles  = np.size(input_signal_to_noise_ratio, 0)
    nhts        = np.size(input_signal_to_noise_ratio, 1)

    ## Process each Profile ##
    for r in range(m_profiles):

        for c in range(nhts):

            # process signal_to_noise_ratio
            if(np.isnan(input_signal_to_noise_ratio[r,c])):
                input_signal_to_noise_ratio[r,c] = input_bad_flag

            # process reflectivity_factor
            if(np.isnan(input_reflectivity_factor[r,c])):
                input_reflectivity_factor[r,c] = input_bad_flag
                
            # process radial_velocity
            if(np.isnan(input_radial_velocity[r,c])):
                input_radial_velocity[r,c] = input_bad_flag

            # process spectrum_width
            if(np.isnan(input_spectrum_width[r,c])):
                input_spectrum_width[r,c] = input_bad_flag

            # process skewness
            if(np.isnan(input_spectrum_skewness[r,c])):
                input_spectrum_skewness[r,c] = input_bad_flag

            # process kurtosis
            if(np.isnan(input_spectrum_kurtosis[r,c])):
                input_spectrum_kurtosis[r,c] = input_bad_flag

            # process mean_noise (in profile)
            if(np.isnan(input_spectrum_mean_noise_level[r,c])):
                input_spectrum_mean_noise_level[r,c] = input_bad_flag

        # end for c loop

    # end for prof_num loop

    #%% Create new empty netCDF4 file

    rootgrp = netCDF4.Dataset(
        output_netCDF_filename, mode='w', format='NETCDF4')
            
    #%% Define Global Attributes
    
    rootgrp.input_source                = input_source_str
    rootgrp.site_id                     = input_site_id_str
    rootgrp.facility_id                 = input_facility_id_str
    rootgrp.data_level                  = input_data_level_str
    rootgrp.location_description        = input_location_description_str
    rootgrp.datastream                  = input_datastream_str
    rootgrp.num_radar_pulses_per_dwell  = "(num_time_domain_integrations)(num_fft_points)(num_frequency_domain_integrations)"
    rootgrp.dwell_time                  = "(num_radar_pulses_per_dwell)(interpulse_period)"
    rootgrp.history                     = input_history_str    

    #%% Create dimensions

    #time_dim        = rootgrp.createDimension('time', len(input_time))
    #gate_dim        = rootgrp.createDimension('gate', len(input_range_along_beam)
    #single_dim      = rootgrp.createDimension('scalar', 1)
    rootgrp.createDimension('time', len(input_time_offset))
    rootgrp.createDimension('gate', len(input_range_along_beam))
    rootgrp.createDimension('scalar', 1)

    #%% Create variables 

    # Define the time dependent variables
    
    # base_time (seconds since 1970-01-01 0:00:00 0:00)
    base_time                                   = rootgrp.createVariable('base_time', np.single, ('scalar'))
    base_time.units                             = 'seconds since 1970-1-1 0:00:00 0:00'
    base_time.long_nane                         = 'Base time in Epoch (01-January-1970)'
    base_time.ancillary_variables               = 'time_offset'

    # time_offset (seconds since base_time)
    time_offset                                 = rootgrp.createVariable('time_offset', np.single, ('time'))
    time_offset.units                           = 'seconds since base_time' 
    time_offset.long_nane                       = 'Time offset from base_time. Time is beginning of dwell.'
    time_offset.ancillary_variables             = 'base_time'

    # time  (seconds since midnight)
    time                                        = rootgrp.createVariable('time', np.single, ('time'))
    time.units                                  = 'seconds since midnight' 
    time.long_nane                              = 'Time since midnight. Time is beginning of dwell.'
    time.standard_name                          = 'time'
    time.calendar                               = 'gregorian'

    # Define the scalers
    lat                                         = rootgrp.createVariable('lat', np.single, ('scalar'))
    lat.units                                   = 'degrees North Latitude'
    lat.long_nane                               = 'Site latitude'
    
    lon                                         = rootgrp.createVariable('lon', np.single, ('scalar'))
    lon.units                                   = 'degrees East Longitude'
    lon.long_nane                               = 'Site longitude'

    alt                                         = rootgrp.createVariable('alt', np.single, ('scalar'))
    alt.units                                   = 'm'
    alt.long_nane                               = 'Site Altitude above sea level'
    
    r_calib_radar_constant                      = rootgrp.createVariable('r_calib_radar_constant', np.single, ('scalar'))
    r_calib_radar_constant.units                = 'dB'
    r_calib_radar_constant.long_nane            = 'Calibration constant relative to reference of 0 dB'

    num_time_domain_integrations                = rootgrp.createVariable('num_time_domain_integrations', np.single, ('scalar'))
    num_time_domain_integrations.units          = 'integer'
    num_time_domain_integrations.long_nane      = 'Number of integrated time domain samples'
    
    num_frequency_domain_integrations           = rootgrp.createVariable('num_frequency_domain_integrations', np.single, ('scalar'))
    num_frequency_domain_integrations.units     = 'integer'
    num_frequency_domain_integrations.long_nane = 'Number of integrated frequency domain spectra.'

    num_fft_points                              = rootgrp.createVariable('num_fft_points', np.single, ('scalar'))
    num_fft_points.units                        = 'integer'
    num_fft_points.long_nane                    = 'Number of FFT points in frequency domain spectra.'

    interpulse_period                           = rootgrp.createVariable('interpulse_period', np.single, ('scalar'))
    interpulse_period.units                     = 'seconds'
    interpulse_period.long_nane                 = 'Interpulse Period, time between transmitted pulses'

    operating_frequency_MHz                     = rootgrp.createVariable('operating_frequency_MHz', np.single, ('scalar'))
    operating_frequency_MHz.units               = 'MHz'
    operating_frequency_MHz.long_nane           = 'Radar operating frequency'

    nyquist_velocity                            = rootgrp.createVariable('nyquist_velocity', np.single, ('scalar'))
    nyquist_velocity.units                      = 'm/s'
    nyquist_velocity.long_nane                  = 'Unambiguous radial velocity, also known as Nyquist velocity = ((299.79)/operating_frequency_MHz)/(4*num_time_domain_integrations*interpulse_period)'

    pulse_width                                 = rootgrp.createVariable('pulse_width', np.single, ('scalar'))
    pulse_width.units                           = 'seconds'
    pulse_width.long_nane                       = 'Duration of transmitted pulse'
    
    range_resolution                            = rootgrp.createVariable('range_resolution', np.single, ('scalar'))
    range_resolution.units                      = 'meters'
    range_resolution.long_nane                  = 'If num_code_bits < 1, Range resolution = (3x10^8)*(pulse_width)/2. If num_code_bits >= 1, Range resolution = (3x10^8)*(pulse_width)/(2*num_code_bits)'

    num_code_bits                               = rootgrp.createVariable('num_code_bits', np.single, ('scalar'))
    num_code_bits.units                         = 'integer'
    num_code_bits.long_nane                     = 'Number of phase bits in transmitted coded pulse. Value of zero indicates no pulse phase coding.'

    reference_noise_power                       = rootgrp.createVariable('reference_noise_power', np.single, ('scalar'))
    reference_noise_power.units                 = 'dB'
    reference_noise_power.long_nane             = 'Reference noise power used in calibration = 10log(refernce_spectrum_mean_noise_level [linear units] * num_fft_points)'
    
    beam_azimuth_angle                          = rootgrp.createVariable('beam_azimuth_angle', np.single, ('scalar'))
    beam_azimuth_angle.units                    = 'degrees'
    beam_azimuth_angle.long_nane                = 'Beam azimuth angle from North (North = 0, East = 90)'

    beam_elevation_angle                        = rootgrp.createVariable('beam_elevation_angle', np.single, ('scalar'))
    beam_elevation_angle.units                  = 'degrees'
    beam_elevation_angle.long_nane              = 'Beam elevation angle from horizontal (horizontal = 0, vertical = 90)'

    # Define the range dependent variables
    # range along radial direction
    range_along_beam                            = rootgrp.createVariable('range_along_beam', np.single, ('gate'))
    range_along_beam.units                      = 'meters'
    range_along_beam.positive                   = 'up'
    range_along_beam.long_nane                  = 'Range in the radial direction from the radar to center of the range gate'

    height_above_radar                          = rootgrp.createVariable('height_above_radar', np.single, ('gate'))
    height_above_radar.units                    = 'meters'
    height_above_radar.positive                 = 'up'
    height_above_radar.long_nane                = 'Height above the radar without regard to horizontal displacement'

    # Define the range-time dependent variabiles
    
    # uaz_snr
    signal_to_noise_ratio                       = rootgrp.createVariable('signal_to_noise_ratio', np.single, ('time', 'gate'))
    signal_to_noise_ratio.units                 = 'dB'
    signal_to_noise_ratio.missing_value         = input_bad_flag
    signal_to_noise_ratio.long_nane             = 'Doppler velocity power spectrum Signal-to-Noise Ratio with mean noise adjustment, also known as 0th moment'
    
    # uaz_zdb
    reflectivity_factor                         = rootgrp.createVariable('reflectivity_factor', np.single, ('time', 'gate'))
    reflectivity_factor.units                   = 'dBZ'
    reflectivity_factor.missing_value           = input_bad_flag
    reflectivity_factor.long_nane               = 'Radar Reflectivity factor: zdb [dBZ] = snr [dB] + 20log10(range_along_beam) [dB] + r_calib_radar_constant [dB]'
    
    # uaz_Vmean
    radial_velocity                             = rootgrp.createVariable('radial_velocity', np.single, ('time', 'gate'))
    radial_velocity.units                       = 'm/s'
    radial_velocity.missing_value               = input_bad_flag
    radial_velocity.positive                    = 'approaching the radar'
    radial_velocity.long_nane                   = 'Doppler velocity power spectrum mean radial velocity, also known as 1st moment, positive values approaching the radar'
    
    # uaz_Vsig*2
    spectrum_width                              = rootgrp.createVariable('spectrum_width', np.single, ('time', 'gate'))
    spectrum_width.units                        = 'm/s'
    spectrum_width.missing_value                = input_bad_flag
    spectrum_width.positive                     = 'approaching the radar'
    spectrum_width.long_nane                    = 'Doppler velocity power spectrum width, defined as 2*sqrt(spectrum variance)'
    
    # uaz_Vskew
    spectrum_skewness                           = rootgrp.createVariable('spectrum_skewness', np.single, ('time', 'gate'))
    spectrum_skewness.units                     = '(m/s)^3'
    spectrum_skewness.missing_value             = input_bad_flag
    spectrum_skewness.positive                  = 'approaching the radar'
    spectrum_skewness.long_nane                 = 'Doppler velocity power spectrum skewness, also known as 3rd moment, positive values approaching the radar'
    
    # uaz_Vkurt
    spectrum_kurtosis                           = rootgrp.createVariable('spectrum_kurtosis', np.single, ('time', 'gate'))
    spectrum_kurtosis.units                     = '(m/s)^4'
    spectrum_kurtosis.missing_value             = input_bad_flag
    spectrum_kurtosis.positive                  = 'approaching the radar'
    spectrum_kurtosis.long_nane                 = 'Doppler velocity power spectrum kurtosis, also known as 4th moment, positive values approaching the radar'

    # uaz_nos_mean
    spectrum_mean_noise_level                   = rootgrp.createVariable('spectrum_mean_noise_level', np.single, ('time', 'gate'))
    spectrum_mean_noise_level.units             = 'dB'
    spectrum_mean_noise_level.missing_value     = input_bad_flag    
    spectrum_mean_noise_level.long_nane         = 'Doppler velocity power spectrum mean noise level determined from Hildebrand & Sekhon (1974) method'
    
    #%% Write data to netCDF file
    
    # Time dependent variables
    
    base_time[:]                            = input_base_time
    time_offset[:]                          = input_time_offset
    time[:]                                 = input_time
           
    # Scalers
    lat[:]                                  = input_lat
    lon[:]                                  = input_lon
    alt[:]                                  = input_alt

    r_calib_radar_constant[:]               = input_r_calib_radar_constant

    num_time_domain_integrations[:]         = input_num_time_domain_integrations
    num_frequency_domain_integrations[:]    = input_num_frequency_domain_integrations
    num_fft_points[:]                       = input_num_fft_points
    interpulse_period[:]                    = input_interpulse_period
    
    operating_frequency_MHz[:]              = input_operating_frequency_MHz
    nyquist_velocity[:]                     = input_nyquist_velocity
    
    pulse_width[:]                          = input_pulse_width
    range_resolution[:]                     = input_range_resolution
    num_code_bits[:]                        = input_num_code_bits
    reference_noise_power[:]                = input_reference_noise_power

    beam_azimuth_angle[:]                   = input_beam_azimuth_angle
    beam_elevation_angle[:]                 = input_beam_elevation_angle
    
    # Range dependent variables
    range_along_beam[:]                     = input_range_along_beam        
    height_above_radar[:]                   = input_height_above_radar
    
    # Range-Time dependent variabiles
    signal_to_noise_ratio[:,:]              = input_signal_to_noise_ratio
    reflectivity_factor[:,:]                = input_reflectivity_factor
    radial_velocity[:,:]                    = input_radial_velocity
    spectrum_width[:,:]                     = input_spectrum_width
    spectrum_skewness[:,:]                  = input_spectrum_skewness
    spectrum_kurtosis[:,:]                  = input_spectrum_kurtosis
    spectrum_mean_noise_level[:,:]          = input_spectrum_mean_noise_level                      
    
    #%% Close the netCDF file
    
    rootgrp.close()
