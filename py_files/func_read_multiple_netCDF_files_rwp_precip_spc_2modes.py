#%% Cell: Define the imports and needed functions

import numpy as np
import glob

from func_read_rwp_precip_scp_2modes_netCDF import \
    func_read_rwp_precip_spc_2modes_netCDF

#%% Cell: Start the function

def func_read_multiple_netCDF_files_rwp_precip_spc_2modes(
        input_netCDF_directory, input_root_filename, single_date, 
        input_suffix, rwp_freq_MHz):

    #%% Cell: Define the initial values
    
    # Reset the Valid data flag
    process_spc_flag = 0

    # Get date string in file format YearMonthDay 20190401
    date_str = single_date.strftime('%Y%m%d')
    date_plot = single_date.strftime('%Y-%m-%d')

    # Determine whether there are any files for this day.
    #files_to_process = glob.glob(
    #    f'{input_netCDF_directory}{input_root_filename}{date_str}{input_suffix}')
    files_to_process = glob.glob(
        f'{input_netCDF_directory}\\{input_root_filename}{date_str}{input_suffix}')

    # How many files are there for this date
    num_files_to_process = len(files_to_process)

    # % There may be one or there could be multiple files for this day.
    # % The logic to account for multiple files is as follows:
    # % 1. Load data from the first file.
    # % 2. If more than one file, then read and concatenate data.
    # % Code assumes operating parameters do not change during the day.

    # % The operating parameters need to be the same for each file.
    # % Also, the number of range gates and the number of spectral
    # % bins need to be the same for all files. The first file
    # % defines operating parameters for the day.
    # % Any file with different operating parameters in a single day
    # % must be processed separately.

    # Desplay how many files will be processed
    if(num_files_to_process < 0.5):
        print(
            f"There are not any raw data files for day: {date_plot}...")
        
        # Set the output variables to dumby values
        uaz_lo  = 0
        uaz_hi  = 0
        
    # end if(num_files_to_process M< 0.5)
        
    #%% Cell: Process the first file for this day
    
    # Just process the first file for this day
    if(num_files_to_process > .5):

        # Set the file_loop to the first record
        file_loop = 0

        # Get the next filename
        input_netCDF_filename = files_to_process[file_loop]

        # Validate file
        try:
            f = open(input_netCDF_filename)
            f.close()
        except IOError:
            print(
                f"Reference file {input_netCDF_filename} not accessible")

        # Read data
        print(
            f'...Reading file #{file_loop + 1} out of {num_files_to_process} files for this day...')

        # Read the UAZR spectra
        print(f'...Loading netCDF file: {input_netCDF_filename}')

        uaz_lo, uaz_hi = func_read_rwp_precip_spc_2modes_netCDF(
            input_netCDF_filename, rwp_freq_MHz)

        # for item, value in uaz_lo.__dict__.items():
        #     print(item, ':', np.shape(value))

        # The objects in uaz_lo and uaz_hi are:
        # shape of () = scalar

        # Ncoh : ()
        # ipp : ()
        # Nspc : ()
        # Npts : (128,)
        # pulse_length : ()
        # Vd : (128,)
        # spc_pop : (1664, 150, 128)
        # timestamp : (1664, 9)
        # range : (150,)
        # lat : ()
        # lon : ()
        # alt : ()

        process_spc_flag = 1

        # At this point:
        # if process_spc_flag == 1, then there is spc data to process

        # Check to see if there are multiple files for this day.
        # If no, then continue with spc processing
        # If yes, then read each file and concatenate data.

    # end if(num_files_to_process > .5):
        
    #%% Cell: Process multiple files, if present
    
    # process all other files for this day, if they are present
    if(num_files_to_process > 1.5):

        # Define keeper matricies
        keep_uaz_hi_pop_spc = uaz_hi.spc_pop
        keep_uaz_hi_timestamp = uaz_hi.timestamp
        keep_uaz_lo_pop_spc = uaz_lo.spc_pop
        keep_uaz_lo_timestamp = uaz_lo.timestamp

        # keep track of the operating parameters and only keep new
        # observations that have the same parameters. Days with
        # different parameters will need to be processed separately.
        # Logic is the that first file of the day defines the
        # operating parameters and all other files are not kept.

        keep_uaz_hi_Ncoh = uaz_hi.Ncoh
        keep_uaz_hi_Npts = uaz_hi.Npts
        keep_uaz_hi_Nspc = uaz_hi.Nspc
        keep_uaz_hi_ipp = uaz_hi.ipp
        keep_uaz_hi_code_bits = uaz_hi.code_bits
        keep_uaz_hi_pulse_length = uaz_hi.pulse_length
        keep_uaz_hi_Vd = uaz_hi.Vd
        keep_uaz_hi_range = uaz_hi.range

        keep_uaz_lo_Ncoh = uaz_lo.Ncoh
        keep_uaz_lo_Npts = uaz_lo.Npts
        keep_uaz_lo_Nspc = uaz_lo.Nspc
        keep_uaz_lo_ipp = uaz_lo.ipp
        keep_uaz_lo_code_bits = uaz_lo.code_bits
        keep_uaz_lo_pulse_length = uaz_lo.pulse_length
        keep_uaz_lo_Vd = uaz_lo.Vd
        keep_uaz_lo_range = uaz_lo.range
        keep_uaz_site_LatLonAlt = [uaz_lo.lat, uaz_lo.lon, uaz_lo.alt]

        # Read file number 2 and higher, and concatenate data
        for file_loop in range(1, num_files_to_process):

            # Get the next filename
            input_netCDF_filename = files_to_process[file_loop]

            # Validate file
            try:
                f = open(input_netCDF_filename)
                f.close()
            except IOError:
                print(
                    f"Reference file {input_netCDF_filename} not accessible")

            # Read data
            print(
                f'Processing file: {file_loop + 1} out of {num_files_to_process}...')

            uaz_lo, uaz_hi = func_read_rwp_precip_spc_2modes_netCDF(
                input_netCDF_filename, rwp_freq_MHz)

            # Concatenate the data...hi mode
            # are the operating parameter the same as before?
            if((keep_uaz_hi_Ncoh == uaz_hi.Ncoh) & (keep_uaz_hi_Npts == uaz_hi.Npts)
                & (keep_uaz_hi_Nspc == uaz_hi.Nspc) & (keep_uaz_hi_ipp == uaz_hi.ipp)
                    & (keep_uaz_hi_code_bits == uaz_hi.code_bits)
                    & (keep_uaz_hi_pulse_length == uaz_hi.pulse_length)):

                # Get sizes of old and current matrices
                [m_old, n_old, p_old] = np.shape(keep_uaz_hi_pop_spc)
                [m_current, n_current, p_current] = np.shape(
                    uaz_hi.spc_pop)

                # Require the number of range gates and number of
                # spectral points to be the same
                if((n_old == n_current) & (p_old == p_current)):

                    # Concatenate the data
                    keep_uaz_hi_pop_spc = np.concatenate(
                        (keep_uaz_hi_pop_spc, uaz_hi.spc_pop), axis=0)

                    keep_uaz_hi_timestamp = np.concatenate(
                        (keep_uaz_hi_timestamp, uaz_hi.timestamp), axis=0)
                # end if((n_old == n_current) & (p_old == p_current)):
                        
            # end if((keep_uaz_hi_Ncoh == uaz_hi.Ncoh) & (keep_uaz_hi_Npts == uaz_hi.Npts)
                
            # Concatenate the data...lo mode
            # are the operating parameter the same as before?
            if((keep_uaz_lo_Ncoh == uaz_lo.Ncoh) & (keep_uaz_lo_Npts == uaz_lo.Npts)
                & (keep_uaz_lo_Nspc == uaz_lo.Nspc) & (keep_uaz_lo_ipp == uaz_lo.ipp)
                    & (keep_uaz_lo_code_bits == uaz_lo.code_bits)
                    & (keep_uaz_lo_pulse_length == uaz_lo.pulse_length)):

                # Get sizes of old and current matrices
                [m_old, n_old, p_old] = np.shape(keep_uaz_lo_pop_spc)
                [m_current, n_current, p_current] = np.shape(
                    uaz_lo.spc_pop)

                # Require the number of range gates and number of
                # spectral points to be the same
                if((n_old == n_current) & (p_old == p_current)):

                    # Concatenate the data
                    keep_uaz_lo_pop_spc = np.concatenate(
                        (keep_uaz_lo_pop_spc, uaz_lo.spc_pop), axis=0)

                    keep_uaz_lo_timestamp = np.concatenate(
                        (keep_uaz_lo_timestamp, uaz_lo.timestamp), axis=0)
                # end if((n_old == n_current) & (p_old == p_current)):
                    
            # end if((keep_uaz_lo_Ncoh == uaz_lo.Ncoh) & (keep_uaz_lo_Npts == uaz_lo.Npts)
        
        # end for file_loop in range(1, num_files_to_process):
        
        # Change variable names back to original names
        
        # change variables back to original name
        uaz_hi.spc_pop = keep_uaz_hi_pop_spc
        uaz_hi.timestamp = keep_uaz_hi_timestamp
        uaz_lo.spc_pop = keep_uaz_lo_pop_spc
        uaz_lo.timestamp = keep_uaz_lo_timestamp

        # set operating parameter back to first file's value
        uaz_hi.Ncoh = keep_uaz_hi_Ncoh
        uaz_hi.Npts = keep_uaz_hi_Npts
        uaz_hi.Nspc = keep_uaz_hi_Nspc
        uaz_hi.ipp = keep_uaz_hi_ipp
        uaz_hi.code_bits = keep_uaz_hi_code_bits
        uaz_hi.pulse_length = keep_uaz_hi_pulse_length
        uaz_hi.Vd = keep_uaz_hi_Vd
        uaz_hi.range = keep_uaz_hi_range
        [uaz_hi.lat, uaz_hi.lon, uaz_hi.alt] = keep_uaz_site_LatLonAlt

        uaz_lo.Ncoh = keep_uaz_lo_Ncoh
        uaz_lo.Npts = keep_uaz_lo_Npts
        uaz_lo.Nspc = keep_uaz_lo_Nspc
        uaz_lo.ipp = keep_uaz_lo_ipp
        uaz_lo.code_bits = keep_uaz_lo_code_bits
        uaz_lo.pulse_length = keep_uaz_lo_pulse_length
        uaz_lo.Vd = keep_uaz_lo_Vd
        uaz_lo.range = keep_uaz_lo_range
        [uaz_lo.lat, uaz_lo.lon, uaz_lo.alt] = keep_uaz_site_LatLonAlt

    # end if(num_files_to_process > 1.5):
    
    #%% Return values from function call
    
    return uaz_lo, uaz_hi, process_spc_flag
