#%% Cell: Define the imports and needed functions

import numpy as np

#%% Cell: Describe the function

# [Vmean_corrected, index_lt_corrected, index_rt_corrected] = ...
#                   func_time_dealias_on_profile(Vmean_input, index_lt_input, index_rt_input, ...
#                   V_Nyquist, Npts, num_prior_profiles, V_dif_threshold, num_range_gate_threshold);
                  

# This function performs a temporal dealiasing on a time-ht velocity matrix.

# Aliasing is examined over the whole profile.
# If enough range gates have aliasing, then perform a dealiasing on those
# range gates.

# Inputs:
# Vmean_input        - Velocities (bad values filled with NaNs)
# index_lt_input     - left integration limit
# index_rt_input     - right integration limit
# V_Nyquist          - Nyquist velocity
# Npts               - number of points in velocity spectrum
# num_prior_profiles - number of prior profiles used to estimate
#                      Vmean_prior.
# V_dif_threshold    - velocity threshold used to indicate aliasing. Should
#                       be a function of V_Nyquist, e.g. 1.5*V_Nyquist.
# num_range_gate_threshold - number of aliased range gates in profile
#                       needed before performing dealiasing procedure on 
#                       each range gate.

# Outputs:
# Vmean_corrected    - velocities after correction
# index_lt_corrected - left integration limit after correction
# index_rt_corrected - right integration limit after correction

# updated: 11-November-2021
# ========================================================================

#%% Cell: Start the function

def func_time_dealias_on_profile(Vmean_input, index_lt_input, index_rt_input,
                                 V_Nyquist, Npts, num_prior_profiles, 
                                 V_dif_threshold, num_range_gate_threshold):
    
    #%% Define the output variables

    m_profiles  = np.size(Vmean_input, 0)
    nhts        = np.size(Vmean_input, 1)

    Vmean_corrected     = 1 * Vmean_input
    index_lt_corrected  = 1 * index_lt_input
    index_rt_corrected  = 1 * index_rt_input
    

    #%% Process with time moving forward, from left to right through the matrix

    # Move from left to right through the data.
    # If the velocity difference is greater than V_dif_threshold,
    # then adjust the velocity by 2*V_Nyquist.

    for r in range(m_profiles):
   
        # make sure there are prior profiles
        if(r >= num_prior_profiles):
        
            # reset the profile arrays
            profile_need_to_fix_flag   = np.zeros(nhts)
            profile_dif                = np.zeros(nhts)
            profile_test               = np.zeros(nhts)
            profile_prior              = np.zeros(nhts)
   
            # process each range gate of this profile
            for c in range(nhts):
      
                # set the test value
                Vmean_test   = Vmean_corrected[r,c]
      
                # Procede if Vmean_test is valid
                if(~np.isnan(Vmean_test)):
   
                    # get the previous samples, determine median or set to zero
                    start_index    = r - num_prior_profiles
                    end_index      = r - 1
         
                    if((start_index - end_index) < 1):
                        Vmean_prior = Vmean_corrected[start_index,c]
                    
                    else:
                        array_prior = Vmean_corrected[np.arange(start_index,end_index+1,1),c]
                        f              = ~np.isnan(array_prior)
                        # three different cases: 
                        # based on number of valid values in array
                        if(np.sum(f) == 0):
                            Vmean_prior    = 0
                        # end if(sum(f) == 0)
                        if(np.sum(f) == 1):
                            Vmean_prior    = array_prior[f]
                        # end if(sum(f) == 0)            
                        if(np.sum(f) > 1):
                            Vmean_prior    = np.median(array_prior[f])
                        # end if(sum(f) > 1)

                        # if neighbor pixel is NaN, then do not do dealiasing
                        if(np.isnan(f[-1])):
                            Vmean_prior    = 0
                        # end if(isnan(f(end)))
                   
                    # end if((start_index - end_index) < 1):
                
                    # calculate the difference between the test and prior.
                    Vmean_dif  = Vmean_test - Vmean_prior
         
                    profile_dif[c]    = Vmean_dif
                    profile_test[c]   = Vmean_test
                    profile_prior[c]  = Vmean_prior
         
                    if(Vmean_dif > V_dif_threshold):
                        profile_need_to_fix_flag[c]   = 1
                    # end if(Vmean_dif > V_dif_threshold)
         
                    if(Vmean_dif < (-1) * V_dif_threshold):
                        profile_need_to_fix_flag[c]   = -1;
                    # end if(Vmean_dif < (-1) * V_dif_threshold)
         
                # end if(~isnan(Vmean_test))
      
            # end for c loop
   
            #%% Count the number of consecutive positive offset samples
            
            running_positive_cnt    = 0
            current_cnt             = 0
   
            for c in range(nhts):
      
                if(profile_need_to_fix_flag[c] > 0.5):
                    current_cnt = current_cnt + 1
         
                    if(current_cnt > running_positive_cnt):
                        running_positive_cnt = current_cnt
                    # end if(current_cnt < running_cnt)
                else:
                    current_cnt = 0
                # end if(pofile_need_to_fix_flag(c) > 0.5)
            
            # end for c loop

            #%% Count the number of consecutive negative offset samples
            
            running_negative_cnt    = 0
            current_cnt             = 0
   
            for c in range(nhts):
      
                if(profile_need_to_fix_flag[c] < -0.5):
                    current_cnt = current_cnt + 1
         
                    if(current_cnt > running_negative_cnt):
                        running_negative_cnt = current_cnt
                    # end if(current_cnt < running_cnt)
                else:
                    current_cnt = 0
                # end if(pofile_need_to_fix_flag(c) > 0.5)
        
            # end for c loop
   
            #%% Process this profile if there are enough aliased range gates    

            # Check the positive offsets
            if( running_positive_cnt > num_range_gate_threshold):
      
                # if here, then remove 2*V_Nyquist to the needed ranges gates
                for c in range(nhts):
         
                    # is there aliasing in this range gate?
                    if(profile_need_to_fix_flag[c] > 0.5):
            
                        # correct the velocities
                        Vmean_corrected[r,c]      = Vmean_corrected[r,c] - 2*V_Nyquist
            
                        # correct the indices
                        index_lt_corrected[r,c]   = index_lt_corrected[r,c] - Npts
                        index_rt_corrected[r,c]   = index_rt_corrected[r,c] - Npts
            
                    # end if(profile_need_to_fix_flag(c) > 0.5)
         
                # end for c loop
      
            # end if( running_positive_cnt > num_range_gate_threshold)
   
            # Check the negative offsets
            if( running_negative_cnt > num_range_gate_threshold):
      
                # if here, then add 2*V_Nyquist to the needed ranges gates
                for c in range(nhts):
         
                    # is there aliasing in this range gate?
                    if(profile_need_to_fix_flag[c] < -0.5):
            
                        # correct the velocities
                        Vmean_corrected[r,c]      = Vmean_corrected[r,c] + 2*V_Nyquist
            
                        # correct the indices
                        index_lt_corrected[r,c]   = index_lt_corrected[r,c] + Npts
                        index_rt_corrected[r,c]   = index_rt_corrected[r,c] + Npts
            
                    # end if( running_negative_cnt > num_range_gate_threshold)
         
                # end for c loop
      
            # end if( sum(f_negative) > num_range_gate_threshold)
   
        # end if(r >= num_prior_profiles):
    # end for r loop

    #%% Process with time moving backward, from right to left through the matrix

    # Move from right to left through the data.
    # If the velocity difference is greater than V_dif_threshold,
    # then adjust the velocity by 2*V_Nyquist.

    # Logic of this section...
    # Get the previous samples (which to the right of the current sample), 
    # and determine the median of these value or set the median to zero.
    # Since moving from right to left, the previous samples are 
    #   start_index  = r + 1; # first sample on right
    #   end_index    = r + num_prior_profiles; # last sample on right

    for r in range(m_profiles-1,-1,-1):

        # make sure there are post-profiles
        if(r < m_profiles - num_prior_profiles - 1):
        
            # reset the profile arrays
            profile_need_to_fix_flag    = np.zeros(nhts)
            profile_dif                 = np.zeros(nhts)
            profile_test                = np.zeros(nhts)
            profile_prior               = np.zeros(nhts)
   
            # process each range gate of this profile
            for c in range(nhts):
      
                # set the test value
                Vmean_test   = Vmean_corrected[r,c]

                # Procede if Vmean_test is valid
                if(~np.isnan(Vmean_test)):
   
                    # get the previous samples, determine median or set to zero
                    start_index    = r + 1
                    end_index      = r + num_prior_profiles
         
                    if((start_index - end_index) < 1):
                        Vmean_prior = Vmean_corrected[start_index,c]
                    
                    else:
                        array_prior = Vmean_corrected[np.arange(start_index,end_index+1,1),c]
                        f              = ~np.isnan(array_prior)
                        # three different cases: 
                        # based on number of valid values in array
                        if(np.sum(f) == 0):
                            Vmean_prior    = 0
                        # end if(sum(f) == 0)
                        if(np.sum(f) == 1):
                            Vmean_prior    = array_prior[f]
                        # end if(sum(f) == 0)            
                        if(np.sum(f) > 1):
                            Vmean_prior    = np.median(array_prior[f])
                        # end if(sum(f) > 1)

                        # if neighbor pixel is NaN, then do not do dealiasing
                        if(np.isnan(f[-1])):
                            Vmean_prior    = 0
                        # end if(isnan(f(end)))
                   
                    # end if((start_index - end_index) < 1):
                
                    # calculate the difference between the test and prior.
                    Vmean_dif  = Vmean_test - Vmean_prior
         
                    profile_dif[c]    = Vmean_dif
                    profile_test[c]   = Vmean_test
                    profile_prior[c]  = Vmean_prior
         
                    if(Vmean_dif > V_dif_threshold):
                        profile_need_to_fix_flag[c]   = 1
                    # end if(Vmean_dif > V_dif_threshold)
         
                    if(Vmean_dif < (-1) * V_dif_threshold):
                        profile_need_to_fix_flag[c]   = -1;
                    # end if(Vmean_dif < (-1) * V_dif_threshold)
         
                # end if(~isnan(Vmean_test))
      
            # end for c loop
   
            #%% Count the number of consecutive positive offset samples
            
            running_positive_cnt    = 0
            current_cnt             = 0
   
            for c in range(nhts):
      
                if(profile_need_to_fix_flag[c] > 0.5):
                    current_cnt = current_cnt + 1
         
                    if(current_cnt > running_positive_cnt):
                        running_positive_cnt = current_cnt
                    # end if(current_cnt < running_cnt)
                else:
                    current_cnt = 0
                # end if(pofile_need_to_fix_flag(c) > 0.5)
            
            # end for c loop

            #%% Count the number of consecutive negative offset samples
            
            running_negative_cnt    = 0
            current_cnt             = 0
   
            for c in range(nhts):
      
                if(profile_need_to_fix_flag[c] < -0.5):
                    current_cnt = current_cnt + 1
         
                    if(current_cnt > running_negative_cnt):
                        running_negative_cnt = current_cnt
                    # end if(current_cnt < running_cnt)
                else:
                    current_cnt = 0
                # end if(pofile_need_to_fix_flag(c) > 0.5)
        
            # end for c loop
   
            #%% Process this profile if there are enough aliased range gates    

            # Check the positive offsets
            if( running_positive_cnt > num_range_gate_threshold):
      
                # if here, then remove 2*V_Nyquist to the needed ranges gates
                for c in range(nhts):
         
                    # is there aliasing in this range gate?
                    if(profile_need_to_fix_flag[c] > 0.5):
            
                        # correct the velocities
                        Vmean_corrected[r,c]      = Vmean_corrected[r,c] - 2*V_Nyquist
            
                        # correct the indices
                        index_lt_corrected[r,c]   = index_lt_corrected[r,c] - Npts
                        index_rt_corrected[r,c]   = index_rt_corrected[r,c] - Npts
            
                    # end if(profile_need_to_fix_flag(c) > 0.5)
         
                # end for c loop
      
            # end if( running_positive_cnt > num_range_gate_threshold)
   
            # Check the negative offsets
            if( running_negative_cnt > num_range_gate_threshold):
      
                # if here, then add 2*V_Nyquist to the needed ranges gates
                for c in range(nhts):
         
                    # is there aliasing in this range gate?
                    if(profile_need_to_fix_flag[c] < -0.5):
            
                        # correct the velocities
                        Vmean_corrected[r,c]      = Vmean_corrected[r,c] + 2*V_Nyquist
            
                        # correct the indices
                        index_lt_corrected[r,c]   = index_lt_corrected[r,c] + Npts
                        index_rt_corrected[r,c]   = index_rt_corrected[r,c] + Npts
            
                    # end if( running_negative_cnt > num_range_gate_threshold)
         
                # end for c loop
      
            # end if( sum(f_negative) > num_range_gate_threshold)
   
        # end if(r >= num_prior_profiles):
            
    # end for r loop

    #%% Cell: Return spc_output
        
    return Vmean_corrected, index_lt_corrected, index_rt_corrected

# end def func_time_dealias_on_profile