import time
import numpy as np


def func_convert_netCDF_time(base_time, time_offset):

    # combine all time offsets with base time (in sec from epoch)
    time_combine_secs_epoch = np.array(base_time + time_offset)

    # organize time as array (year month day hour min sec dayofweek dayofyear daylightsaving)
    time_vector = np.array(time.gmtime(base_time))

    # pad zeros
    time_profile = np.zeros(
        (len(time_combine_secs_epoch), len(time_vector)))

    # loop through all time_combine_secs_epoch and convert to time vectors
    for idx, time_sec in np.ndenumerate(time_combine_secs_epoch):
        time_profile[idx, :] = np.array(time.gmtime(time_sec), ndmin=2)

    # get time as fraction and millisecs
    time_hour_fraction = time_profile[:, 3] + \
        time_profile[:, 4]/60 + time_profile[:, 5]/(60*60)

    time_millisec = time_profile[:, 5] * 0  # No millisec info

    # UzarTimeStamp format [year month day dayofyear hour min sec millisec timefraction]
    uzar_time_stamp = ([
        time_profile[:, 0],  # year
        time_profile[:, 1],  # month
        time_profile[:, 2],  # day
        time_profile[:, 7],  # day of year
        time_profile[:, 3],  # hour
        time_profile[:, 4],  # min
        time_profile[:, 5],  # sec
        time_millisec,
        time_hour_fraction]
    )
    uzar_time_stamp_t = np.transpose(uzar_time_stamp)
    return uzar_time_stamp_t
