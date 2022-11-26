# % program: func_load_all_netCDF_variables
# % updated: 06-June-2021

# % This routine loads all netCDF variables into a dictionary of numpy arrays:

# % Function needs the variable named: input_netCDF_filename
import netCDF4
import numpy as np


def func_load_all_netCDF_variables(filename):
    my_data = {}
    cdf_data = netCDF4.Dataset(filename)
    variables = list(cdf_data.variables)
    for name in variables:
        my_data[name] = np.array(cdf_data[name][:])
    cdf_data.close()

    return my_data
