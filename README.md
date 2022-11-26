# RWP_Python_moments
This repository contains Python codes to estimate spectrum moments from Radar Wind Profiler (RWP) Doppler velcoity power spectra. 
Details of the processing steps are described in Williams et al. "Calibrating radar wind profiler reflectivity factor using surface disdrometer observations", submitted to Atmospheric Measurement Techniques, 02-Dec-2022.

# Directory Structure
On your system, create five directories:
   1) images_precip_day
   2) image_precip_event
   3) nc_precip_crw0
   4) py_files
   5) rwp_precip_spec

# Input Soure Code
Place the Python code from this repository into the directory: py_files

# Input Source Data
Go to the US Department of Energy (DOE) Atmospheric Radiation Measurement (ARM) program Data Discovery (https://adc.arm.gov/discovery/#/) and download the precipitation mode RWP spectra.  

The link (http://dx.doi.org/10.5439/1025129) will take you to Data Discovery and will select the correct data type. 

Place the downloaded input soure files into the directory: rwp_precip_spec

# Output netCDF files
The generated netCDF files will be placed into the directory: nc_precip_crw0

# Output Daily Images
The code will generate daily images and will place those images into the directory: images_precip_day

# Output Event Images
The code will generate event images (if there is a precipitation event) and will place those images into the directory: images_precip_event

# Running the Source Code
1) The user should only need to modify the "main" routine. All other functions are labeled "func" and should not need to be modified.

2) The user will need to verify the paths defined for your system.

3) Change the start and end dates in the "main" routine to include the date range of the source data files. 

