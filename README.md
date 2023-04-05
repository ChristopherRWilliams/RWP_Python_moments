# RWP_Python_moments
This repository contains Python codes to estimate spectrum moments from Radar Wind Profiler (RWP) Doppler velocity power spectra. 
Details of the processing steps are described in Williams et al., "Calibrating radar wind profiler reflectivity factor using surface disdrometer observations", published in Atmospheric Measurement Techniques, https://doi.org/10.5194/egusphere-2022-1405, 2023.

# Directory Structure
On your system, create five directories:
   1) images_precip_day
   2) image_precip_event
   3) nc_precip_crw0
   4) py_files
   5) rwp_precip_spec

# Input Source Code
Place the Python code from this repository into the directory: py_files

# Input Source Data
Go to the US Department of Energy (DOE) Atmospheric Radiation Measurement (ARM) program Data Discovery (https://adc.arm.gov/discovery/#/) and download at least one precipitation mode RWP spectra data file from the Southern Great Plains (SGP) Central Facility (C1). An example data file has the format: sgp915rwpprecipspecC1.a0.yyyyMMdd.hhmmss.cdf, where yyyyMMdd.hhmmss is the date and time at the beginning of the data file. 

The link (http://dx.doi.org/10.5439/1025129) will take you to Data Discovery and will select the correct data type. 

Place the downloaded input source files into the directory: rwp_precip_spec

# Output netCDF files
The generated netCDF files will be placed into the directory: nc_precip_crw0

# Output Daily Images
The code will generate daily images and will place those images into the directory: images_precip_day

# Output Event Images
The code will generate event images (if there is precipitation during a day) and will place those images into the directory: images_precip_event

# Running the Source Code
1) The user should only need to modify the "main" routine. All other functions are labeled "func" and should not need to be modified.

2) The user will need to change the paths so that they work for your system.

3) Change the start and end dates in the "main" routine to include the date range of the source data files that have been downloaded. 

