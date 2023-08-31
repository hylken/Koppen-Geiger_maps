
# Title

![2071_2099_ssp245_World_koppen_geiger_0p01](https://github.com/hylken/Koppen-Geiger_maps/assets/122258032/41fedb32-30a5-4d43-a167-9a60e6063b31)


This repository contains XXX for paper XXX.

## Climate Model Analysis

Scripts related to processing the CMIP6 model data and assessing the sensitivities of the models.

### `climate_model_step1_data_to_npz.py`

This script processes the climate model data. It loops over various scenarios, models, and ensemble members, loading the data from NetCDF files. The data is then resampled, and various statistics such as mean, minimum, and maximum are computed for both global and land surface time series. The processed data is saved in NPZ format, and optional verification figures are generated to visualize the time series and spatial distribution of the variables.

### `climate_model_step2_compute_stats.py`

This script loads the HadCRUT observed air temperature ensemble data, estimates internal climate variability using CMIP6 data, and compares both simulated and observed historical air temperature trends. It computes the Equilibrium Climate Sensitivity (ECS) and Transient Climate Response (TCR) for the CMIP6 models. The script also screens models based on historical trends and TCR. If a TCR value is not available for a model, only the historical trend is used. Additionally, it creates a .txt file that summarizes the likely historical trend, ECS, and TCR ranges for the paper and generates a dataframe detailing projected mean changes in temperature and precipitation for the land surface, both with and without model screening.

### `climate_model_step3_generate_figs_tables.py`

This script generates a scatterplot illustrating the historical trend versus TCR for the paper. It also produces a .tex table detailing the number of ensemble members, trends, and sensitivities. The script also prints the mean global warming for each scenario. Additionally, it creates box plots that display warming trends for the land surface under different scenarios, considering both all models and a the screened subset of models. Finally, it produces figures depicting projected global mean changes in temperature and precipitation, along with their associated uncertainties.

## Climatologies and Köppen-Geiger maps

Scripts related to generating the final high-resolution precipitation and air temperature climatologies and the Köppen-Geiger maps for the historical and future periods and future scenarios.

### `climatologies_step1_historical.py`

This script loads the CRU air temperature and GPCC precipitation datasets. It then iterates over each combination of station-based air temperature and precipitation climatic datasets, generating adjusted high-resolution climatologies for each historical period using change offsets and factors derived from the CRU and GPCC data. Lastly, based on the ensemble of adjusted station-based climatologies, the script calculates final climatologies, their associated uncertainty estimates, and produces Köppen-Geiger maps for the historical periods.

### `climatologies_step2_future.py`
The script loads climate model data and averages it over ensemble members. It then then loops over the variables and ensemble members, creating a dictionary to store the data. For each variable and ensemble member, it checks if both historical and projection data exist. If not, it skips processing. It then loads the historical and projection data from .npz files and stores it in the dictionary. It checks the temporal completeness of the data and if it is incomplete, it sets the data to NaN. It then computes the average over all ensemble members for each variable.

It then generates high-resolution future climatologies by looping over the variables and months, loading a high-resolution historic reference climatology for each variable and month. It produces a change map to adjust the climatology to the target period and month based on the monthly climate model data. The code then saves the temporally adjusted, high-resolution future climatology to an nc file. After processing all models, the code computes the ensemble mean and standard deviation and computes Köppen-Geiger classification maps and uncertainty.

### `climatologies_step3_resample_and_package.py`

The script resamples high-resolution netCDF files to create low-resolution, upscaled versions and saves them as netCDF files. It loops over all files in the climatologies directory and checks if the file is an ensemble mean or standard deviation file or a Köppen-Geiger map file. For each file, it loops over the upscale map sizes specified in the configuration file and creates a new file name with the appropriate suffix. If the file is an ensemble mean or standard deviation file, the code resamples the precipitation and air temperature climatologies for each month and variable and saves the resampled data to a new netCDF file.

The code then converts the Köppen-Geiger maps from netCDF files to geoTIFF files. It first creates a Köppen-Geiger colormap for the geoTIFF using the koppen_table data. The code then loops over all files in the climatologies directory and checks if the file is a Köppen-Geiger map file. For each file, it creates a new file name with the .tif extension and loads the kg_class data from the netCDF file. The code then saves the data to a geoTIFF file using the tools.write_to_geotiff function and the previously created colormap.

### `climatologies_step4_validation.py`

### `climatologies_step5_generate_figs_tables.py`

## License

(Provide information on the licensing of the code)

## Usage

## Dependencies

# Overview

# Data

# System requirements


# Instructions

Clone the repository:
```
git clone https://github.com/hylken/Koppen-Geiger_maps
cd Koppen-Geiger_maps
```
Produce a configuration file with the correct paths and folders based on the provided template (`config.ini`).

Create the environment and install the packages as follows:
```
conda update conda
conda install -n base -c conda-forge mamba
mamba create -n Koppen-Geiger_maps -c conda-forge scipy pandas numpy netcdf4 matplotlib pymatreader statsmodels adjusttext seaborn scikit-image basemap rasterio pyshp plotly python-kaleido
pip install pyshp
```
