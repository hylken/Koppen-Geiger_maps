# High-resolution (1 km) Köppen-Geiger maps for 1901–2099 based on constrained CMIP6 projections

This repository contains the code for the paper titled "High-resolution (1 km) Köppen-Geiger maps for 1901–2099 based on constrained CMIP6 projections," published in Scientific Data.

## Part 1: Climate Model Analysis

Scripts for processing CMIP6 model data and evaluating model sensitivities.

### `climate_model_step1_data_to_npz.py`

This script processes the climate model data. It loops over various scenarios, models, and ensemble members, loading the data from NetCDF files. The data is then resampled, and various statistics such as mean, minimum, and maximum are computed for both global and land surface time series. The processed data is saved in NPZ format, and optional verification figures are generated to visualize the time series and spatial distribution of the variables.

### `climate_model_step2_compute_stats.py`

This script loads the HadCRUT observed air temperature ensemble data, estimates internal climate variability using CMIP6 data, and compares both simulated and observed historical air temperature trends. It computes the Equilibrium Climate Sensitivity (ECS) and Transient Climate Response (TCR) for the CMIP6 models. The script also screens models based on historical trends and TCR. If a TCR value is not available for a model, only the historical trend is used. Additionally, it creates a .txt file that summarizes the likely historical trend, ECS, and TCR ranges for the paper and generates a dataframe detailing projected mean changes in temperature and precipitation for the land surface, both with and without model screening.

### `climate_model_step3_generate_figs_tables.py`

This script generates a scatterplot illustrating the historical trend versus TCR for the paper. It also produces a .tex table detailing the number of ensemble members, trends, and sensitivities. The script also prints the mean global warming for each scenario. Additionally, it creates box plots that display warming trends for the land surface under different scenarios, considering both all models and a the screened subset of models. Finally, it produces figures depicting projected global mean changes in temperature and precipitation, along with their associated uncertainties.

## Part 2: Climatologies and Köppen-Geiger maps

These scripts pertain to the generation of high-resolution precipitation and air temperature climatologies. They also produce Köppen-Geiger maps for both historical periods and future scenarios.

### `climatologies_step1_historical.py`

This script loads the CRU air temperature and GPCC precipitation datasets. It then iterates over each combination of station-based air temperature and precipitation climatic datasets, generating adjusted high-resolution climatologies for each historical period using change offsets and factors derived from the CRU and GPCC data. Lastly, based on the ensemble of adjusted station-based climatologies, the script calculates final climatologies, their associated uncertainty estimates, and produces Köppen-Geiger maps for the historical periods.

### `climatologies_step2_future.py`

This script processes climate model data to generate high-resolution future climatologies. It first loads the climate model data and averages it over ensemble members. The script then loops over the variables and ensemble members, creating a dictionary to store the data. For each variable and ensemble member, it checks if both historical and projection data exist. If not, it skips processing. The script then loads the historical and projection data from npz files and stores it in the dictionary. It checks the temporal completeness of the data, and if it is incomplete, it sets the data to NaN. The script then computes the average over all ensemble members for each variable.

Next, the script generates high-resolution future climatologies by looping over the variables and months, loading a high-resolution historic reference climatology for each variable and month. It produces a change map to adjust the climatology to the target period and month based on the monthly climate model data. The code then saves the temporally adjusted, high-resolution future climatology to an nc file. After processing all models, the script computes the ensemble mean and standard deviation and the Köppen-Geiger classification maps and uncertainty.

### `climatologies_step3_resample_and_package.py`

For each high-resolution Köppen-Geiger map or climatology in netCDF format, this script produces a low-resolution, upscaled version. It then converts the Köppen-Geiger maps to geoTIFF format, including a color map. Finally, the script creates six zip files to facilitate data distribution.

### `climatologies_step4_validation.py`

This script calculates the areas covered by major Köppen-Geiger classes and the transitions between periods. It also computes the Köppen-Geiger classifications using station data. Subsequently, the script evaluates the classification accuracy for each historical period, considering both the 30 Köppen-Geiger classes and the major classes, by loading the global Köppen-Geiger map for each respective period.

### `climatologies_step5_generate_figs_tables.py`

This script generates figures of Köppen-Geiger maps for various periods and scenarios. Additionally, it creates a LaTeX table detailing the classification accuracy and constructs node and link lists for the Sankey diagrams. The script then produces figures of the Sankey diagrams for all periods and scenarios.

# Instructions

Clone the repository:
```
git clone https://github.com/hylken/Koppen-Geiger_maps
cd Koppen-Geiger_maps
```
Produce a configuration file with the correct paths and folders based on the provided template (`config_kaust.py`).

Create the environment and install the packages as follows:
```
conda update conda
conda config --set solver libmamba
conda create -n Koppen-Geiger_maps -c conda-forge scipy pandas numpy netcdf4 matplotlib pymatreader statsmodels adjusttext seaborn scikit-image basemap rasterio pyshp plotly python-kaleido
```