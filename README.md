# High-resolution (1 km) Köppen-Geiger maps for 1901–2099 based on constrained CMIP6 projections

This repository contains the code for the paper titled [High-resolution (1 km) Köppen-Geiger maps for 1901–2099 based on constrained CMIP6 projections](https://doi.org/10.1038/s41597-023-02549-6), published in Scientific Data.

## Part 1: Climate Model Analysis

Scripts for processing CMIP6 climate model data and evaluating model sensitivities.

### `climate_model_step1_data_to_npz.py`

This script processes the CMIP6 climate model data. It loops over various scenarios, models, and ensemble members, loading the data from raw NetCDF files. The data is then resampled, and various statistics such as mean, minimum, and maximum are computed for both global and land surface time series. The processed data is saved in NPZ format.

### `climate_model_step2_compute_stats.py`

This script loads the HadCRUT observed air temperature ensemble data, estimates internal climate variability using CMIP6 data, and compares both simulated and observed historical air temperature trends. It computes the Equilibrium Climate Sensitivity (ECS) and Transient Climate Response (TCR) for the CMIP6 models. The script also screens models based on historical trends and TCR. 

### `climate_model_step3_generate_figs_tables.py`

This script generates a scatterplot illustrating the historical trend versus TCR for the paper. It also produces a LaTeX table detailing the number of ensemble members, trends, and sensitivities. Additionally, it creates box plots that show warming trends for the land surface under different scenarios, considering both all models and the screened subset of models.

## Part 2: Climatologies and Köppen-Geiger maps

These scripts pertain to the generation of high-resolution precipitation and air temperature climatologies. They also produce Köppen-Geiger maps for both historical periods and future scenarios.

### `climatologies_step1_historical.py`

This script loads the CRU air temperature and GPCC precipitation datasets. It then iterates over each combination of station-based air temperature and precipitation climatic datasets, generating adjusted high-resolution climatologies for each historical period using change offsets and factors derived from the CRU and GPCC data. Lastly, the script calculates final climatologies, their associated uncertainty estimates, and produces Köppen-Geiger maps for the historical periods.

### `climatologies_step2_future.py`

This script loads the climate model data (both historical data and future projections) and computes the average for all ensemble members for each model. It then produces high-resolution future climatologies, maps of ensemble mean and standard deviation, and the final Köppen-Geiger maps for the future periods.

### `climatologies_step3_resample_and_package.py`

For each high-resolution Köppen-Geiger map or climatology in NetCDF format, this script produces a low-resolution, upscaled version. It then converts the Köppen-Geiger maps to geoTIFF format, including a color map. Finally, the script creates six zip files to facilitate data distribution.

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