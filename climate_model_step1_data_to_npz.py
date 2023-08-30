#!/usr/bin/env python  
# -*- coding: utf-8 -*-

__author__ = "Hylke E. Beck"
__email__ = "hylke.beck@gmail.com"
__date__ = "November 2022"

import os
import sys
import pdb
import time
import glob
import random
import pandas as pd
import numpy as np
import tools
from datetime import datetime
from netCDF4 import Dataset 
import matplotlib.pyplot as plt
import scipy.io
import gc
from skimage.transform import resize_local_mean

# For some models the resolution changes during the simulation 
# For example, see the following two files:
# pr_Amon_EC-Earth3-LR_piControl_r1i1p1f1_gr_199101-199101.nc
# pr_Amon_EC-Earth3-LR_piControl_r1i1p1f1_gr_221901-221912.nc
# This issue should be addressed at some point

def main():


    #==============================================================================
    #   Settings
    #==============================================================================

    config = tools.load_config(sys.argv[1])

    scenarios = ['1pctCO2','abrupt-4xCO2','piControl','historical','ssp119','ssp126','ssp245','ssp370','ssp434','ssp460','ssp585']
    random.shuffle(scenarios)
    variables = ['tas','pr','rsdt','rsut','rlut']

    generate_figures = True

    # Load land surface mask
    mfile = scipy.io.loadmat(os.path.join('assets','GSHHS_oceanmask.mat'))
    mask_full = mfile['DATA'].astype(bool)
    
    
    #==============================================================================
    #   Loop over scenarios, models, and members and load data and generate 
    #   figures for verication.
    #==============================================================================

    if os.path.isdir(os.path.join(config['folder_out'],'climate_model_data','figures'))==False:
        os.makedirs(os.path.join(config['folder_out'],'climate_model_data','figures'))

    # Loop over variables
    for vv in np.arange(len(variables)):
        variable = variables[vv]
        
        # Loop over scenarios        
        for ss in np.arange(len(scenarios)):
            scenario = scenarios[ss]
            
            if ((variable==('rsdt')) | (variable==('rsut')) | (variable==('rlut'))) & ((scenario!='abrupt-4xCO2') & (scenario!='piControl')):
                continue
            
            # Generate list of models
            files = glob.glob(os.path.join(config['folder_dataraw'],'CMIP6', scenario, variable+'_Amon*.nc'))
            models = sorted(np.unique([os.path.basename(x).split('_')[2] for x in files]).tolist())
            try:
                models.remove('historical')
            except:
                pass
            random.shuffle(models)
            
            # Loop over models
            for mm in np.arange(len(models)):
                model = models[mm]
                
                #if (model!='EC-Earth3') | (scenario!='1pctCO2') | (variable!='tas'):
                #    continue
                    
                # Generate list of ensemble members
                files = glob.glob(os.path.join(config['folder_dataraw'],'CMIP6',scenario,variable+'_Amon_'+model+'_*.nc'))
                members = sorted(np.unique([os.path.basename(x).split('_')[4] for x in files]).tolist())
                random.shuffle(members)
                
                # Loop over ensemble members
                for ee in np.arange(len(members)):
                    member = members[ee]
                    
                    # Check if already processed
                    if os.path.isfile(os.path.join(config['folder_out'],'climate_model_data',scenario+'_'+model+'_'+member+'_'+variable+'.npz')):
                        continue
                    
                    # List of files to be loaded
                    files = glob.glob(os.path.join(config['folder_dataraw'],'CMIP6', scenario,variable+'_Amon_'+model+'_'+scenario+'_'+member+'_*.nc'))
                    files = sorted(files)
                    if len(files)==0:
                        continue
                        
                    # Make date array
                    date_start = datetime.strptime(os.path.basename(files[0]).split('_')[6][:6], '%Y%m').replace(month=1,day=1) 
                    date_end = datetime.strptime(os.path.basename(files[-1]).split('_')[6][7:13], '%Y%m').replace(month=12,day=31)
                    DatesMon = pd.period_range(start=date_start, end=date_end, freq='M')
                    Years = np.unique(DatesMon.year)                    
                    
                    print('-------------------------------------------------------------------------------')
                    print('Processing '+variable+' '+scenario+' '+model+' '+member) 
                    t0 = time.time()
                    
                    # Load first file to determine model grid dimensions
                    try:
                        print('Loading '+files[0])
                        dset = Dataset(files[0])
                        ncdata = np.array(dset[variable][:])
                        dset.close()
                    except:
                        print('Unable to read '+files[-1]+', skipping')
                        continue
                    
                    # Initialize data array
                    max_memory_usage = 64 # GB
                    try:    
                        nelem = ncdata.shape[1]*ncdata.shape[2]*len(DatesMon)
                    except:
                        print('Skipping because just two dimensions')
                        continue
                    mem_usage = nelem*4/10**9 # GB
                    if mem_usage>max_memory_usage:
                        print('Skipping because too much data for memory')
                        continue
                    data = np.zeros((ncdata.shape[1],ncdata.shape[2],len(DatesMon)),dtype=np.single)*np.NaN
                    
                    # Resample mask to compute land surface mean time series
                    mask_small = resize_local_mean(mask_full.astype(np.single),(ncdata.shape[1],ncdata.shape[2]))>0.5
                    
                    # Loop over files
                    for ff in np.arange(len(files)):
                        try:
                            dset = Dataset(files[ff])
                            ncdata = np.array(dset[variable][:])
                            dset.close()
                        except:
                            print('Unable to read '+files[ff])
                            continue
                  
                        date_start = datetime.strptime(os.path.basename(files[ff]).split('_')[6][:6], '%Y%m')
                        date_end = datetime.strptime(os.path.basename(files[ff]).split('_')[6][7:13], '%Y%m')
                        time_arr = pd.period_range(start=date_start, end=date_end, freq='M')

                        # Check if time and data array are same length
                        if len(time_arr)!=ncdata.shape[0]:
                            print(files[ff]+' netCDF data and time fields inconsistent, skipping')
                            continue
                        
                        # Ingest data into array
                        # Requires try statement because some models change resolution in middle of simulation, for example: 
                        # pr_Amon_EC-Earth3-LR_piControl_r1i1p1f1_gr_199101-199101.nc vs pr_Amon_EC-Earth3-LR_piControl_r1i1p1f1_gr_241901-241912.nc
                        for dd in np.arange(len(time_arr)):
                            try:
                                ind = np.where(DatesMon==time_arr[dd])[0][0]
                                data[:,:,ind] = np.roll(np.flipud(ncdata[dd,:,:]),int(ncdata.shape[2]/2),axis=1)
                            except:
                                continue
                        del ncdata
                    
                    # Fix units
                    if variable=='pr':
                        data = data*10**6 # mm/month
                    elif variable=='tas':
                        data = data-273.15 # Degrees Celsius
                    
                    # Compute global mean time series
                    resy = 180/data.shape[0]
                    resx = 360/data.shape[1]
                    lat = 90-np.arange(180/resy)*resy-resy/2
                    lon = -180+np.arange(360/resx)*resx+resx/2 
                    xi, yi = np.meshgrid(lon, lat)
                    area_map = (40075*resx/360)**2*np.cos(np.deg2rad(yi)) # Grid-cell area in km2
                    data_mean_yr = np.zeros((data.shape[0],data.shape[1],len(Years)),dtype=np.single)*np.NaN
                    data_min_yr = np.zeros((data.shape[0],data.shape[1],len(Years)),dtype=np.single)*np.NaN
                    data_max_yr = np.zeros((data.shape[0],data.shape[1],len(Years)),dtype=np.single)*np.NaN
                    ts_mean_yr = np.zeros((len(Years),),dtype=np.single)*np.NaN
                    ts_min_yr = np.zeros((len(Years),),dtype=np.single)*np.NaN
                    ts_max_yr = np.zeros((len(Years),),dtype=np.single)*np.NaN
                    ts_mean_yr_land = np.zeros((len(Years),),dtype=np.single)*np.NaN
                    ts_min_yr_land = np.zeros((len(Years),),dtype=np.single)*np.NaN
                    ts_max_yr_land = np.zeros((len(Years),),dtype=np.single)*np.NaN
                    for yy in np.arange(len(Years)):
                        sel = DatesMon.year==Years[yy]
                        if sum(sel)==12:
                            data_mean_yr[:,:,yy] = np.round(np.mean(data[:,:,sel],axis=2),2)
                            data_min_yr[:,:,yy] = np.round(np.min(data[:,:,sel],axis=2),2)
                            data_max_yr[:,:,yy] = np.round(np.max(data[:,:,sel],axis=2),2)
                            ts_mean_yr[yy] = np.round(np.mean(data_mean_yr[:,:,yy]*area_map)/np.mean(area_map),2)
                            ts_min_yr[yy] = np.round(np.mean(data_min_yr[:,:,yy]*area_map)/np.mean(area_map),2)
                            ts_max_yr[yy] = np.round(np.mean(data_max_yr[:,:,yy]*area_map)/np.mean(area_map),2)
                            ts_mean_yr_land[yy] = np.round(np.mean(data_mean_yr[:,:,yy][mask_small]*area_map[mask_small])/np.mean(area_map[mask_small]),2)
                            ts_min_yr_land[yy] = np.round(np.mean(data_min_yr[:,:,yy][mask_small]*area_map[mask_small])/np.mean(area_map[mask_small]),2)
                            ts_max_yr_land[yy] = np.round(np.mean(data_max_yr[:,:,yy][mask_small]*area_map[mask_small])/np.mean(area_map[mask_small]),2)
                    print("Time elapsed is "+str(time.time()-t0)+" sec")
                
                    # Verify output time series
                    if sum(np.isnan(ts_mean_yr)==False)<10:
                        print('Time series too short, skipping')
                        continue

                    # Saving data to npz with default compression level (=6)
                    # Unfortunately, the default compression level cannot be changed, so this is slow
                    print('Saving data to npz')
                    t0 = time.time()
                    np.savez_compressed(os.path.join(config['folder_out'],'climate_model_data',scenario+'_'+model+'_'+member+'_'+variable+'.npz'),\
                        data=data,data_mean_yr=data_mean_yr,data_min_yr=data_min_yr,data_max_yr=data_max_yr,\
                        ts_mean_yr=ts_mean_yr,ts_min_yr=ts_min_yr,ts_max_yr=ts_max_yr,\
                        ts_mean_yr_land=ts_mean_yr_land,ts_min_yr_land=ts_min_yr_land,ts_max_yr_land=ts_max_yr_land,\
                        Years=Years,DatesMon=DatesMon,mask_small=mask_small)
                    print("Time elapsed is "+str(time.time()-t0)+" sec")
                    
                    if (ee==0) & (generate_figures==True):
                        print('Generating verification figures')
                        t0 = time.time()
                        plt.figure(1)
                        plt.plot(np.squeeze(data[50,50,:]))
                        plt.savefig(os.path.join(config['folder_out'],'climate_model_data','figures',variable+'_'+scenario+'_'+model+'_'+member+'_ts.png'))                    
                        plt.figure(2)
                        if variable=='pr':
                            plt.imshow(data[:,:,200],vmin=-20,vmax=200)
                        elif variable=='tas':
                            plt.imshow(data[:,:,200],vmin=-40,vmax=40)
                        else:
                            plt.imshow(data[:,:,200])
                        plt.colorbar()
                        plt.savefig(os.path.join(config['folder_out'],'climate_model_data','figures',variable+'_'+scenario+'_'+model+'_'+member+'_map.png'))
                        print("Time elapsed is "+str(time.time()-t0)+" sec")
                        plt.close('all')
                    
                    del data
                    gc.collect()
                    
    pdb.set_trace()
    
    
if __name__ == '__main__':
    main()