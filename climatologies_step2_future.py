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
import shutil
import pandas as pd
import numpy as np
import tools
from datetime import datetime
from netCDF4 import Dataset
import warnings
import gc
import tables

def main():


    #==============================================================================
    #   Settings
    #==============================================================================
    
    np.set_printoptions(suppress=True)
    config = tools.load_config(sys.argv[1])
    koppen_table = pd.read_csv(os.path.join('assets','koppen_table.csv'))
    scenarios = ['ssp119','ssp126','ssp245','ssp370','ssp434','ssp460','ssp585']
    
    # Create land-sea mask
    filepath = os.path.join(config['folder_maps'],'WorldClim_V21_Pmean_01.mat')
    f = h5py.File(filepath)
    mask = np.isnan(np.transpose(np.array(f['DATA'],dtype=np.single)))
    f.close()
    
    # Loop over periods and scenarios
    for period_future in config['periods_future']:
        for scenario in scenarios:
            out_dir = os.path.join(config['folder_out'],'climatologies', \
                str(period_future[0])+'_'+str(period_future[1]),scenario)
            if os.path.isdir(out_dir)==False:
                os.makedirs(out_dir)
            
            # Load model subset information
            with np.load(os.path.join(config['folder_stats'],'projected_change.npz')) as dset_projected_change:
                models = dset_projected_change['models']
                model_subset = dset_projected_change['model_subset']
            
            # Loop over models
            for mm in np.arange(len(models)):
                model = models[mm]
                ncout = os.path.join(out_dir,model+'.nc')
                
                print('-------------------------------------------------------------------------------') 
                print(str(period_future[0])+'-'+str(period_future[1])+' '+scenario+' '+model)

                # Only process model if in Model Subset, if not in Model Subset, 
                # delete output file (if created by previous run)
                if model_subset[mm]==False:
                    print('Not in model subset, skipping')
                    if os.path.isfile(ncout):
                        print('Deleting output from previous run')
                        os.remove(ncout)
                    continue
                
                # Check if output file already exists
                if os.path.isfile(ncout):
                    print('Already processed, skipping')
                    continue
                                
                
                #==============================================================================
                #   Load climate model data, average over ensemble members
                #==============================================================================

                t0 = time.time()
                
                # Get list of ensemble members
                files = glob.glob(os.path.join(config['folder_out'], 'climate_model_data', scenario+'_'+model+'_*_tas.npz'))
                members = [os.path.basename(x).split('_')[2] for x in files]
                if len(members)==0:
                    print('No ensemble members available, skipping')
                    continue
                    
                # Loop over variables and ensemble members
                DatesMon = pd.date_range(start=datetime(1991,1,1), end=datetime(2099,12,31), freq='MS')
                tmp = np.load(files[0])
                dims = tmp['data'].shape
                sim_data = {}
                for vv in np.arange(len(config['vars'])):
                    varname = config['vars'][vv][0]
                    sim_data[varname] = np.zeros((dims[0],dims[1],len(DatesMon),len(members)),dtype=np.single)*np.NaN
                    for ee in np.arange(len(members)):
                        member = members[ee]
                        print(member)
                        
                        # Check if both historical and projection data exist
                        cmip6_hist_file = os.path.join(config['folder_out'], 'climate_model_data', 'historical_'+model+'_'+member+'_'+config['vars'][vv][3]+'.npz')
                        cmip6_future_file = os.path.join(config['folder_out'], 'climate_model_data', scenario+'_'+model+'_'+member+'_'+config['vars'][vv][3]+'.npz')
                        if (os.path.isfile(cmip6_hist_file)==False) | (os.path.isfile(cmip6_future_file)==False):
                            print(varname+": historical or future CMIP6 file doesn't exist, skipping")
                            continue
                        
                        # Load historical data
                        with np.load(cmip6_hist_file,allow_pickle=True) as hist:
                            months = tools.periods_to_datetimes(hist['DatesMon'])
                            data = hist['data']
                            for yy in np.arange(len(DatesMon)):
                                sel = months==DatesMon[yy]
                                if sum(sel)==1:
                                    sim_data[varname][:,:,yy,ee] = np.squeeze(data[:,:,sel])
                            
                        # Load projection 
                        with np.load(cmip6_future_file,allow_pickle=True) as future:
                            months = tools.periods_to_datetimes(future['DatesMon'])
                            data = future['data']
                            for yy in np.arange(len(DatesMon)):
                                sel = months==DatesMon[yy]
                                if sum(sel)==1:
                                    sim_data[varname][:,:,yy,ee] = np.squeeze(data[:,:,sel])
                            
                        # Check data temporal completeness 
                        # If incomplete, make everything NaN
                        if sum(np.isnan(sim_data[varname][0,0,:,ee]))>0:
                            print(varname+' time series incomplete, skipping')
                            sim_data[varname][:,:,:,ee] = np.NaN
                            continue
                            
                    # Compute average over all ensemble members
                    sim_data[varname] = np.nanmean(sim_data[varname],axis=3)
                
                # Check if data for both P and T
                if np.isnan(sim_data['P'][0,0,0]) | np.isnan(sim_data['Temp'][0,0,0]):
                    print('No P or Temp data available, skipping')
                    continue
                
                
                #==============================================================================
                #   Generate high-resolution future climatologies, ensemble mean and standard 
                #   deviation, and the final Koppen-Geiger maps
                #==============================================================================
                
                # Loop over variables and months
                for vv in np.arange(len(config['vars'])):
                    varname = config['vars'][vv][0]
                    for month in np.arange(1,13):
                    
                        # Load high-res historic reference climatology
                        period_historic = config['periods_historical'][-1]
                        suffix = str(180/config['mapsize'][0]).replace('.','p')[:10]
                        dset = Dataset(os.path.join(config['folder_out'],'climatologies', \
                            str(period_historic[0])+'_'+str(period_historic[1]),'ensemble_mean_'+suffix+'.nc'))
                        data = np.array(dset.variables[config['vars'][vv][1]][month-1,:,:],dtype=np.single)
                        dset.close()
                        
                        # Produce change map to adjust climatology to period in 
                        # question based on the monthly climate model data
                        change = tools.produce_change_map(
                            reference_map = data,
                            reference_period = period_historic,
                            monthly_data = sim_data[varname],
                            monthly_dates = DatesMon,
                            target_period = period_future,
                            target_month = month,
                            varname = varname,
                            change_offset = 3,
                            change_limits = config['vars'][vv][4]
                            )
                        
                        if np.isnan(np.nanmean(change['target_map'])):
                            raise ValueError('Only NaNs, something is wrong')
                        
                        # Save temporally-adjusted high-res future climatology
                        tools.write_to_netcdf_3d(ncout,change['target_map'],config['vars'][vv][1],config['vars'][vv][2],month,1)
                        del data

                        # Attempt to fix random "There are 150 HDF5 objects open!" errors
                        tables.file._open_files.close_all()
                                
                del sim_data    
                
                print("Time elapsed is "+str(time.time()-t0)+" sec")
                
            gc.collect()
        
        
            print('-------------------------------------------------------------------------------') 
            print(str(period_future[0])+'-'+str(period_future[1])+' '+scenario+' ensemble mean and std')
            t0 = time.time()
            tools.compute_ens_mean_std(out_dir,config['vars'],config['mapsize'],3600,True,mask)
            print("Time elapsed is "+str(time.time()-t0)+" sec") 
        
    
            print('-------------------------------------------------------------------------------') 
            print(str(period_future[0])+'-'+str(period_future[1])+' '+scenario+' KG classification and uncertainty')
            t0 = time.time()
            tools.compute_kg_maps(out_dir,koppen_table,config['mapsize'],3600,True,mask)
            print("Time elapsed is "+str(time.time()-t0)+" sec") 
            
    pdb.set_trace()
    
    
if __name__ == '__main__':
    main()