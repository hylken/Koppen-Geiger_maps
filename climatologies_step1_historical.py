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
import pandas as pd
import numpy as np
import tools
from datetime import datetime
from netCDF4 import Dataset
import warnings
from skimage.transform import resize
import h5py
import gc
import tables
import matplotlib.pyplot as plt

def main():
    
    
    #==============================================================================
    #   Settings
    #==============================================================================

    np.set_printoptions(suppress=True)
    config = tools.load_config(sys.argv[1])
    koppen_table = pd.read_csv(os.path.join('assets','koppen_table.csv'))
    
    # Create land-sea mask
    filepath = os.path.join(config['folder_maps'],'WorldClim_V21_Pmean_01.mat')
    f = h5py.File(filepath)
    mask = np.isnan(np.transpose(np.array(f['DATA'],dtype=np.single)))
    f.close()
    
    # Loop over periods
    for period_historical in config['periods_historical']:
        out_dir = os.path.join(config['folder_out'],'climatologies',str(period_historical[0])+'_'+str(period_historical[1]))
        if os.path.isdir(out_dir)==False:
            os.makedirs(out_dir)
        
        
        #==============================================================================
        #   Load observed P and T data, to adjust the high-resolution climatologies
        #==============================================================================

        print('-------------------------------------------------------------------------------')
        print('Loading monthly CRU T data') 
        t0 = time.time()
        DatesMon = pd.date_range(start=datetime(1850,1,1), end=datetime(2021,12,31), freq='MS')
        obs_data = {}
        obs_data['Temp'] = np.zeros((360,720,len(DatesMon)),dtype=np.single)*np.NaN
        filepath = glob.glob(os.path.join(config['folder_dataraw'], 'CRU_TS', 'cru_ts*.1901.20*.tmp.dat.nc'))[0]
        with Dataset(filepath) as dset:
            data = np.flip(np.array(dset.variables['tmp'][:],dtype=np.single),axis=1)
            data[data>100] = np.NaN
            dates = pd.date_range(start=datetime(1901,1,1), end=datetime(2049,12,1), freq='MS')
            dates = dates[:data.shape[0]]
            for yy in np.arange(len(DatesMon)):
                sel = dates==DatesMon[yy]
                if sum(sel)==1:
                    obs_data['Temp'][:,:,yy] = np.squeeze(data[sel,:,:])
        print("Time elapsed is "+str(time.time()-t0)+" sec")
        
        print('-------------------------------------------------------------------------------')
        print('Loading monthly GPCC data') 
        t0 = time.time()
        obs_data['P'] = np.zeros((720,1440,len(DatesMon)),dtype=np.single)*np.NaN
        files = glob.glob(os.path.join(config['folder_dataraw'],'GPCC','*.nc'))
        for ii in np.arange(len(files)):
            dset = Dataset(files[ii])
            precip = np.array(dset.variables['precip'][:],dtype=np.single)
            dset.close()
            precip[precip<0] = np.NaN
            dates = pd.date_range(
                start=datetime(int(os.path.basename(files[ii]).split('_')[4]),1,1),
                end=datetime(int(os.path.basename(files[ii]).split('_')[5]),12,1),
                freq='MS')
            for jj in np.arange(len(dates)):
                try:
                    ind = np.where(DatesMon==dates[jj])[0][0]
                    obs_data['P'][:,:,ind] = precip[jj,:,:]
                except:
                    pass
        del precip
        print("Time elapsed is "+str(time.time()-t0)+" sec")
        
        
        #==============================================================================
        #   Generate high-resolution historical climatologies; each P and T dataset 
        #   combination is an ensemble member
        #==============================================================================
        
        count = 0
        for Pdataset in config['Pdatasets'].keys():
            for Tdataset in config['Tdatasets'].keys():
                count += 1
                print('-------------------------------------------------------------------------------') 
                print(str(period_historical[0])+'-'+str(period_historical[1])+' ensemble member '+str(count).zfill(3))

                # Check if output file already exists
                ncout = os.path.join(out_dir,str(count).zfill(3)+'.nc')
                if os.path.isfile(ncout):
                    print('Already processed')
                    continue
                t0 = time.time()
                
                # Loop over variables and months
                for vv in np.arange(len(config['vars'])):
                    varname = config['vars'][vv][0]
                    for month in np.arange(1,13):
                        
                        # Load high-resolution reference climatology
                        if varname=='P':
                            filepath = os.path.join(config['folder_maps'],Pdataset+'_Pmean_'+str(month).zfill(2)+'.mat')
                            print('Loading '+filepath)
                            f = h5py.File(filepath)
                            data = np.transpose(np.array(f['DATA'],dtype=np.single))
                            reference_map = resize(data,config['mapsize'],order=1,mode='constant',anti_aliasing=False)
                            reference_map[mask] = np.NaN
                            f.close()
                            reference_period = config['Pdatasets'][Pdataset]
                        elif varname=='Temp':
                            filepath = os.path.join(config['folder_maps'],Tdataset+'_Tmean_'+str(month).zfill(2)+'.mat')
                            print('Loading '+filepath)
                            f = h5py.File(filepath)
                            data = np.transpose(np.array(f['DATA'],dtype=np.single))
                            reference_map = resize(data,config['mapsize'],order=1,mode='constant',anti_aliasing=False)
                            reference_map[mask] = np.NaN
                            f.close()
                            reference_period = config['Tdatasets'][Tdataset]
                                 
                        # Produce change map to adjust climatology to period in 
                        # question based on the monthly observed climate data                            
                        change = tools.produce_change_map(
                            reference_map = reference_map,
                            reference_period = reference_period,
                            monthly_data = obs_data[varname],
                            monthly_dates = DatesMon,
                            target_period = period_historical,
                            target_month = month,
                            varname = varname,
                            change_offset = 3,
                            change_limits = config['vars'][vv][4]
                            )
                            
                        if np.isnan(np.nanmean(change['target_map'])):
                            raise ValueError('Only NaNs, something is wrong')
                        
                        # Save temporally-adjusted, high-resolution climatology
                        print('Writing to '+ncout)
                        tools.write_to_netcdf_3d(ncout,change['target_map'],config['vars'][vv][1],config['vars'][vv][2],month,1)
                        del reference_map
                        
                        # Attempt to fix random "There are 150 HDF5 objects open!" errors
                        tables.file._open_files.close_all()
                        try:
                            del f
                        except:
                            pass
                            
                print("Time elapsed is "+str(time.time()-t0)+" sec")
                
                gc.collect()
            
        del obs_data
        
                
        print('-------------------------------------------------------------------------------') 
        print(str(period_historical[0])+'-'+str(period_historical[1])+' ensemble mean and std')
        t0 = time.time()
        tools.compute_ens_mean_std(out_dir,config['vars'],config['mapsize'],3600,True,mask)
        print("Time elapsed is "+str(time.time()-t0)+" sec") 
              
              
        print('-------------------------------------------------------------------------------') 
        print(str(period_historical[0])+'-'+str(period_historical[1])+' KG classification and uncertainty')
        t0 = time.time()
        tools.compute_kg_maps(out_dir,koppen_table,config['mapsize'],3600,True,mask)
        print("Time elapsed is "+str(time.time()-t0)+" sec") 
        
    pdb.set_trace()
    
    
if __name__ == '__main__':
    main()