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
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import warnings
import zipfile
import rasterio

def main():


    #==============================================================================
    #   Settings
    #==============================================================================
    
    np.set_printoptions(suppress=True)
    config = tools.load_config(sys.argv[1])
    script_name = os.path.basename(sys.argv[0]).replace('.py', '')
    koppen_table = pd.read_csv(os.path.join('assets','koppen_table.csv'))


    #==============================================================================
    #   For each high-resolution netCDF, make low resolution, upscaled version, and
    #   save as netCDF
    #==============================================================================

    folders = [
        os.path.join(config['folder_out'], 'climatologies_step1_historical'),
        os.path.join(config['folder_out'], 'climatologies_step2_future')
    ]

    subfolders = [
        os.path.join(parent, name)
        for parent in folders
        for name in os.listdir(parent)
        if os.path.isdir(os.path.join(parent, name))
    ]

    # Loop over all files and convert to netCDF
    for subfolder in subfolders:
        files = glob.glob(os.path.join(subfolder,'*.nc'))
        for file in files:
            suffix_full = str(180/config['mapsize'][0]).replace('.','p')[:10]
            for upscale_mapsize in config['upscale_mapsizes']:
                suffix_resample = str(180/upscale_mapsize[0]).replace('.','p')                
                
                # Resample precipitation and air temperature climatologies
                if (file=='ensemble_mean_'+suffix_full+'.nc') | (file=='ensemble_std_'+suffix_full+'.nc'): 
                    t0 = time.time()
                    file_new = file.replace(suffix_full,suffix_resample)
                    ncout = os.path.join(config['folder_out'],script_name,os.path.basename(os.path.normpath(subfolder)),file_new)
                    os.makedirs(os.path.dirname(ncout), exist_ok=True)
                    if (os.path.isfile(ncout)) & (config['skip_existing']==True):
                        continue
                    elif (os.path.isfile(ncout)) & (config['skip_existing']==False):
                        os.remove(ncout)
                    print('Creating '+ncout)
                    dset = Dataset(os.path.join(root,file))
                    for month in np.arange(1,13):
                        for vv in np.arange(len(config['vars'])):
                            data = np.array(dset.variables[config['vars'][vv][1]][month-1,:,:])
                            data = tools.mapresize(data,upscale_mapsize,measure='mean',nantol=0.75)
                            tools.write_to_netcdf_3d(ncout,data,config['vars'][vv][1],config['vars'][vv][2],month,1)
                    dset.close()                    
                    print("Time elapsed is "+str(time.time()-t0)+" sec")
                
                # Resample Koppen-Geiger maps
                if file=='koppen_geiger_'+suffix_full+'.nc':
                    t0 = time.time()
                    file_new = file.replace(suffix_full,suffix_resample)
                    ncout = os.path.join(config['folder_out'],script_name,os.path.basename(os.path.normpath(subfolder)),file_new)
                    if (os.path.isfile(ncout)) & (config['skip_existing']==True):
                        continue
                    elif (os.path.isfile(ncout)) & (config['skip_existing']==False):
                        os.remove(ncout)
                    print('Creating '+ncout)
                    dset = Dataset(os.path.join(root,file))
                    kg_class = np.array(dset.variables['kg_class'][:])
                    kg_class = tools.mapresize(kg_class,upscale_mapsize,measure='mode',nantol=0.75,nanint=kg_class[0,0])
                    kg_confidence = np.array(dset.variables['kg_confidence'][:])
                    kg_confidence = tools.mapresize(kg_confidence,upscale_mapsize,measure='mean',nantol=0.75)
                    dset.close()
                    tools.write_to_netcdf_2d(ncout,kg_class,'kg_class','',1)
                    tools.write_to_netcdf_2d(ncout,kg_confidence,'kg_confidence','%',1)
                    print("Time elapsed is "+str(time.time()-t0)+" sec")
                
    
    #==============================================================================
    #   Convert the Koppen-Geiger maps to geoTIFF
    #==============================================================================

    # Create Koppen-Geiger colormap for geoTIFF
    cmap = {}
    cmap[0] = (255,255,255)
    for ii in np.arange(koppen_table.shape[0]):
        cmap[koppen_table['Class'][ii]] = tuple(koppen_table[['Red','Green','Blue']].iloc[ii])

    # Loop over all files and convert to geoTIFF
    for subfolder in subfolders:
        files = glob.glob(os.path.join(subfolder,'*koppen_geiger*.nc'))
        for file in files:
            t0 = time.time()
            ncout = os.path.join(config['folder_out'],script_name,os.path.basename(os.path.normpath(subfolder)),file.replace('.nc','.tif'))
            if (os.path.isfile(ncout)) & (config['skip_existing']==True):
                continue
            elif (os.path.isfile(ncout)) & (config['skip_existing']==False):
                os.remove(ncout)
            print('Creating '+ncout)
            dset = Dataset(os.path.join(root,file))
            data = np.array(dset.variables['kg_class'][:])
            dset.close()
            tools.write_to_geotiff(ncout,data,cmap,0)
            print("Time elapsed is "+str(time.time()-t0)+" sec")
   
    
    #==============================================================================
    #   Make zip files
    #==============================================================================

    # Delete existing zip files
    filelist = glob.glob(os.path.join(config['folder_out'],script_name,'*.zip'))
    for filepath in filelist:
        print('Deleting '+filepath)
        os.remove(filepath)
        
    # Create zip file with netCDF Koppen-Geiger maps
    t0 = time.time()
    zip_file = os.path.join(config['folder_out'],script_name,'koppen_geiger_nc.zip')
    folder = os.path.join(config['folder_out'],script_name)
    pattern = 'koppen_geiger*.nc'
    compress_type = zipfile.ZIP_DEFLATED
    print('Creating '+zip_file)
    tools.zip_folder(zip_file,folder,pattern,compress_type)
    print("Time elapsed is "+str(time.time()-t0)+" sec")
    
    # Create zip file with geoTIFF Koppen-Geiger maps
    t0 = time.time()
    zip_file = os.path.join(config['folder_out'],script_name,'koppen_geiger_tif.zip')
    folder = os.path.join(config['folder_out'],script_name)
    pattern = 'koppen_geiger*.tif'
    compress_type = zipfile.ZIP_DEFLATED
    print('Creating '+zip_file)
    tools.zip_folder(zip_file,folder,pattern,compress_type)
    print("Time elapsed is "+str(time.time()-t0)+" sec")

    # Create zip file with precipitation and air temperature climatologies
    for mapsize in config['upscale_mapsizes']+[config['mapsize']]:
        suffix = str(180/mapsize[0]).replace('.','p')[:10]
        t0 = time.time()
        zip_file = os.path.join(config['folder_out'],script_name,'climate_data_'+suffix+'.zip')
        folder = os.path.join(config['folder_out'],script_name)
        pattern = 'ensemble_*'+suffix+'*.nc'
        compress_type = zipfile.ZIP_DEFLATED
        print('Creating '+zip_file)
        tools.zip_folder(zip_file,folder,pattern,compress_type)
        print("Time elapsed is "+str(time.time()-t0)+" sec")
    
    # Upload zip files to server
    if config['perform_sync']:
        print('Syncing to remote')
        tools.sync_data(
            config['sync_cmd'],
            dir_local = os.path.join(config['folder_out'],script_name),
            dir_remote = '')
    
    print("IMPORTANT: Don't forget to manually add the legend.txt file to the two Koppen-Geiger zip files!")
    
    pdb.set_trace()
    
    
if __name__ == '__main__':
    main()