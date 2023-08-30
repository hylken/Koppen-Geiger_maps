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
import pickle
import matplotlib as plt

def main():


    #==============================================================================
    #   Settings
    #==============================================================================
    
    warnings.filterwarnings('ignore')
    np.set_printoptions(suppress=True)
    config = tools.load_config(sys.argv[1])
    koppen_table = pd.read_csv(os.path.join('assets','koppen_table.csv'))


    #==============================================================================
    #   Compute areas covered by major KG classes and transitions
    #==============================================================================

    scenarios = ['ssp119','ssp126','ssp245','ssp370','ssp434','ssp460','ssp585']
    periods = config['periods_historical']+config['periods_future']
    mapsize = config['upscale_mapsizes'][0]    
    
    # Area map (million km2)
    res = 180/config['upscale_mapsizes'][0][0]
    xi, yi = np.meshgrid(np.arange(-180+res/2,180+res/2,res), np.arange(90-res/2,-90-res/2,-res))
    area_map = 10**-6*(40075*res/360)**2*np.cos(np.deg2rad(yi))

    # Loop over scenarios and periods
    df_kg_major_change_prct = pd.DataFrame(np.zeros((len(scenarios),2))*np.NaN,index=scenarios,columns=['1901-1930 to 1991-2020','1991-2020 to 2071-2100'])    
    for scenario in scenarios:
        print('===============================================================================')
        print('Compute areas covered by major KG classes and transitions for '+scenario)
        kg_maps = np.zeros((mapsize[0],mapsize[1],len(periods)),dtype=np.single)*np.NaN
        for pp in np.arange(len(periods)):
            period = periods[pp]
            
            # Load global Koppen-Geiger map
            suffix = str(180/config['upscale_mapsizes'][0][0]).replace('.','p')
            ncfile1 = os.path.join(config['folder_out'],'climatologies',str(period[0])+'_'+str(period[1]),'koppen_geiger_'+suffix+'.nc')
            ncfile2 = os.path.join(config['folder_out'],'climatologies',str(period[0])+'_'+str(period[1]),scenario,'koppen_geiger_'+suffix+'.nc')
            if os.path.isfile(ncfile1):
                ncfile = ncfile1
            elif os.path.isfile(ncfile2):
                ncfile = ncfile2
            else:
                raise Exception("Unable to load map")
            print('loading '+ncfile)
            dset = Dataset(ncfile)
            kg_class = np.array(dset.variables['kg_class'][:]).astype(int)
            dset.close()
            
            # Compute map of major classes
            kg_major = np.zeros(kg_class.shape,dtype=int)
            for ii in np.arange(koppen_table.shape[0]):
                mask = kg_class==koppen_table['Class'][ii]
                kg_major[mask] = koppen_table['Major'][ii]
            
            # Insert map into array
            kg_maps[:,:,pp] = kg_major
        
        # Discard Antarctica
        mask = np.min(kg_maps,axis=2)==0
        for pp in np.arange(len(periods)):
            kg_maps[:,:,pp][mask] = 0
        
        # Make tables with area for major classes and transitions
        df_kg_major_area_pct = pd.DataFrame(np.zeros((len(periods),6))*np.NaN,index=None,columns=['Period','A','B','C','D','E'])
        df_kg_major_area_mm2 = pd.DataFrame(np.zeros((len(periods),6))*np.NaN,index=None,columns=['Period','A','B','C','D','E'])
        df_transitions_mm2 = pd.DataFrame(np.zeros((1000,5))*np.NaN,index=None,columns=['From','To','Source','Target','Area'])
        count = 0
        for ll in np.arange(len(periods)):
        
            # Compute areas
            mask_land = kg_maps[:,:,0]!=0
            df_kg_major_area_pct.iloc[ll,0] = str(periods[ll])
            df_kg_major_area_mm2.iloc[ll,0] = str(periods[ll])
            for cl in np.arange(1,6):
                mask_cl = kg_maps[:,:,ll]==cl
                df_kg_major_area_pct.iloc[ll,cl] = 100*np.sum(area_map[mask_cl])/np.sum(area_map[mask_land])
                df_kg_major_area_mm2.iloc[ll,cl] = np.sum(area_map[mask_cl])
                
            # Don't compute transitions for last period
            if ll==len(periods)-1:
                continue
            
            # Loop over transitions
            for source in np.arange(1,6):
                for target in np.arange(1,6): 
                    mask = ((kg_maps[:,:,ll]==source)==True) & ((kg_maps[:,:,ll+1]==target)==True)
                    mask_area = np.sum(area_map[mask])
                    
                    # Add transition to table
                    if mask_area>0.05:
                        df_transitions_mm2.iloc[count,:] = np.array([str(periods[ll]),str(periods[ll+1]),source,target,mask_area])
                        count +=1
                    
        # Save results
        df_kg_major_area_pct.to_csv(os.path.join(config['folder_stats'],'climatologies',scenario+'_kg_major_area_pct.csv'),index=False)
        df_kg_major_area_mm2.to_csv(os.path.join(config['folder_stats'],'climatologies',scenario+'_kg_major_area_mm2.csv'),index=False)
        df_transitions_mm2.to_csv(os.path.join(config['folder_stats'],'climatologies',scenario+'_transitions_mm2.csv'),index=False)
        
        # Compute percentage of land surface that changes
        mask_land = kg_maps[:,:,0]!=0        
        ll_1901 = np.array(periods)[:,0]==1901
        ll_1991 = np.array(periods)[:,0]==1991
        ll_2071 = np.array(periods)[:,0]==2071
        diff1 = np.squeeze(np.abs(kg_maps[:,:,ll_1991]-kg_maps[:,:,ll_1901]))>0
        diff2 = np.squeeze(np.abs(kg_maps[:,:,ll_2071]-kg_maps[:,:,ll_1991]))>0
        diff1 = 100*np.mean(diff1[mask_land]*area_map[mask_land])/np.mean(area_map[mask_land])
        diff2 = 100*np.mean(diff2[mask_land]*area_map[mask_land])/np.mean(area_map[mask_land])
        df_kg_major_change_prct.loc[scenario] = [diff1,diff2]
        
    # Save results
    df_kg_major_change_prct.to_csv(os.path.join(config['folder_stats'],'climatologies','kg_major_change_prct.csv'),index=False)
   
    
    #==============================================================================
    #   Load station data and compute KG classes
    #==============================================================================

    if os.path.isfile(os.path.join(config['folder_out'],'climatologies_validation','station_data.pickle')):
        print('===============================================================================')
        print('Loading existing station data file')        
        t = time.time()
        station_data = pickle.load(open(os.path.join(config['folder_out'],'climatologies_validation','station_data.pickle'),'rb'))
        print("Time elapsed is "+str(time.time()-t)+" sec")

    if os.path.isfile(os.path.join(config['folder_out'],'climatologies_validation','station_data.pickle'))==False:
        print('===============================================================================')
        print('Loading station data')
        dates_daily = pd.date_range(start=pd.to_datetime(datetime(1900,1,1)),end=pd.to_datetime('today'), freq='D')
        station_files = glob.glob(os.path.join(config['folder_station'],'*.mat'))
        station_data = {}
        station_data['lat'] = np.zeros((len(station_files)),dtype=np.single)*np.NaN
        station_data['lon'] = np.zeros((len(station_files)),dtype=np.single)*np.NaN
        station_data['name'] = np.zeros((len(station_files)),dtype=object)*np.NaN
        station_data['T_monthly_clim'] = np.zeros((len(station_files),len(config['periods_historical']),12),dtype=np.single)*np.NaN
        station_data['P_monthly_clim'] = np.zeros((len(station_files),len(config['periods_historical']),12),dtype=np.single)*np.NaN
        station_data['Class'] = np.zeros((len(station_files),len(config['periods_historical'])),dtype=np.single)*np.NaN
        station_data['Major'] = np.zeros((len(station_files),len(config['periods_historical'])),dtype=np.single)*np.NaN
       
        # Loop over stations
        for ii in np.arange(len(station_files)): 
            print('Loading station '+str(ii))
            t0 = time.time()
            if ii in np.linspace(0,len(station_files),20).astype(int): 
                print(str(np.round(100*ii/len(station_files)))+' % completed')
            
            # Read latitude, longitude, and name from mat file
            station_data['lat'][ii] = tools.readmatfile(station_files[ii],'StationCoords/Lat')[0]
            station_data['lon'][ii] = tools.readmatfile(station_files[ii],'StationCoords/Lon')[0]
            station_data['name'][ii] = os.path.basename(station_files[ii])[:-4]
            if np.isnan(station_data['lat'][ii]+station_data['lon'][ii]) | (abs(station_data['lat'][ii])>89.5) | (abs(station_data['lon'][ii])>179.5): 
                continue
            
            # Read data from mat file
            vars = ['PRCP','TMIN','TMAX','TAVG']
            statdata = {}
            for var in vars:                                
                statdata[var] = np.zeros((len(dates_daily),1))*np.NaN
                try:                    
                    statdata[var] = tools.readmatfile(station_files[ii],var).flatten().reshape(-1,1)
                    if len(statdata[var])<len(dates_daily):
                        statdata[var] = np.concatenate((statdata[var],np.zeros((len(dates_daily),1))*np.NaN),axis=0)
                    statdata[var] = statdata[var][:len(dates_daily)] 
                    if (var=='PRCP'):
                        statdata[var] = statdata[var]*30.4 # Compute monthly total
                    if (var=='PRCP') & ("GSOD" in station_name[ii]):
                        statdata[var] = tools.eliminate_trailing_zeros(statdata[var]) # Eliminate erroneous zero precipitation in GSOD stations                    
                except: 
                    continue                
            
            # If TAVG not available, compute mean daily air temperature from TMIN and TMAX
            sel = np.isnan(statdata['TAVG'])
            statdata['TAVG'][sel] = ((statdata['TMIN']+statdata['TMAX'])/2)[sel]
            
            # Compute monthly precipitation and air temperature climatologies and Koppen-Geiger class
            for pp in np.arange(len(config['periods_historical'])):
                period = config['periods_historical'][pp]
                sel = (dates_daily>=datetime(period[0],1,1)) & (dates_daily<=datetime(period[1],12,31))
                T_monthly_clim = tools.compute_monthly_climatology(statdata['TAVG'][sel],dates_daily[sel])
                P_monthly_clim = tools.compute_monthly_climatology(statdata['PRCP'][sel],dates_daily[sel])
                station_data['T_monthly_clim'][ii,pp,:] = T_monthly_clim
                station_data['P_monthly_clim'][ii,pp,:] = P_monthly_clim
                KG_dict = tools.koppen_geiger(T_monthly_clim.reshape(12,1,1),P_monthly_clim.reshape(12,1,1),koppen_table)
                station_data['Class'][ii,pp] = KG_dict['Class']
                station_data['Major'][ii,pp] = KG_dict['Major']
                
            print(station_data['Major'][ii,:])
            
            print("Time elapsed is "+str(time.time()-t0)+" sec")
            
        if os.path.isdir(os.path.join(config['folder_out'],'climatologies_validation'))==False:
            os.makedirs(os.path.join(config['folder_out'],'climatologies_validation'))        
        with open(os.path.join(config['folder_out'],'climatologies_validation','station_data.pickle'), 'wb') as f:
            pickle.dump(station_data, f)
        

    #==============================================================================
    #   Count number of stations for each provider
    #==============================================================================

    # Get list of providers
    provider = []
    for station_name in station_data['name']:
        provider.append(station_name.split('_')[0])
    providers = np.unique(provider)
    print(providers)
    
    # Get number of stations for each provider
    for provider in providers:
        res = [i for i in station_data['name'] if provider in i]
        print(provider+' '+str(len(res))+' stations')
    
    
    #==============================================================================
    #   Compute classification accuracy for each period for both 30 classes and 
    #   for the major classes
    #==============================================================================
   
    print('===============================================================================')
    print('Computing accuracy for historical periods')
    nperiods = len(config['periods_historical'])
    df_accuracy = pd.DataFrame(np.zeros((nperiods,6))*np.NaN,index=None,columns=['Period','nobs','Class','Major','conf_correct','conf_incorrect'])
    for pp in np.arange(nperiods):
        period = config['periods_historical'][pp]
        print(period)
        
        # Load global Koppen-Geiger map
        suffix = str(180/config['mapsize'][0]).replace('.','p')
        ncfile = os.path.join(config['folder_out'],'climatologies',str(period[0])+'_'+str(period[1]),'koppen_geiger_'+suffix+'.nc')
        dset = Dataset(ncfile)
        kg_class = np.array(dset.variables['kg_class'][:]).astype(int)
        kg_confidence = np.array(dset.variables['kg_confidence'][:]).astype(int)
        dset.close()
        
        # Compute map of major classes
        kg_major = np.zeros(kg_class.shape,dtype=int)
        for ii in np.arange(koppen_table.shape[0]):
            mask = kg_class==koppen_table['Class'][ii]
            kg_major[mask] = koppen_table['Major'][ii]
        
        # Convert lat/lon to row/col
        ys = config['mapsize'][0]*(90-station_data['lat'])/180-0.5
        ys = np.round(ys).astype(int)
        ys[(ys<0) | (ys>=config['mapsize'][0])] = 0
        xs = config['mapsize'][1]*(180+station_data['lon'])/360-0.5
        xs = np.round(xs).astype(int)
        xs[(xs<0) | (xs>=config['mapsize'][1])] = 0
              
        # Compute accuracy using station data
        valid = (np.isnan(station_data['Class'][:,pp])==False) \
            & (kg_class[ys,xs]!=0) \
            & (np.abs(station_data['lat'])<89.5) \
            & (np.abs(station_data['lon'])<179.5) \
            & (np.isnan(station_data['lat'])==False) \
            & (np.isnan(station_data['lon'])==False)
        df_accuracy['Period'][pp] = str(period)
        df_accuracy['nobs'][pp] = np.sum(valid)
        df_accuracy['Class'][pp] = 100*np.sum(station_data['Class'][:,pp][valid]==kg_class[ys,xs][valid])/np.sum(valid)
        df_accuracy['Major'][pp] = 100*np.sum(station_data['Major'][:,pp][valid]==kg_major[ys,xs][valid])/np.sum(valid)
        
        # Compute confidence level of correct and incorrect classifications
        correct = (valid) & (station_data['Class'][:,pp]==kg_class[ys,xs])        
        incorrect = (valid) & (station_data['Class'][:,pp]!=kg_class[ys,xs])
        df_accuracy['conf_correct'][pp] = np.mean(kg_confidence[ys,xs][correct])
        df_accuracy['conf_incorrect'][pp] = np.mean(kg_confidence[ys,xs][incorrect])
        
        # Performance comparison inside and outside US
        US = (station_data['lat']<50) & (station_data['lat']>30) & (station_data['lon']<-66) & (station_data['lon']>-126)
        sel1 = (valid) & (US)
        sel2 = (valid) & (~US)
        acc_US = 100*np.sum(station_data['Class'][:,pp][sel1]==kg_class[ys,xs][sel1])/np.sum(sel1)
        acc_nonUS = 100*np.sum(station_data['Class'][:,pp][sel2]==kg_class[ys,xs][sel2])/np.sum(sel2)
        print('Accuracy in US: '+str(acc_US)+' (n='+str(sum(sel1))+') Outside US: '+str(acc_nonUS)+' ('+str(sum(sel2))+')')
        
        sel3 = (valid) & (np.isnan(np.mean(station_data['Class'],axis=1))==False)
        acc_same = 100*np.sum(station_data['Class'][:,pp][sel3]==kg_class[ys,xs][sel3])/np.sum(sel3)    
        print('Accuracy same: '+str(acc_same))
        
        # Clear memory
        del kg_class
        del kg_major
        del kg_confidence
     
    # Save to csv
    print(df_accuracy)
    df_accuracy.set_index('Period',inplace=True)
    if os.path.isdir(os.path.join(config['folder_stats'],'validation'))==False:
        os.makedirs(os.path.join(config['folder_stats'],'validation'))
    df_accuracy.to_csv(os.path.join(config['folder_stats'],'validation','accuracy.csv'))
    
    pdb.set_trace()

    
if __name__ == '__main__':
    main()