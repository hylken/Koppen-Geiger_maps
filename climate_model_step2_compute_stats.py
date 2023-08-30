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
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.weightstats import ztest as ztest
import warnings
from skimage.transform import resize_local_mean

def main():


    #==============================================================================
    #   Settings
    #==============================================================================

    warnings.filterwarnings('ignore')
    np.set_printoptions(suppress=True)
    config = tools.load_config(sys.argv[1])


    #==============================================================================
    #   Preparation
    #==============================================================================

    DatesMon = pd.date_range(start=datetime(1850,1,1), end=datetime(2099,12,31), freq='MS')
    Years = np.unique(DatesMon.year)

    if os.path.isdir(os.path.join(config['folder_stats'],'sim_vs_obs'))==False:
        os.makedirs(os.path.join(config['folder_stats'],'sim_vs_obs'))
    if os.path.isdir(os.path.join(config['folder_stats'],'sim_vs_obs','sim'))==False:
        os.makedirs(os.path.join(config['folder_stats'],'sim_vs_obs','sim'))
    if os.path.isdir(os.path.join(config['folder_stats'],'sensitivity'))==False:
        os.makedirs(os.path.join(config['folder_stats'],'sensitivity'))

    
    #==============================================================================
    #   Load HadCRUT ensemble
    #==============================================================================

    if os.path.isfile(os.path.join(config['folder_stats'],'sim_vs_obs','HadCRUT.csv'))==False: 
        print('--------------------------------------------------------------------------------')
        print('Processing observed data')
        t = time.time()
        
        # Grid-cell area map
        res = 5
        lat = 90-np.arange(180/res)*res-res/2
        lon = -180+np.arange(360/res)*res+res/2 
        _, yi = np.meshgrid(lon, lat)
        area_map = (40075*res/360)**2*np.cos(np.deg2rad(yi))
        
        # Load HadCRUT monthly T data ensemble
        HadCRUT = np.zeros((len(Years),200),dtype=np.single)*np.NaN
        for ee in np.arange(1,201):
            filepath = glob.glo(os.path.join(config['folder_dataraw'], 'HadCRUT', 'HadCRUT.*.analysis.anomalies.'+str(ee)+'.nc'))[0]
            ncfile = Dataset(filepath)
            tmp = np.array(ncfile.variables['tas'][:]) 
            ncfile.close()
            tmp[tmp==-1e30] = np.NaN
            dates = pd.date_range(start=datetime(1850,1,1), end=datetime(2099,12,1), freq='MS')
            dates = dates[:tmp.shape[0]]
            HADCRU_T = np.zeros((tmp.shape[1],tmp.shape[2],len(DatesMon)),dtype=np.single)*np.NaN
            for ii in np.arange(len(DatesMon)):
                jj = (dates.year==DatesMon[ii].year) & (dates.month==DatesMon[ii].month)
                if np.sum(jj)>0:
                    map = np.flipud(np.squeeze(tmp[jj,:,:]))
                    map[np.isnan(map)] = np.nanmean(map)
                    HADCRU_T[:,:,ii] = map
            del tmp

            # Compute global means
            for yy in np.arange(len(Years)):
                sel = DatesMon.year==Years[yy]
                tmp = np.mean(HADCRU_T[:,:,sel],axis=2)
                HadCRUT[yy,ee-1] = np.mean(tmp*area_map)/np.mean(area_map) 
        
        HadCRUT = pd.DataFrame(HadCRUT,index=Years,columns=np.arange(1,201))
        HadCRUT.to_csv(os.path.join(config['folder_stats'],'sim_vs_obs','HadCRUT.csv'))
        print('Time elapsed is '+str(time.time() - t)+' sec')
    
    else:
        HadCRUT = pd.read_csv(os.path.join(config['folder_stats'],'sim_vs_obs','HadCRUT.csv'),index_col=0)
   
      
    #==============================================================================
    #   Estimate internal variability using CMIP6 and compare simulated and 
    #   observed trends
    #==============================================================================

    print('--------------------------------------------------------------------------------')
    print('Estimating internal variability and comparing simulated and observed trends')
    t = time.time()
    
    # Select trend period
    start_year = 1980
    end_year = 2014
            
    # Get list of models
    files = glob.glob(os.path.join(config['folder_out'], 'climate_model_data', 'historical_*_tas.npz'))
    models = np.unique([os.path.basename(x).split('_')[1] for x in files]).tolist()
    try:
        models.remove('historical')
    except:
        pass
    
    # Initialize results arrays
    pd_pvals = pd.DataFrame(index=models,columns=np.arange(1,201))
    pd_members = pd.DataFrame(index=models,columns=np.arange(20))
    pd_nmembers = pd.Series(index=models,dtype=object)
    pd_sim_trends = pd.DataFrame(index=models,columns=np.arange(20))
    pd_obs_trends = pd.Series(index=np.arange(1,201),dtype=object)
    pd_int_var = pd.Series(index=models,dtype=object)
    
    # Loop over models
    for mm in np.arange(len(models)):
        model = models[mm]
        
        #if model!='MPI-ESM1-2-HR':
        #    continue
        
        # Get list of model ensemble members
        files = glob.glob(os.path.join(config['folder_out'], 'climate_model_data', 'historical_'+model+'_*_tas.npz'))
        members = [os.path.basename(x).split('_')[2] for x in files]
        
        # Loop over model ensemble members
        sim_trends = np.zeros((len(members),))*np.NaN
        for ee in np.arange(np.min((len(members),20))):
            member = members[ee]
            pd_members.iloc[mm,ee] = member            
        
            # Load historical simulation            
            tas = np.zeros((len(Years),))*np.NaN
            with np.load(os.path.join(config['folder_out'], 'climate_model_data', 'historical_'+model+'_'+member+'_tas.npz')) as historical:
                data = historical['ts_mean_yr']
                yrs = historical['Years']
                for yy in np.arange(len(Years)):
                    sel = yrs==Years[yy]
                    if sum(sel)==1:
                        tas[yy] = data[sel]
                
            # Load SSP245 simulation (if it exists) to extend historical data
            # Not a good idea; almost 30% of the models has no SSP245 data
            #try:
            #    with np.load(os.path.join(config['folder_out'], 'climate_model_data', 'ssp245_'+model+'_'+member+'_tas.npz')) as future:
            #        for yy in np.arange(len(Years)):
            #            sel = future['Years']==Years[yy]
            #            if sum(sel)==1:
            #                tas[yy] = future['ts_mean_yr'][sel]
            #except:
            #    pass
                
            # Save to csv
            pd.Series(tas,index=Years.astype(int)).to_csv(os.path.join(config['folder_stats'],'sim_vs_obs','sim',model+'_'+member+'.csv'),header=False)
            
            # Simulated trend
            sel = (Years>=start_year) & (Years<=end_year) # & (~np.isnan(tas)) & (tas>-100) & (tas<100)
            Y, X = tas[sel], Years[sel]
            res = stats.linregress(X, Y)
            sim_trends[ee] = res.slope*10
            pd_sim_trends.iloc[mm,ee] = res.slope*10
            #plt.plot(Y)
        
        
        # Estimate internal variability
        # We need reasonable sample of trends to approximate the real distribution
        pd_nmembers.iloc[mm] = np.sum(np.isnan(sim_trends)==False)
        if pd_nmembers.iloc[mm]>=10:
            pd_int_var.iloc[mm] = np.nanstd(sim_trends)
            #plt.title(model)
            #plt.show()
            #plt.close()
        
        # Loop over observation ensemble members
        for ee in np.arange(200):
            sel = (HadCRUT.index>=start_year) & (HadCRUT.index<=end_year)
            Y, X = HadCRUT.iloc[sel,ee].values, HadCRUT.index[sel].values
            Y = Y+np.arange(len(Y))*0.0014 # Correction due to use of SST data over oceans (blending effect; see https://www.science.org/doi/10.1126/sciadv.aaz9549)            
            res = stats.linregress(X, Y)
            pd_obs_trends.iloc[ee] = res.slope*10
            
            # Compute p-value to test whether observed trend comes from 
            # distribution of simulated trends. The statsmodels ztest gives very
            # small p-values for some reason, so let's emulate the Matlab ztest 
            # function instead.
            #model_val['pvals'][mm,ee] = ztest(sim_trends, value=res.slope*10)[1] 
            z = (res.slope*10-np.nanmean(sim_trends))/(np.nanstd(sim_trends)/np.sqrt(1)) 
            pd_pvals.iloc[mm,ee] = stats.norm.cdf(z,0,1)*2
    
    # Output to csv files
    pd_pvals.to_csv(os.path.join(config['folder_stats'],'sim_vs_obs','pvals.csv'))
    pd_members.to_csv(os.path.join(config['folder_stats'],'sim_vs_obs','members.csv'))
    pd_nmembers.to_csv(os.path.join(config['folder_stats'],'sim_vs_obs','nmembers.csv'))
    pd_sim_trends.to_csv(os.path.join(config['folder_stats'],'sim_vs_obs','sim_trends.csv'))
    pd_obs_trends.to_csv(os.path.join(config['folder_stats'],'sim_vs_obs','obs_trends.csv'))
    pd_int_var.to_csv(os.path.join(config['folder_stats'],'sim_vs_obs','int_var.csv'))

    # Compute median pvals
    pvals_median = np.median(pd_pvals.values,axis=1)
    pd.Series(pvals_median,index=models).to_csv(os.path.join(config['folder_stats'],'sim_vs_obs','pvals_median.csv'))
    
    print('Time elapsed is '+str(time.time() - t)+' sec')
    
    
    #==============================================================================
    #   Compute Equilibrium Climate Sensitivity (ECS) and 
    #   Transient Climate Response (TCR)
    #==============================================================================

    print('--------------------------------------------------------------------------------')
    print('Computing ECS and TCR values')
    t = time.time()

    # Initialize results arrays
    pd_members = pd.DataFrame(index=models,columns=np.arange(20))
    pd_ECS = pd.DataFrame(index=models,columns=np.arange(20))
    pd_TCR = pd.DataFrame(index=models,columns=np.arange(20))    
    
    # Loop over models
    for mm in np.arange(len(models)):
        model = models[mm]
        
        # Get list of ensemble members
        files = glob.glob(os.path.join(config['folder_out'], 'climate_model_data', 'piControl_'+model+'_*_tas.npz'))
        members = [os.path.basename(x).split('_')[2] for x in files]
            
        #------------------------------------------------------------------------------
        #   Compute Equilibrium Climate Sensitivity (ECS)
        #   Approach: Gregory et al. (2004; https://doi.org/10.1029/2003GL018747)
        #   Note: We don't use rtmt (TOA net radiative flux) because it's available for 
        #   fewer models. We only use years 20 to 150 as the resulting ECS values agree
        #   better with slab models and long simulations (Dunne et al., 2020). Equation
        #   for computing regression line: 
        #   https://www4.stat.ncsu.edu/~dickey/summer_institute/formulas
        #------------------------------------------------------------------------------
        
        # Loop over ensemble members
        for ee in np.arange(len(members)):
            member = members[ee]
            pd_members.iloc[mm,ee] = member        
                    
            try:

                # Load data
                abrupt4xCO2_tas = np.load(os.path.join(config['folder_out'], 'climate_model_data', 'abrupt-4xCO2_'+model+'_'+member+'_tas.npz'))
                abrupt4xCO2_rsdt = np.load(os.path.join(config['folder_out'], 'climate_model_data', 'abrupt-4xCO2_'+model+'_'+member+'_rsdt.npz'))
                abrupt4xCO2_rsut = np.load(os.path.join(config['folder_out'], 'climate_model_data', 'abrupt-4xCO2_'+model+'_'+member+'_rsut.npz'))
                abrupt4xCO2_rlut = np.load(os.path.join(config['folder_out'], 'climate_model_data', 'abrupt-4xCO2_'+model+'_'+member+'_rlut.npz'))
                piControl_tas = np.load(os.path.join(config['folder_out'], 'climate_model_data', 'piControl_'+model+'_'+member+'_tas.npz'))
                piControl_rsdt = np.load(os.path.join(config['folder_out'], 'climate_model_data', 'piControl_'+model+'_'+member+'_rsdt.npz'))
                piControl_rsut = np.load(os.path.join(config['folder_out'], 'climate_model_data', 'piControl_'+model+'_'+member+'_rsut.npz'))
                piControl_rlut = np.load(os.path.join(config['folder_out'], 'climate_model_data', 'piControl_'+model+'_'+member+'_rlut.npz'))
                                 
                # Regression
                sel = np.arange(20,150)
                deltaN = (abrupt4xCO2_rsdt['ts_mean_yr'][sel]-abrupt4xCO2_rsut['ts_mean_yr'][sel]-abrupt4xCO2_rlut['ts_mean_yr'][sel]) \
                    -(piControl_rsdt['ts_mean_yr'][sel]-piControl_rsut['ts_mean_yr'][sel]-piControl_rlut['ts_mean_yr'][sel]) # Top of atmosphere (or model) net radiative flux anomaly (W/m2)
                deltaT = abrupt4xCO2_tas['ts_mean_yr'][sel]-piControl_tas['ts_mean_yr'][sel] # Temperature anomaly (degrees C)
                Y, X = deltaN, deltaT
                X[(X<-100) | (X>100)] = np.NaN            
                slope = np.sum((Y-np.mean(Y))*(X-np.mean(X)))/np.sum((X-np.mean(X))**2)
                intercept = np.mean(Y)-slope*np.mean(X)
                
                if sum(np.isnan(X))>10:
                    print(model+' '+member+' gaps in time series')
                
                # Compute ECS
                pd_ECS.iloc[mm,ee] = (-intercept/slope)/2
                
            except:
                pass
        
        # Try to compute ECS for models without both experiments for any member
        if np.isnan(np.nanmean(pd_ECS.iloc[mm,:])):
            files_abrupt4xCO2_tas = glob.glob(os.path.join(config['folder_out'], 'climate_model_data', 'abrupt-4xCO2_'+model+'_*_tas.npz'))
            files_abrupt4xCO2_rsdt = glob.glob(os.path.join(config['folder_out'], 'climate_model_data', 'abrupt-4xCO2_'+model+'_*_rsdt.npz'))
            files_abrupt4xCO2_rsut = glob.glob(os.path.join(config['folder_out'], 'climate_model_data', 'abrupt-4xCO2_'+model+'_*_rsut.npz'))
            files_abrupt4xCO2_rlut = glob.glob(os.path.join(config['folder_out'], 'climate_model_data', 'abrupt-4xCO2_'+model+'_*_rlut.npz'))
            files_piControl_tas = glob.glob(os.path.join(config['folder_out'], 'climate_model_data', 'piControl_'+model+'_*_tas.npz'))
            files_piControl_rsdt = glob.glob(os.path.join(config['folder_out'], 'climate_model_data', 'piControl_'+model+'_*_rsdt.npz'))
            files_piControl_rsut = glob.glob(os.path.join(config['folder_out'], 'climate_model_data', 'piControl_'+model+'_*_rsut.npz'))
            files_piControl_rlut = glob.glob(os.path.join(config['folder_out'], 'climate_model_data', 'piControl_'+model+'_*_rlut.npz'))
            try:
                abrupt4xCO2_tas = np.load(files_abrupt4xCO2_tas[0])
                abrupt4xCO2_rsdt = np.load(files_abrupt4xCO2_rsdt[0])
                abrupt4xCO2_rsut = np.load(files_abrupt4xCO2_rsut[0])
                abrupt4xCO2_rlut = np.load(files_abrupt4xCO2_rlut[0])
                piControl_tas = np.load(files_piControl_tas[0])
                piControl_rsdt = np.load(files_piControl_rsdt[0])
                piControl_rsut = np.load(files_piControl_rsut[0])
                piControl_rlut = np.load(files_piControl_rlut[0])                        
            
                # Regression
                sel = np.arange(20,150)
                deltaN = (abrupt4xCO2_rsdt['ts_mean_yr'][sel]-abrupt4xCO2_rsut['ts_mean_yr'][sel]-abrupt4xCO2_rlut['ts_mean_yr'][sel])-(piControl_rsdt['ts_mean_yr'][sel]-piControl_rsut['ts_mean_yr'][sel]-piControl_rlut['ts_mean_yr'][sel]) # Top of atmosphere (or model) net radiative flux anomaly (W/m2)
                deltaT = abrupt4xCO2_tas['ts_mean_yr'][sel]-piControl_tas['ts_mean_yr'][sel] # Temperature anomaly (degrees C)
                Y, X = deltaN, deltaT
                X[(X<-100) | (X>100)] = np.NaN            
                slope = np.sum((Y-np.mean(Y))*(X-np.mean(X)))/np.sum((X-np.mean(X))**2)
                intercept = np.mean(Y)-slope*np.mean(X)
                
                # Compute ECS
                pd_ECS.iloc[mm,0] = (-intercept/slope)/2
                        
            except:
                pass


        #------------------------------------------------------------------------------
        #   Compute Transient Climate Response (TCR)
        #   Approach: IPCC AR5 (2014)
        #   Note: Lauer et al. (2020; https://doi.org/10.5194/gmd-13-4205-2020) uses
        #   regression over control to calculate TCR, to account for residual model 
        #   drift. However, most studies don't seem to do this. And wouldn't model 
        #   drift also be present in onepctCO2?
        #------------------------------------------------------------------------------

        # Loop over ensemble members
        for ee in np.arange(len(members)):
            member = members[ee]
            pd_members.iloc[mm,ee] = member
                    
            try:
            
                # Load data
                with np.load(os.path.join(config['folder_out'], 'climate_model_data', '1pctCO2_'+model+'_'+member+'_tas.npz')) as onepctCO2_tas, \
                    np.load(os.path.join(config['folder_out'], 'climate_model_data', 'piControl_'+model+'_'+member+'_tas.npz')) as piControl_tas:
                    
                    # Regression
                    #Y = piControl_tas['ts_mean_yr'][1:140]
                    #Y[(Y<-100) | (Y>100)] = np.NaN
                    #X = piControl_tas['Years'][1:140]            
                    #slope = np.sum((Y-np.mean(Y))*(X-np.mean(X)))/np.sum((X-np.mean(X))**2)
                    #intercept = np.mean(Y)-slope*np.mean(X)
                    
                    # Compute TCR
                    #anomT = onepctCO2_tas['ts_mean_yr'][1:150]-(slope*piControl_tas['Years'][1:150]+intercept) # Temperature anomaly (degrees C)
                    #anomT = onepctCO2_tas['ts_mean_yr'][1:150]-(slope*piControl_tas['Years'][1:150]+intercept) # Temperature anomaly (degrees C)
                    #anomT[(anomT<-100) | (anomT>100)] = np.NaN
                    #pd_TCR.iloc[mm,ee] = np.mean(anomT[60:80])
                    
                    # Compute TCR
                    anomT = onepctCO2_tas['ts_mean_yr'][60:80]-piControl_tas['ts_mean_yr'][60:80]
                    anomT[(anomT<-100) | (anomT>100)] = np.NaN
                    pd_TCR.iloc[mm,ee] = np.mean(anomT)
                    
            except:
                pass
                        
        # Try to compute TCR for models without both experiments for any member (for example EC-Earth3)
        if np.isnan(np.nanmean(pd_TCR.iloc[mm,:])):
            files_1pctCO2 = glob.glob(os.path.join(config['folder_out'], 'climate_model_data', '1pctCO2_'+model+'_*_tas.npz'))
            files_piControl = glob.glob(os.path.join(config['folder_out'], 'climate_model_data', 'piControl_'+model+'_*_tas.npz'))
            try:
                onepctCO2_tas = np.load(files_1pctCO2[0])
                piControl_tas = np.load(files_piControl[0])
                
                # Compute TCR
                anomT = onepctCO2_tas['ts_mean_yr'][60:80]-piControl_tas['ts_mean_yr'][60:80]
                anomT[(anomT<-100) | (anomT>100)] = np.NaN
                pd_TCR.iloc[mm,0] = np.mean(anomT)
                    
            except:
                pass
        
    # Output ECS and TCR values to csv
    pd_members.to_csv(os.path.join(config['folder_stats'],'sensitivity','members.csv'))
    pd_ECS.to_csv(os.path.join(config['folder_stats'],'sensitivity','ECS.csv'))
    pd_TCR.to_csv(os.path.join(config['folder_stats'],'sensitivity','TCR.csv'))

    print('Time elapsed is '+str(time.time() - t)+' sec')
        

    #==============================================================================
    #   Screen models using historical trend and TCR. If TCR value not available 
    #   for a model, only use historical trend. Also create .txt file summarizing
    #   the sensitivities for the paper.
    #==============================================================================
    
    # Backup original stdout
    original_stdout = sys.stdout 
    
    with open(os.path.join(config['folder_stats'],'sensitivity','model_sensitivity_summary.txt'), 'w') as f:
        sys.stdout = f  # Change the standard output to the file we created.

        mean, std = np.mean(pd_obs_trends.values), np.std(pd_obs_trends.values)
        int_var = np.nanmean(pd_int_var.values)
        int_var_std = np.nanstd(pd_int_var.values)
        lo, hi = mean-np.sqrt(std**2+2*int_var**2), mean+np.sqrt(std**2+2*int_var**2)
        print('Observational trend (corrected for blending bias): '+str(mean)+' degrees C/decade')
        print('Observational trend uncertainty: '+str(std)+' degrees C/decade')
        print('Internal variability: '+str(int_var)+' degrees C/decade')
        print('Internal variability std: '+str(int_var_std)+' degrees C/decade')
        print('Observational trend range: '+str(lo)+' to '+str(hi)+' degrees C/decade')
        
        trend_vals = np.nanmean(pd_sim_trends,axis=1).astype(np.single)
        tcr_vals = np.nanmean(np.array(pd_TCR.values).astype(np.single),axis=1)
        ecs_vals = np.nanmean(np.array(pd_ECS.values).astype(np.single),axis=1)
        
        trend_within_range = ((trend_vals>=lo) & (trend_vals<=hi)).astype(np.single)
        trend_within_range[np.isnan(trend_vals)] = np.NaN
        tcr_within_range = ((tcr_vals>=1.4) & (tcr_vals<=2.2)).astype(np.single)
        tcr_within_range[np.isnan(tcr_vals)] = np.NaN
        ecs_within_range = ((ecs_vals>=2.5) & (ecs_vals<=4)).astype(np.single)
        ecs_within_range[np.isnan(ecs_vals)] = np.NaN
        
        model_subset = (trend_within_range==1) | (tcr_within_range==1)
        print('Number of models in Model Subset: '+str(sum(model_subset)))
        print('Included models:')
        print(*np.array(models)[model_subset], sep =', ')
        
        print('Number of models with air temperature trends: '+str(np.sum(np.isnan(trend_vals)==False)))
        print('Number of models with air temperature trends within the likely range: '+str(np.nansum(trend_within_range)))
        print('Number of models with TCR values: '+str(np.sum(np.isnan(tcr_vals)==False)))
        print('Number of models with TCR values within the likely range: '+str(np.nansum(tcr_within_range)))
        print('Number of models with ECS values: '+str(np.sum(np.isnan(ecs_vals)==False)))
        print('Number of models with ECS values within the likely range: '+str(np.nansum(ecs_within_range)))
        
        all_within_range = ((trend_within_range==1) | (np.isnan(trend_vals))) & ((tcr_within_range==1) | (np.isnan(tcr_vals))) & ((ecs_within_range==1) | (np.isnan(ecs_vals)))
        print('Number of models with all metrics within the likely range: '+str(np.sum(all_within_range)))

        all_three_available = (np.isnan(trend_vals)==False) & (np.isnan(tcr_vals)==False) & (np.isnan(ecs_vals)==False)
        print('Number of models with values for all three metrics: '+str(np.sum(all_three_available)))
        
        one_outside_range = (trend_within_range==0) | (tcr_within_range==0) | (ecs_within_range==0)
        print('Number of models with at least one metric outside of the likely range: '+str(np.sum(one_outside_range)))

        both_available = (np.isnan(trend_vals)==False) & (np.isnan(tcr_vals)==False)
        print('Number of models with air temperature trends and TCR values: '+str(np.sum(both_available)))
        print('Number of models with air temperature trends and TCR values and both outside the likely range: '+str(np.sum((both_available) & (trend_within_range==0) & (tcr_within_range==0))))
        
        only_trend_available = (np.isnan(trend_vals)==False) & (np.isnan(tcr_vals))
        print('Number of models with only air temperature trends: '+str(np.sum(only_trend_available)))
        print('Number of models with only air temperature trends and outside the likely range: '+str(np.sum((only_trend_available) & (trend_within_range==0))))
        
    sys.stdout = original_stdout  # Reset the standard output to its original value
    
    
    #==============================================================================
    #   Generate dataframe with projected mean changes in temperature and 
    #   precipitation for the land surface with and without model screening
    #==============================================================================
    
    print('--------------------------------------------------------------------------------')
    print('Generating dataframe with temperature and precipitation projections for the land surface')
    t = time.time()

    # Scenarios and variables to process
    scenarios = ['ssp119','ssp126','ssp245','ssp370','ssp434','ssp460','ssp585']
    reference_period = (1991,2020)
    future_period = (2071,2099)
    variables = ['tas','pr']
    statistics = ['ts_mean_yr','ts_min_yr','ts_max_yr','ts_mean_yr_land','ts_min_yr_land','ts_max_yr_land']

    # Loop over variables, statistics, scenarios, and models
    df = pd.DataFrame()
    for vv in np.arange(len(variables)):        
        variable = variables[vv]
        for aa in np.arange(len(statistics)):
            statistic = statistics[aa]
            for ss in np.arange(len(scenarios)):
                scenario = scenarios[ss]
                for mm in np.arange(len(models)):
                    model = models[mm]
                    
                    # Get list of model ensemble members
                    files = glob.glob(os.path.join(config['folder_out'], 'climate_model_data', scenario+'_'+model+'_*_'+variable+'.npz'))
                    members = [os.path.basename(x).split('_')[2] for x in files]
                    if len(members)==0:
                        continue
                    
                    # Loop over model ensemble members
                    sim_data = np.zeros((20,len(Years)),dtype=np.single)*np.NaN
                    for ee in np.arange(len(members)):
                        member = members[ee]            
                        
                        # Load historical simulation
                        try:
                            with np.load(os.path.join(config['folder_out'], 'climate_model_data', 'historical_'+model+'_'+member+'_'+variable+'.npz')) as historical:
                                yrs = historical['Years']
                                data = historical[statistic]
                                for yy in np.arange(len(Years)):
                                    sel = yrs==Years[yy]
                                    if sum(sel)==1:
                                        sim_data[ee,yy] = data[sel]
                        except:
                            continue
                            
                        # Load projection
                        with np.load(os.path.join(config['folder_out'], 'climate_model_data', scenario+'_'+model+'_'+member+'_'+variable+'.npz')) as future:
                            yrs = future['Years']
                            data = future[statistic]
                            for yy in np.arange(len(Years)):
                                sel = yrs==Years[yy]
                                if sum(sel)==1:
                                    sim_data[ee,yy] = data[sel]
                
                    # Mean over all ensemble members
                    sim_data = np.nanmean(sim_data,axis=0) 
                
                    # Compute difference between future and reference periods
                    sel_ref = (Years>=reference_period[0]) & (Years<=reference_period[1])
                    sel_fut = (Years>=future_period[0]) & (Years<=future_period[1])                
                    if variable=='tas':
                        change = np.mean(sim_data[sel_fut])-np.mean(sim_data[sel_ref])
                    elif variable=='pr':
                        change = 100*np.mean(sim_data[sel_fut])/np.mean(sim_data[sel_ref])-100
                        
                    # Append to dataframe
                    if ~np.isnan(change):
                        newrow = {'Variable':variable,'Statistic':statistic,'Scenario':scenario,'Model':model,'Strategy':'All Models','Value':change}
                        df = pd.concat([df, pd.DataFrame([newrow])], ignore_index=True)
                    
                    # Append to dataframe
                    if (~np.isnan(change)) & (model_subset[mm]==True):
                        newrow = {'Variable':variable,'Statistic':statistic,'Scenario':scenario,'Model':model,'Strategy':'Model Subset','Value':change}
                        df = pd.concat([df, pd.DataFrame([newrow])], ignore_index=True)
          
    # Output to csv file
    df.to_csv(os.path.join(config['folder_stats'],'projected_change.csv'))

    print('Time elapsed is '+str(time.time() - t)+' sec')
    
    
    #==============================================================================
    #   Generate maps of projected global mean changes in temperature and 
    #   precipitation
    #==============================================================================
    
    print('--------------------------------------------------------------------------------')
    print('Generating maps of temperature and precipitation projections')
    
    statistics = ['data_mean_yr','data_min_yr','data_max_yr']
    
    change = np.zeros((len(variables),len(statistics),len(scenarios),len(models),180,360),dtype=np.single)*np.NaN
    for vv in np.arange(len(variables)):
        variable = variables[vv]
        for aa in np.arange(len(statistics)):
            statistic = statistics[aa]
            for ss in np.arange(len(scenarios)):
                scenario = scenarios[ss]
                
                if scenario!='ssp245':
                    continue
                    
                print(variable+' '+statistic+' '+scenario)
                
                t = time.time()
                
                for mm in np.arange(len(models)):
                    model = models[mm]
                   
                    # Get list of model ensemble members
                    files = glob.glob(os.path.join(config['folder_out'], 'climate_model_data', scenario+'_'+model+'_*_'+variable+'.npz'))
                    members = [os.path.basename(x).split('_')[2] for x in files]
                    if len(members)==0:
                        continue
                    
                    # Loop over model ensemble members
                    with np.load(os.path.join(config['folder_out'], 'climate_model_data', scenario+'_'+model+'_'+members[0]+'_'+variable+'.npz')) as future:
                        dims = future[statistic].shape
                    sim_data = np.zeros((dims[0],dims[1],20,len(Years)),dtype=np.single)*np.NaN
                    for ee in np.arange(len(members)):
                        member = members[ee]            
                        
                        # Load simulation data
                        try:
                            with np.load(os.path.join(config['folder_out'], 'climate_model_data', 'historical_'+model+'_'+member+'_'+variable+'.npz')) as historical:
                                data = historical[statistic]
                                yrs = historical['Years']
                                for yy in np.arange(len(Years)):
                                    sel = yrs==Years[yy]
                                    if sum(sel)==1:
                                        sim_data[:,:,ee,yy] = np.squeeze(data[:,:,sel])
                                        
                            with np.load(os.path.join(config['folder_out'], 'climate_model_data', scenario+'_'+model+'_'+member+'_'+variable+'.npz')) as future:
                                data = future[statistic]
                                yrs = future['Years']
                                for yy in np.arange(len(Years)):
                                    sel = yrs==Years[yy]
                                    if sum(sel)==1:
                                        sim_data[:,:,ee,yy] = np.squeeze(data[:,:,sel])
                        except:
                            continue
                        
                    # Mean over all ensemble members
                    sim_data = np.nanmean(sim_data,axis=2) 
                    
                    # Compute difference between future and reference periods
                    sel_ref = (Years>=reference_period[0]) & (Years<=reference_period[1])
                    sel_fut = (Years>=future_period[0]) & (Years<=future_period[1])                
                    if variable=='tas':
                        change[vv,aa,ss,mm,:,:] = resize_local_mean(np.nanmean(sim_data[:,:,sel_fut],axis=2)-np.nanmean(sim_data[:,:,sel_ref],axis=2),(180,360))
                    elif variable=='pr':
                        change[vv,aa,ss,mm,:,:] = resize_local_mean(100*np.nanmean(sim_data[:,:,sel_fut],axis=2)/np.nanmean(sim_data[:,:,sel_ref],axis=2)-100,(180,360))
                print('Time elapsed is '+str(time.time() - t)+' sec')                
            
    np.savez_compressed(os.path.join(config['folder_stats'],'projected_change.npz'),\
        change=change,variables=variables,statistics=statistics,scenarios=scenarios,models=models,\
        model_subset=model_subset)

    pdb.set_trace()
    
    
if __name__ == '__main__':
    main()