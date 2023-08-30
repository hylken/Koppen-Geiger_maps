#!/usr/bin/env python  
# -*- coding: utf-8 -*-

__author__ = "Hylke E. Beck"
__email__ = "hylke.beck@gmail.com"
__date__ = "November 2022"

import os
import sys
import pdb
import glob
import pandas as pd
import numpy as np
import tools
import matplotlib.pyplot as plt
from adjustText import adjust_text
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.patches import Rectangle
import seaborn as sns
import warnings
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
import scipy

def main():


    #==============================================================================
    #   Settings
    #==============================================================================

    config = tools.load_config(sys.argv[1])
    dirout = os.path.join(sys.path[0], config['folder_stats'])

    warnings.filterwarnings('ignore')
    np.set_printoptions(suppress=True)
    
    # Download font from https://fontsgeek.com/fonts/Myriad-Pro-Regular
    # Place in ~/.fonts
    # Rebuild cache with: fc-cache -f -v
    # Verify installation with: fc-list | egrep -i "myriad"
    # Delete matplotlib cache (Ubuntu only?): rm ~/.cache/matplotlib -rf
    plt.rcParams["font.family"] = "Myriad Pro"


    #==============================================================================
    #   Generate historical trend vs TCR scatterplot
    #==============================================================================
    
    print('Generating historical trend vs TCR scatterplot')
    
    # Color theme of figures: 
    # https://paletton.com/#uid=71x1q1kllllaFw0g0qFqFg0w0aFkllllaFw0g0qFqFg0w0aFkllllaFw0g0qFqFg0w0aFkllllaFw0g0qFqFg0w0aF

    # Load data
    pd_pvals_median = pd.read_csv(os.path.join(dirout,'sim_vs_obs','pvals_median.csv'),index_col=0)
    pd_nmembers = pd.read_csv(os.path.join(dirout,'sim_vs_obs','nmembers.csv'),index_col=0)
    pd_sim_trends = pd.read_csv(os.path.join(dirout,'sim_vs_obs','sim_trends.csv'),index_col=0)
    pd_obs_trends = pd.read_csv(os.path.join(dirout,'sim_vs_obs','obs_trends.csv'),index_col=0)
    pd_int_var = pd.read_csv(os.path.join(config['folder_stats'],'sim_vs_obs','int_var.csv'),index_col=0)
    pd_ECS = pd.read_csv(os.path.join(dirout,'sensitivity','ECS.csv'),index_col=0)
    pd_TCR = pd.read_csv(os.path.join(dirout,'sensitivity','TCR.csv'),index_col=0)
    dset_projected_change = np.load(os.path.join(config['folder_stats'],'projected_change.npz'))
    
    models = np.array(pd_sim_trends.index.tolist())

    fig, ax = plt.subplots()
    
    # Plot TCR range and best estimate and historical trend range and best estimate
    # Historical trend range includes internal variability 
    mean, std = np.mean(pd_obs_trends.values), np.std(pd_obs_trends.values)
    int_var = np.nanmean(pd_int_var.values)
    lo, hi = mean-std-np.sqrt(2*int_var**2), mean+std+np.sqrt(2*int_var**2)
    ax.axvspan(lo, hi, color='#807115', alpha=0.1, lw=0) 
    ax.plot(np.array([mean,mean]), np.array([0,5]), color='#807115',lw=1, alpha=0.8,ls=':')
    ax.axhspan(1.4, 2.2, color='#2B1657', alpha=0.1, lw=0) # AR6 TCR likely range 
    ax.plot(np.array([0,5]),np.array([1.8,1.8]),color='#2B1657',lw=1, alpha=0.7,ls=':')
    ax.set_xlabel('Air temperature trend 1980–2014 (°C decade$^{-1}$)')
    ax.set_ylabel('Transient Climate Response (TCR; °C)')

    # Colormap for ECS values
    colors=["#297C46","#AA9B39","#AA4A39"]
    nodes = [0.0,0.5,1.0]
    my_cmap = LinearSegmentedColormap.from_list("mycmap", list(zip(nodes, colors)))
    
    # Generate scatterplot
    ECS_min, ECS_max = 2, 6
    x, y, c = np.nanmean(pd_sim_trends,axis=1), np.nanmean(pd_TCR,axis=1), (np.nanmean(pd_ECS,axis=1)-ECS_min)/(ECS_max-ECS_min)
    c = c.clip(0,1)
    sel = ~np.isnan(x+y+c)
    x, y, c, l = x[sel], y[sel], c[sel], models[sel]
    ax.scatter(x,y,s=13,c=c,vmin=0,vmax=1,cmap=my_cmap,zorder=9999)
    
    # Add grid lines
    ax.grid(True)
    ax.grid(color=(0.95,0.95,0.95))
    ax.set_axisbelow(True)    
    
    # Set figure layout and plot model names
    plt.xlim([0.1, 0.4])
    plt.ylim([1.25, 3.25])    
    x0,x1 = ax.get_xlim()
    y0,y1 = ax.get_ylim()
    ax.set_aspect(abs(x1-x0)/abs(y1-y0))
    fig.set_size_inches(4.5, 4.5)
    ax.set_position([0.13,0.13,0.83,0.83]) 
    texts = [plt.text(x[i], y[i],l[i],size=7,color='k',alpha=0.3) for i in range(len(x))]
    adjust_text(texts, x, y, arrowprops=dict(arrowstyle="-", lw=0.3,alpha=0.3,color='k'),
        expand_objects=(1.2,1.2),expand_text=(1.2,1.2),lim=600000,force_text=0.05,force_objects=0.05)
    
    # Set figure layout and print    
    fig.savefig(os.path.join(dirout,'scatter_sim_trend_vs_TCR.png'), dpi=300) #,bbox_inches='tight'
    plt.close('all')
    
    # Separate colorbar figure
    a = np.array([[ECS_min, ECS_max]])
    plt.figure(figsize=(4.5, 0.7))
    img = plt.imshow(a, cmap=my_cmap)
    plt.gca().set_visible(False)
    cax = plt.axes([0.3, 0.7, 0.5, 0.25]) 
    cb = plt.colorbar(orientation="horizontal", cax=cax)
    cb.set_label('Equilibrium Climate Sensitivity (ECS; °C)')
    plt.savefig(os.path.join(dirout,'scatter_sim_trend_vs_TCR_colorbar.png'), dpi=300) #,bbox_inches='tight'
    plt.close('all')
    
    
    #==============================================================================
    #   Generate .tex table with model statistics
    #==============================================================================
    
    print('Generating table with model statistics')
    
    # Generate tex file
    with open(os.path.join(dirout,'model_values.tex'), 'w') as f:
        for mm in np.arange(len(models)):
            f.write(models[mm]+' & '+str(pd_nmembers.iloc[mm,0])+' & ')
            f.write('$'+"{:.3f}".format(np.nanmean(pd_sim_trends.iloc[mm,:]))+'$ & ') 
            f.write('$'+"{:.3f}".format(pd_int_var.iloc[mm,0])+'$ & ') 
            f.write('$'+"{:.2f}".format(np.nanmean(pd_TCR.iloc[mm,:]))+'$ & ') 
            f.write('$'+"{:.2f}".format(np.nanmean(pd_ECS.iloc[mm,:]))+'$ & ')
            included = 'No'
            if dset_projected_change['model_subset'][mm]==True:
                included = 'Yes'
            f.write(included+'\\\\\n')
        
    # Replace missing values in tex file with dash
    with open(os.path.join(dirout,'model_values.tex'), 'r') as file:
        filedata = file.read()
    filedata = filedata.replace('$nan$', '--')
    with open(os.path.join(dirout,'model_values.tex'), 'w') as file:
        file.write(filedata)


    #==============================================================================
    #   Print mean global warming anomaly in future for each scenario
    #==============================================================================
    
    # Load future change estimates
    df = pd.read_csv(os.path.join(dirout,'projected_change.csv'),index_col=0)
    scenarios = np.unique(df['Scenario']).tolist()
    
    # Load HadCRUT data and compute difference between industrial and reference 
    HadCRUT = pd.read_csv(os.path.join(config['folder_stats'],'sim_vs_obs','HadCRUT.csv'),index_col=0)
    HadCRUT_mean = np.mean(HadCRUT,axis=1)
    reference_period1 = (1850,1900)
    reference_period2 = (1961,1990)
    reference_period3 = (1991,2020)
    sel1 = (HadCRUT.index>=reference_period1[0]) & (HadCRUT.index<=reference_period1[1])
    sel2 = (HadCRUT.index>=reference_period2[0]) & (HadCRUT.index<=reference_period2[1])
    sel3 = (HadCRUT.index>=reference_period3[0]) & (HadCRUT.index<=reference_period3[1])
    anomaly = np.mean(HadCRUT_mean[sel3])-np.mean(HadCRUT_mean[sel1])
    
    # Compute future anomaly
    for ss in np.arange(len(scenarios)):
        scenario = scenarios[ss]
        sel = (df['Variable']=='tas') & (df['Statistic']=='ts_mean_yr') & (df['Scenario']==scenario) & (df['Strategy']=='Model Subset')
        change = anomaly+np.mean(df['Value'][sel])
        print(f'{change:.2f}'+' '+scenario)
    
    
    #==============================================================================
    #   Generate box plots showing warming for the land surface for different 
    #   scenarios, both using all models and the model subset
    #==============================================================================
    
    # Variables to plot
    variables = ['tas','pr']
    variable_labels = ['air temperature change (°C)','precipitation change (%)']
    statistics = ['ts_mean_yr_land','ts_min_yr_land','ts_max_yr_land']
    statistic_labels = ['Mean ','Annual monthly minimum\n','Annual monthly maximum\n']

    # Loop over variables
    for vv in np.arange(len(variables)):
        variable = variables[vv]
        for aa in np.arange(len(statistics)):
            statistic = statistics[aa]
            print('Generating box plots '+variable+' '+statistic)
            fig, ax = plt.subplots()
            
            # Plot boxes
            sns.set_palette(["#655091", "#4E9C68"]) 
            ax = sns.boxplot(x='Scenario',y='Value',data=df[(df["Variable"]==variable) & (df["Statistic"]==statistic)],\
                hue='Strategy', showfliers=False, linewidth=1, whis=[5,95], saturation=1, medianprops=dict(color='#000000',alpha=0.4))
            tools.adjust_box_widths(fig,0.8)
            
            # Overlay values
            sns.set_palette(["#8B7AAE", "#7DBB92"]) 
            sns.swarmplot(x='Scenario',y='Value',data=df[(df["Variable"]==variable) & (df["Statistic"]==statistic)],\
                hue='Strategy', alpha=0.4, dodge=True, ec='k', linewidth=0, size=3, ax=ax)

            # Add horizontal grid lines
            ax.yaxis.grid(True)
            ax.yaxis.grid(color=(0.9,0.9,0.9))
            ax.set_axisbelow(True)
            sns.despine(ax=ax)
            
            # Fix x-labels
            labels = ax.get_xticklabels()
            labels = [x.get_text().upper() for x in labels]
            labels = [x[:4]+'-'+x[4]+'.'+x[5:] for x in labels]
            ax.set_xticklabels(labels)            
            ax.set(xlabel='', ylabel=statistic_labels[aa]+variable_labels[vv])
            plt.xticks(rotation=22.5, ha='right')
            
            # Draw legend on top panel
            if (variable=='tas') & (statistic=='ts_min_yr_land'):                
                handles, labels = ax.get_legend_handles_labels()
                ax.legend(handles[:2], labels[:2], title='', loc='upper center',frameon=True,facecolor='w',edgecolor='w')
            else:
                plt.legend([],[], frameon=False)
            
            # Remove x-axis except on bottom panel
            if (variable!='pr'):
                ax.spines.bottom.set_visible(False)
                ax.axes.get_xaxis().set_visible(False)
            
            # Set figure layout and print
            ax.set_position([0.13,0.15,0.83,0.83])
            fig.set_size_inches(4.5, 2.6)
            fig.savefig(os.path.join(dirout,'boxplot_'+variable+'_'+statistic+'.png'), dpi=300) #,bbox_inches='tight'
            plt.close('all')

              
    #==============================================================================
    #   Generate figures of projected global mean changes in temperature and 
    #   precipitation and uncertainty
    #==============================================================================
    
    # Variables to plot
    vars = dset_projected_change['variables']    
    var_names = ['air temperature','precipitation']
    var_units = ['°C','%']
    stat_labels = ['Mean','Annual monthly minimum','Annual monthly maximum']
    stats = dset_projected_change['statistics']    
    
    ss = np.where(dset_projected_change['scenarios']=='ssp245')[0][0]
    scenario = dset_projected_change['scenarios'][ss]
    
    # Grid-cell area map
    res = 1
    lat = 90-np.arange(180/res)*res-res/2
    lon = -180+np.arange(360/res)*res+res/2 
    _, yi = np.meshgrid(lon, lat)
    area_map = (40075*res/360)**2*np.cos(np.deg2rad(yi))
    
    # Loop over variables
    for vv in np.arange(len(vars)):
        var = vars[vv]
        for aa in np.arange(len(stats)):
            stat = stats[aa]
            print('Making maps '+var+' '+stat)
            
            tmp = dset_projected_change['change'][vv,aa,ss,:,:,:]
            
            # Mean based on screened model subset
            var_mean_ranges = [(0,5),(-30,30)]
            if stat=='data_min_yr':
                var_mean_ranges = [(0,8),(-30,30)]
            var_mean_cmaps = [sns.color_palette("flare", as_cmap=True),sns.diverging_palette(145, 300, s=60, as_cmap=True)]
            fig, ax = plt.subplots()
            m = Basemap(projection='robin',lon_0=0,resolution='c')
            m.drawcoastlines(linewidth=0.5)
            im = m.imshow(np.flipud(np.nanmean(tmp[dset_projected_change['model_subset'],:,:],axis=0)),var_mean_cmaps[vv],vmin=var_mean_ranges[vv][0],vmax=var_mean_ranges[vv][1])
            m.drawparallels(np.arange(-90,120,30),linewidth=0.5)
            m.drawmeridians(np.arange(0,360,60),linewidth=0.5)
            cb = plt.colorbar(im,orientation='horizontal',fraction=0.046, pad=0.04)
            cb.set_label('Best estimate (mean across Model Subset; '+var_units[vv]+')')
            fig.set_size_inches(4.5, 4)
            fig.savefig(os.path.join(dirout,'map_'+scenario+'_'+var+'_'+stat+'_meansubset.png'), dpi=300,bbox_inches='tight')
            plt.close('all')
            
            # Std based on screened model subset
            var_std_ranges = [(0,2.5),(0,30)]
            if stat=='data_min_yr':
                var_std_ranges = [(0,4),(0,30)]
            var_std_cmaps = [sns.color_palette("crest", as_cmap=True),sns.color_palette("crest", as_cmap=True)]
            fig, ax = plt.subplots()
            m = Basemap(projection='robin',lon_0=0,resolution='c')
            m.drawcoastlines(linewidth=0.5)
            im = m.imshow(np.flipud(np.nanstd(tmp[dset_projected_change['model_subset'],:,:],axis=0)),var_std_cmaps[vv],vmin=var_std_ranges[vv][0],vmax=var_std_ranges[vv][1])
            m.drawparallels(np.arange(-90,120,30),linewidth=0.5)
            m.drawmeridians(np.arange(0,360,60),linewidth=0.5)
            cb = plt.colorbar(im,orientation='horizontal',fraction=0.046, pad=0.04)
            cb.set_label('Uncertainty (standard deviation across Model Subset; '+var_units[vv]+')')
            fig.set_size_inches(4.5, 4)
            fig.savefig(os.path.join(dirout,'map_'+scenario+'_'+var+'_'+stat+'_stdsubset.png'), dpi=300,bbox_inches='tight')
            plt.close('all')
            
            # Difference in mean between subset and all
            var_meandiff_ranges = [(-1,1),(-20,20)]
            var_meandiff_cmaps = [sns.color_palette("vlag", as_cmap=True),sns.diverging_palette(220, 20, as_cmap=True).reversed()]
            fig, ax = plt.subplots()
            m = Basemap(projection='robin',lon_0=0,resolution='c')
            m.drawcoastlines(linewidth=0.5)
            im = m.imshow(np.flipud(np.nanmean(tmp[dset_projected_change['model_subset'],:,:],axis=0)-np.nanmean(tmp,axis=0)),var_meandiff_cmaps[vv],vmin=var_meandiff_ranges[vv][0],vmax=var_meandiff_ranges[vv][1])
            m.drawparallels(np.arange(-90,120,30),linewidth=0.5)
            m.drawmeridians(np.arange(0,360,60),linewidth=0.5)
            cb = plt.colorbar(im,orientation='horizontal',fraction=0.046, pad=0.04)
            cb.set_label('Difference in best estimate between\nModel Subset and All Models ('+var_units[vv]+')')
            fig.set_size_inches(4.5, 4)
            fig.savefig(os.path.join(dirout,'map_'+scenario+'_'+var+'_'+stat+'_meandiff.png'), dpi=300,bbox_inches='tight')
            plt.close('all')
            
            # Difference in std between subset and all 
            var_stddiff_ranges = [(0.25,1.75),(0.25,1.75)]
            var_stddiff_cmaps = [sns.diverging_palette(150, 275, s=80, l=55, n=9, as_cmap=True),sns.diverging_palette(150, 275, s=80, l=55, n=9, as_cmap=True)]
            fig, ax = plt.subplots()
            m = Basemap(projection='robin',lon_0=0,resolution='c')
            m.drawcoastlines(linewidth=0.5)
            uncertainty_ratio = np.nanstd(tmp[dset_projected_change['model_subset'],:,:],axis=0)/np.nanstd(tmp,axis=0)
            im = m.imshow(np.flipud(uncertainty_ratio),var_stddiff_cmaps[vv],vmin=var_stddiff_ranges[vv][0],vmax=var_stddiff_ranges[vv][1])
            m.drawparallels(np.arange(-90,120,30),linewidth=0.5)
            m.drawmeridians(np.arange(0,360,60),linewidth=0.5)
            cb = plt.colorbar(im,orientation='horizontal',fraction=0.046, pad=0.04)
            cb.set_label('Ratio of Model Subset uncertainty\nto All Models uncertainty')
            fig.set_size_inches(4.5, 4)
            fig.savefig(os.path.join(dirout,'map_'+scenario+'_'+var+'_'+stat+'_stddiff.png'), dpi=300,bbox_inches='tight')
            plt.close('all')
    
            print(var+' '+stat+' mean uncertainty reduction '+str(np.nanmean(uncertainty_ratio*area_map)/np.nanmean(area_map)))
            
    # Experiment to determine average underestimation of standard deviation from small samples
    iqrs = np.zeros((10000,))*np.NaN
    for ii in np.arange(len(iqrs)): 
        iqrs[ii] = scipy.stats.iqr(np.random.normal(loc=0.0, scale=1.0, size=(30,)))        
    
    pdb.set_trace()
    
    
if __name__ == '__main__':
    main()