#!/usr/bin/env python  
# -*- coding: utf-8 -*-

__author__ = "Hylke E. Beck"
__email__ = "hylke.beck@gmail.com"
__date__ = "November 2022"

import os
import sys
import pdb
import time
import pandas as pd
import numpy as np
import tools
from datetime import datetime
from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
import warnings
import shapefile
import plotly.express as px
import plotly.graph_objects as go
        
def main():


    #==============================================================================
    #   Settings
    #==============================================================================
    
    np.set_printoptions(suppress=True)
    config = tools.load_config(sys.argv[1])
    plt.rcParams["font.family"] = "Myriad pro"
    
    koppen_table = pd.read_csv(os.path.join('assets','koppen_table.csv'))
    country_shp = shapefile.Reader(os.path.join('assets','TM_WORLD_BORDERS-0.3','TM_WORLD_BORDERS-0.3.shp'),encoding='ISO8859-1')
    treelines_shp = shapefile.Reader(os.path.join('assets','treeline','treeline_simple.shp'),encoding='ISO8859-1')

    rgb_vals = koppen_table[['Red','Green','Blue']].values
    kg_cmap = matplotlib.colors.ListedColormap(rgb_vals/255)
    kg_cmap.set_bad('white',1.)
    
    
    #==============================================================================
    #   Generate figures of Koppen-Geiger maps for periods and scenarios
    #==============================================================================
    
    regions = [
        ['World',(-180, 180, -60, 84),False,(10,5)],
        ['Alps',(6,14,42,50),True,(2.8,2.8)],
        ['Rocky Mts.',(-122.05,-112,44,54),True,(2.8,2.8)]
        ]
        
    # Create output folder if it doesn't exist
    if os.path.isdir(os.path.join(config['folder_stats'],'climatologies'))==False:
        os.makedirs(os.path.join(config['folder_stats'],'climatologies'))        
    
    # Loop over all files
    for root, dirs, files in os.walk(os.path.join(config['folder_out'],'climatologies')):
        for file in files:
            suffix = str(180/config['mapsize'][0]).replace('.','p')
            if file=='koppen_geiger_'+suffix+'.nc': 
                t0 = time.time()
                print('Producing figures of '+file)
                dset = Dataset(os.path.join(root,file))
                data = np.array(dset.variables['kg_class'][:]).astype(np.single)
                dset.close()
                data[data==0] = np.NaN
                for rr in np.arange(len(regions)):
                    fname = os.path.join(root,regions[rr][0]+'_'+file).replace('.nc','.png').replace(os.path.join(config['folder_out'],'climatologies'),'').replace(os.path.sep,'_')[1:]
                    figout = os.path.join(config['folder_stats'],'climatologies',fname)
                    tools.plot_map(
                        data=data,
                        data_extent=(-180, 180, -90, 90),
                        figout=figout,
                        figdims=regions[rr][3],
                        cmap=kg_cmap,
                        plot_extent=regions[rr][1],
                        lims=(0.5,30.5),
                        interpolation='nearest',
                        shp=country_shp,
                        color='k',
                        show_axes=regions[rr][2]
                        )
                print("Time elapsed is "+str(time.time()-t0)+" sec")
    
    pdb.set_trace()
    
    #==============================================================================
    #   Generate figures of Koppen-Geiger maps with tree lines overlay
    #==============================================================================
    '''
    regions = [
        ['Western Canada',(-135,-129,58,64),False,(2.8,2.8)], 
        ['Carpathian Mts',(24,26,44.5,46.5),False,(2.8,2.8)],
        ['Siberian Federal Zone',(108,118,52,62),False,(2.8,2.8)],
        ['Merida',(-74,-68,5,11),False,(2.8,2.8)],
        ['Gaoligong Mountain',(98,102,26,30),False,(2.8,2.8)],
        ['Rivero Island',(-74.5,-72.5,-45.5,-43.5),False,(2.8,2.8)],
        ['Groot Winterhoek Nature Reserve',(18,20,-34,-32),False,(2.8,2.8)],
        ['Snow Mts',(147.5,149.5,-37.5,-35.5),False,(2.8,2.8)],
        ]
        
    # Create output folder if it doesn't exist
    if os.path.isdir(os.path.join(config['folder_stats'],'climatologies'))==False:
        os.makedirs(os.path.join(config['folder_stats'],'climatologies'))        
    
    # Loop over regions
    suffix = str(180/config['mapsize'][0]).replace('.','p')
    filepath = os.path.join(config['folder_out'],'climatologies',
        str(config['periods_historical'][-1][0])+'_'+str(config['periods_historical'][-1][1]),
        'koppen_geiger_'+suffix+'.nc')
    print('Loading '+filepath)
    t0 = time.time()
    dset = Dataset(filepath)
    data = np.array(dset.variables['kg_class'][:]).astype(np.single)
    dset.close()
    data[data==0] = np.NaN
    for rr in np.arange(len(regions)):
        print('Making figure for '+regions[rr][0])
        fname = regions[rr][0].replace(' ','_')+'_treelines.png'
        figout = os.path.join(config['folder_stats'],'climatologies',fname)
        tools.plot_map(
            data=data,
            data_extent=(-180, 180, -90, 90),
            figout=figout,
            figdims=regions[rr][3],
            cmap=kg_cmap,
            plot_extent=regions[rr][1],
            lims=(0.5,30.5),
            interpolation='nearest',
            shp=treelines_shp,
            color='r',
            show_axes=regions[rr][2]
            )
    print("Time elapsed is "+str(time.time()-t0)+" sec")
    '''
    

    #==============================================================================
    #   Generate LaTeX table of classification accuracy
    #==============================================================================
    '''
    df_accuracy = pd.read_csv(os.path.join(config['folder_stats'],'validation','accuracy.csv'),index_col=0)
    
    with open(os.path.join(config['folder_stats'],'validation','accuracy.tex'), 'w') as f:
        for pp in np.arange(df_accuracy.shape[0]):
            f.write(df_accuracy.index[pp].replace('(','').replace(')','').replace(', ','--')+' & ')
            f.write('$'+"{:.0f}".format(df_accuracy.iloc[pp,0])+'$ & ')
            f.write('$'+"{:.1f}".format(df_accuracy.iloc[pp,1])+'$ & ')
            f.write('$'+"{:.1f}".format(df_accuracy.iloc[pp,2])+'$ && ')
            f.write('$'+"{:.1f}".format(df_accuracy.iloc[pp,3])+'$ & ')
            f.write('$'+"{:.1f}".format(df_accuracy.iloc[pp,4])+'$')
            if pp<df_accuracy.shape[0]-1:
                f.write('\\\\\n')
    '''

    #==============================================================================
    #   Generate Sankey diagram
    #==============================================================================
    
    scenarios = ['ssp119','ssp126','ssp245','ssp370','ssp434','ssp460','ssp585']
    periods = config['periods_historical']+config['periods_future']
    trans_thresh_mm2 = 0.6
    
    colors = np.array([ \
        'rgba(50, 50, 255, opacity)', \
        'rgba(255, 0, 0, opacity)', \
        'rgba(90, 220, 70, opacity)', \
        'rgba(50, 150, 150, opacity)', \
        'rgba(100, 100, 100, opacity)' \
        ])
    
    for scenario in scenarios:
    
        # Load area and transition data
        df_kg_major_area_pct = pd.read_csv(os.path.join(config['folder_stats'],'climatologies',scenario+'_kg_major_area_pct.csv'),index_col=0)
        df_kg_major_area_mm2 = pd.read_csv(os.path.join(config['folder_stats'],'climatologies',scenario+'_kg_major_area_mm2.csv'),index_col=0)
        df_transitions_mm2 = pd.read_csv(os.path.join(config['folder_stats'],'climatologies',scenario+'_transitions_mm2.csv'),index_col=None)
        if scenario=='ssp245':
            print('Pausing')
            pdb.set_trace()
            
        # Make node and link lists
        node_colors, node_areas_mm2, node_areas_pct = [], [], []
        link_sources, link_targets = [], []
        for ll in np.arange(len(periods)):
        
            # Node data
            node_colors = node_colors+[c.replace('opacity', '1') for c in colors]
            node_areas_mm2 = node_areas_mm2+df_kg_major_area_mm2.loc[str(periods[ll]),:].tolist()
            node_areas_pct = node_areas_pct+df_kg_major_area_pct.loc[str(periods[ll]),:].tolist()
            
            # Link sources and targets            
            if ll<len(periods)-1:
                sel = np.array(df_transitions_mm2['Area']>trans_thresh_mm2) & np.array(df_transitions_mm2['From']==str(periods[ll]))
                link_sources = link_sources+(df_transitions_mm2['Source'][sel]-1+ll*5).tolist()
                link_targets = link_targets+(df_transitions_mm2['Target'][sel]-1+(ll+1)*5).tolist()
                
        # Select links to plot
        sel = np.array(df_transitions_mm2['Area']>trans_thresh_mm2)
        
        # Link colors
        inds = df_transitions_mm2['Source'][sel].values.astype(int)-1
        link_colors = colors[inds].tolist()
        opacities = (1/(df_transitions_mm2['Area'].values[sel]+0)*0.85+0.15).clip(0,1)
        for ii in np.arange(len(link_colors)): 
            link_colors[ii] = link_colors[ii].replace('opacity',str(opacities[ii]))

        # Link values
        link_values = df_transitions_mm2['Area'][sel].tolist()

        # Create a Sankey instance
        fig = go.Figure(data=[go.Sankey(
            node = dict(
                pad = 20,
                thickness = 30,
                line = dict(color = "black", width = 0),
                color = node_colors
            ),
            link = dict(
                line = dict(color='rgba(1,1,1,0)'),
                source = link_sources,
                target = link_targets,
                value = link_values,
                color = link_colors
            ))])

        # Write figure
        outfile = os.path.join(config['folder_stats'],'climatologies',scenario+'_sankey.pdf')
        fig.write_image(outfile)
        print(outfile)
        
    pdb.set_trace()
    
    
if __name__ == '__main__':
    main()