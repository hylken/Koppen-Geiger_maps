config = {
    'mapsize': (18000,36000),
    'upscale_mapsizes': [(1800,3600),(360,720),(180,360)],
    'periods_historical': [(1901,1930),(1931,1960),(1961,1990),(1991,2020)],
    'periods_future': [(2041,2070),(2071,2099)],
    'version': 'V2',
    'skip_existing': True,
    'Pdatasets': {'WorldClim_V2': (1970,2000), 'CHPclim_V1': (1980,2009), 'CHELSA_V12': (1979,2013), 'CHELSA_V21': (1981,2010)},
    'Tdatasets': {'WorldClim_V2': (1970,2000), 'CHELSA_V12': (1979,2013), 'CHELSA_V21': (1981,2010)},
    'folder_out': '/mnt/datawaha/hyex/beckhe/RESEARCH/Paper_30_New_KG_maps',
    'folder_stats': '/home/beckhe/Koppen-Geiger_maps/stats_figs_tables', 
    'folder_maps': '/mnt/datawaha/hyex/beckhe/RESEARCH/Data/MAPS', 
    'folder_dataraw': '/mnt/datawaha/hyex/beckhe/DATA_RAW',
    'folder_dataproc': '/mnt/datawaha/hyex/beckhe/DATA_PROCESSED',
    'folder_station': '/mnt/datawaha/hyex/beckhe/DATA_PROCESSED/station_data',
    'vars': [
        ['P','precipitation','mm month-1','pr',(0.67,1.5)],
        ['Temp','air_temperature','°C','tas',(-15,15)]
    ],
    'perform_sync': True,
    'sync_cmd': "rclone copy -v --include='*.zip' $dir_local GoogleDrive:temp/Koppen-Geiger_maps",
    }
