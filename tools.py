#!/usr/bin/env python  
# -*- coding: utf-8 -*-

__author__ = "Hylke E. Beck"
__email__ = "hylke.beck@gmail.com"
__date__ = "November 2022"

import os, sys, glob, time, pdb
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
from scipy import ndimage as nd
from scipy import stats
from matplotlib.patches import PathPatch
import matplotlib.pyplot as plt
import importlib 
from netCDF4 import Dataset
from skimage.transform import resize
import warnings
import rasterio
import scipy.io
import h5py
from rasterio.transform import from_origin
import zipfile
import fnmatch
import subprocess
import matplotlib.patches as patches
from matplotlib.ticker import EngFormatter, StrMethodFormatter
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection


def fill(data, invalid=None):
    """This function performs nearest neighbor interpolation gap fill on the input
    data. The input 'data' is expected to be a numpy array and 'invalid' is a 
    boolean numpy array with the same shape as 'data' that marks invalid/missing
    data points. The function uses the distance_transform_edt function from the
    ndimage module to calculate the nearest valid data point for each invalid
    data point, and returns the filled data. If the 'invalid' parameter is not
    provided, the function assumes that nans in 'data' are invalid.
    """
    
    if invalid is None: invalid = np.isnan(data)
    ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True)
    return data[tuple(ind)]  


def load_config(filepath):
    """This function loads a configuration file into a python dictionary. The
    filepath of the configuration file needs to be provided as an input.
    """
    
    config = importlib.import_module(filepath, package=None)
    config = config.config
    return config

    
def adjust_box_widths(g, fac):
    """This function adjusts the widths of a boxplot generated using seaborn.
    """

    # iterating through Axes instances
    for ax in g.axes:

        # iterating through axes artists:
        for c in ax.get_children():

            # searching for PathPatches
            if isinstance(c, PathPatch):
                # getting current width of box:
                p = c.get_path()
                verts = p.vertices
                verts_sub = verts[:-1]
                xmin = np.min(verts_sub[:, 0])
                xmax = np.max(verts_sub[:, 0])
                xmid = 0.5*(xmin+xmax)
                xhalf = 0.5*(xmax - xmin)

                # setting new width of box
                xmin_new = xmid-fac*xhalf
                xmax_new = xmid+fac*xhalf
                verts_sub[verts_sub[:, 0] == xmin, 0] = xmin_new
                verts_sub[verts_sub[:, 0] == xmax, 0] = xmax_new

                # setting new width of median line
                for l in ax.lines:
                    if np.all(l.get_xdata() == [xmin, xmax]):
                        l.set_xdata([xmin_new, xmax_new])
                        

def koppen_geiger(T,P,koppen_table):
    """This function classifies monthly temperature and precipitation
    climatologies according to the Koppen-Geiger climate classification. The
    inputs 'T' and 'P' represent the temperature and precipitation
    climatologies respectively, and should be provided as three-dimensional
    arrays with the third dimension representing time (12 months). The
    temperature data should be in units of degrees Celsius and the
    precipitation data in units of mm/month. The 'koppen_table' input is a
    table used to classify the climate. The function returns the climate class
    and major type of the area. 
    """

    # Make boolean indexing array to select summer months
    # Summer is defined here as the warmest 6-month period between Oct-Mar
    # (ONDJFM) and Apr-Sep (AMJJAS). In the Northern Hemisphere, summer will 
    # usually be from Apr-Sep, while in the Southern Hemisphere, summer will
    # usually be from Oct-Mar.
    T_ONDJFM = np.mean(T[np.array([0,1,2,9,10,11]),:,:],axis=0)
    T_AMJJAS = np.mean(T[np.array([3,4,5,6,7,8]),:,:],axis=0)
    tmp = T_AMJJAS>T_ONDJFM
    sum_index = np.zeros(T.shape,dtype=bool)
    sum_index[np.array([3,4,5,6,7,8]),:,:] = np.tile(tmp,(6,1,1))
    sum_index[np.array([0,1,2,9,10,11]),:,:] = np.tile(~tmp,(6,1,1))
    del tmp
    
    # Total P in summer and winter and P in driest month
    Pw = np.sum(P*(~sum_index).astype(np.single),axis=0)
    Ps = np.sum(P*sum_index.astype(np.single),axis=0)
    Pdry = np.min(P,axis=0)

    # P in wettest and driest months in summer and winter
    tmp = sum_index.astype(np.single)
    tmp[tmp==0] = np.NaN
    Psdry = np.nanmin(P*tmp,axis=0)
    Pswet = np.nanmax(P*tmp,axis=0)
    tmp = (~sum_index).astype(np.single) 
    tmp[tmp==0] = np.NaN
    Pwdry = np.nanmin(P*tmp,axis=0)
    Pwwet = np.nanmax(P*tmp,axis=0)

    # Mean annual temperature and preciptiation (MAT and MAP)
    # Number of months with temperature >10 degrees C
    # Temperature of hottest and coldest months
    MAT = np.mean(T,axis=0)
    MAP = np.sum(P,axis=0)
    Tmon10 = np.sum(T>10,axis=0)
    Thot = np.max(T,axis=0)
    Tcold = np.min(T,axis=0)

    # Threshold P
    Pthresh = 2*MAT+14
    Pthresh[Pw*2.333>Ps] = 2*MAT[Pw*2.333>Ps]
    Pthresh[Ps*2.333>Pw] = 2*MAT[Ps*2.333>Pw]+28

    # Classification of B classes
    B = MAP<10*Pthresh
    BW = (B) & (MAP<5*Pthresh)
    BWh = (BW) & (MAT>=18)
    BWk = (BW) & (MAT<18)
    BS = (B) & (MAP>=5*Pthresh)
    BSh = (BS) & (MAT>=18)
    BSk = (BS) & (MAT<18)

    # Classification of A classes
    A = (Tcold>=18) & (~B) # Added "& ~B"
    Af = (A) & (Pdry>=60)
    Am = (A) & (~Af) & (Pdry>=100-MAP/25)
    Aw = (A) & (~Af) & (Pdry<100-MAP/25)

    # Classification of C classes
    C = (Thot > 10) & (Tcold > 0) & (Tcold<18) & (~B)
    Cs = (C) & (Psdry<40) & (Psdry<Pwwet/3)
    Cw = (C) & (Pwdry<Pswet/10)
    overlap = (Cs) & (Cw)
    Cs[(overlap) & (Ps>Pw)] = 0
    Cw[(overlap) & (Ps<=Pw)] = 0
    Csa = (Cs) & (Thot>=22)
    Csb = (Cs) & (~Csa) & (Tmon10>=4)
    Csc = (Cs) & (~Csa) & (~Csb) & (Tmon10>=1) & (Tmon10<4)
    Cwa = (Cw) & (Thot>=22)
    Cwb = (Cw) & (~Cwa) & (Tmon10>=4)
    Cwc = (Cw) & (~Cwa) & (~Cwb) & (Tmon10>=1) & (Tmon10<4)
    Cf = (C) & (~Cs) & (~Cw)
    Cfa = (Cf) & (Thot>=22)
    Cfb = (Cf) & (~Cfa) & (Tmon10>=4)
    Cfc = (Cf) & (~Cfa) & (~Cfb) & (Tmon10>=1) & (Tmon10<4)

    # Classification of D classes
    D = (Thot>10) & (Tcold<=0) & (~B) # Added "& ~B"
    Ds = (D) & (Psdry<40) & (Psdry<Pwwet/3)
    Dw = (D) & (Pwdry<Pswet/10)
    overlap = (Ds) & (Dw)
    Ds[(overlap) & (Ps>Pw)] = 0
    Dw[(overlap) & (Ps<=Pw)] = 0
    Dsa = (Ds) & (Thot>=22)
    Dsb = (Ds) & (~Dsa) & (Tmon10>=4)
    Dsd = (Ds) & (~Dsa) & (~Dsb) & (Tcold<-38)
    Dsc = (Ds) & (~Dsa) & (~Dsb) & (~Dsd)

    Dwa = (Dw) & (Thot>=22)
    Dwb = (Dw) & (~Dwa) & (Tmon10>=4)
    Dwd = (Dw) & (~Dwa) & (~Dwb) & (Tcold<-38)
    Dwc = (Dw) & (~Dwa) & (~Dwb) & (~Dwd)
    Df = (D) & (~Ds) & (~Dw)
    Dfa = (Df) & (Thot>=22)
    Dfb = (Df) & (~Dfa) & (Tmon10>=4)
    Dfd = (Df) & (~Dfa) & (~Dfb) & (Tcold<-38)
    Dfc = (Df) & (~Dfa) & (~Dfb) & (~Dfd)

    # Classification of E classes
    E = (Thot <= 10) & (~B) # Added "& ~B", and replaced "Thot<10" with "Thot<=10"
    ET = (E) & (Thot>0)
    EF = (E) & (Thot<=0)
    
    # Make maps if final KG subclass and final KG major class
    # Using eval() is bad practice, but the code above is more concise
    # without dict
    Class = np.zeros((T.shape[1],T.shape[2]),dtype=np.single)*np.NaN
    Major = np.zeros((T.shape[1],T.shape[2]),dtype=np.single)*np.NaN
    for ii in np.arange(koppen_table.shape[0]):
        mask = eval(koppen_table['Symbol'][ii])
        Class[mask] = koppen_table['Class'][ii]
        Major[mask] = koppen_table['Major'][ii]

    # Make oceans NaN
    mask = np.isnan(P[0,:,:]+T[0,:,:])
    Class[mask] = np.NaN
    Major[mask] = np.NaN

    return {'Class': Class, 'Major': Major}
    
    
def write_to_netcdf_3d(file, data, varname, varunits, month, least_sig_dig):
    """This function writes a 3-dimensional data array to a netCDF file, or
    updates it if the file already exists. The input 'file' is the filepath of
    the netCDF file to be written to, 'data' is the 3-dimensional data array to
    be written, 'varname' is the variable name for the data, 'varunits' is the
    units of the data, 'month' is the month number of the data and
    'least_sig_dig' is the number of least significant digits that should be
    preserved in the data.
    """
    
    if os.path.isfile(file)==False:

        res = 360/float(data.shape[1])
        lon = np.arange(data.shape[1], dtype=np.single)*res - 180.0 + res/2
        lat = 90-np.arange(data.shape[0], dtype=np.single)*res - res/2
        
        if os.path.isdir(os.path.dirname(file))==False:
            os.makedirs(os.path.dirname(file))
        
        ncfile = Dataset(file, 'w', format='NETCDF4')
        ncfile.history = 'Created on %s' % datetime.utcnow().strftime('%Y-%m-%d %H:%M')

        ncfile.createDimension('lon', len(lon))
        ncfile.createDimension('lat', len(lat))
        ncfile.createDimension('time', 12)

        ncfile.createVariable('lon', 'f4', ('lon',))
        ncfile.variables['lon'][:] = lon
        ncfile.variables['lon'].units = 'degrees_east'
        ncfile.variables['lon'].long_name = 'longitude'

        ncfile.createVariable('lat', 'f4', ('lat',))
        ncfile.variables['lat'][:] = lat
        ncfile.variables['lat'].units = 'degrees_north'
        ncfile.variables['lat'].long_name = 'latitude'

        ncfile.createVariable('time', 'f4', 'time')
        ncfile.variables['time'][:] = np.arange(1,13)
        ncfile.variables['time'].units = 'month'
        ncfile.variables['time'].long_name = 'month of year'
    
    else:
        ncfile = Dataset(file, 'r+', format='NETCDF4')   

    if varname not in ncfile.variables.keys():
        ncfile.createVariable(varname, data.dtype, ('time', 'lat', 'lon'),\
            zlib=True, complevel=1,\
            chunksizes=(1,int(np.minimum(data.shape[0],200)),int(np.minimum(data.shape[1],200)),),\
            fill_value=-9999, least_significant_digit=least_sig_dig)

    ncfile.variables[varname][month-1,:,:] = data
    ncfile.variables[varname].units = varunits

    ncfile.close()
    

def write_to_netcdf_2d(file, data, varname, varunits, least_sig_dig):
    """This function writes a 2-dimensional data array to a netCDF file, or
    updates it if the file already exists. The input 'file' is the filepath of
    the netCDF file to be written to, 'data' is the 2-dimensional data array to
    be written, 'varname' is the variable name for the data, 'varunits' is the
    units of the data, 'least_sig_dig' is the number of least significant
    digits that should be preserved in the data.
    """
    
    if os.path.isfile(file)==False:

        res = 360/float(data.shape[1])
        lon = np.arange(data.shape[1], dtype=np.single)*res - 180.0 + res/2
        lat = 90-np.arange(data.shape[0], dtype=np.single)*res - res/2
        
        if os.path.isdir(os.path.dirname(file))==False:
            os.makedirs(os.path.dirname(file))
        
        ncfile = Dataset(file, 'w', format='NETCDF4')
        ncfile.history = 'Created on %s' % datetime.utcnow().strftime('%Y-%m-%d %H:%M')

        ncfile.createDimension('lon', len(lon))
        ncfile.createDimension('lat', len(lat))

        ncfile.createVariable('lon', 'f4', ('lon',))
        ncfile.variables['lon'][:] = lon
        ncfile.variables['lon'].units = 'degrees_east'
        ncfile.variables['lon'].long_name = 'longitude'

        ncfile.createVariable('lat', 'f4', ('lat',))
        ncfile.variables['lat'][:] = lat
        ncfile.variables['lat'].units = 'degrees_north'
        ncfile.variables['lat'].long_name = 'latitude'

    else:
        ncfile = Dataset(file, 'r+', format='NETCDF4')   

    if varname not in ncfile.variables.keys():
        ncfile.createVariable(varname, data.dtype, ('lat', 'lon'),\
            zlib=True, complevel=1,\
            chunksizes=(int(np.minimum(data.shape[0],200)),int(np.minimum(data.shape[1],200)),),\
            fill_value=-9999, least_significant_digit=least_sig_dig)

    ncfile.variables[varname][:,:] = data
    ncfile.variables[varname].units = varunits

    ncfile.close()
    
    
def produce_change_map(reference_map,reference_period,monthly_data, \
    monthly_dates,target_period,target_month,varname,change_offset, \
    change_limits):
    """This function produces two maps, one showing the change between a 
    reference period and a target period, and the other showing the target 
    map itself. The input 'reference_map' is a reference map used to create 
    the target map, 'reference_period' is a tuple specifying the start and 
    end years of the reference period, 'monthly_data' is a 3-dimensional 
    array of the data, 'monthly_dates' is an array of dates corresponding to 
    the monthly data, 'target_period' is a tuple specifying the start and 
    end years of the target period, 'target_month' is the month for which 
    the change map is to be created, 'varname' is the variable name for the 
    data, 'change_offset' is an offset to be added to the data for the 
    calculation, and 'change_limits' is a tuple specifying the minimum and 
    maximum values for the change map. 
    """
    
    ind1 = (monthly_dates.year>=target_period[0])\
        & (monthly_dates.year<=target_period[1])\
        & (monthly_dates.month==target_month)
    ind2 = (monthly_dates.year>=reference_period[0])\
        & (monthly_dates.year<=reference_period[1])\
        & (monthly_dates.month==target_month)
    
    if varname=='P':
        change_map = (np.mean(monthly_data[:,:,ind1],axis=2)+change_offset) \
            /(np.mean(monthly_data[:,:,ind2],axis=2)+change_offset)
        change_map[np.isnan(change_map)] = 1
    elif varname=='Temp':
        change_map = np.mean(monthly_data[:,:,ind1],axis=2) \
            -np.mean(monthly_data[:,:,ind2],axis=2)
        change_map[np.isnan(change_map)] = 0
    
    change_map = change_map.clip(change_limits[0],change_limits[1]).astype(np.single)
    change_map = resize(change_map,reference_map.shape,order=1,mode='constant', \
        anti_aliasing=False)
    
    if varname=='P':
        target_map = reference_map*change_map
    elif varname=='Temp':
        target_map = reference_map+change_map
    
    return {'target_map':target_map,'change_map':change_map}
    
    
def compute_ens_mean_std(ens_dir,vars,mapsize,subsetsize,skip_existing):
    """This function computes the ensemble mean and standard deviation of a 
    set of netCDF files in a given directory. The input 'ens_dir' is the 
    directory containing the netCDF files, 'vars' is a list of variable 
    names, 'mapsize' is the size of the map, 'subsetsize' is the size of the 
    subsets to be used in the computation, and 'skip_existing' is a boolean 
    value indicating whether to skip the computation if the ensemble mean 
    and standard deviation files already exist. The function first creates 
    two output files, one for the ensemble mean and one for the ensemble 
    standard deviation. It then loops through the variables, months and 
    subsets, reading the data from the input files, computing the mean and 
    standard deviation of each subset, and writing the results to the output 
    files.
    """
    
    suffix = str(180/mapsize[0]).replace('.','p')
    
    if (os.path.isfile(os.path.join(ens_dir,'ensemble_mean_'+suffix+'.nc'))) \
            & (os.path.isfile(os.path.join(ens_dir,'ensemble_std_'+suffix+'.nc'))) \
            & (skip_existing==True):
        return
        
    if os.path.isfile(os.path.join(ens_dir,'ensemble_mean_'+suffix+'.nc')):
        os.remove(os.path.join(ens_dir,'ensemble_mean_'+suffix+'.nc'))
    if os.path.isfile(os.path.join(ens_dir,'ensemble_std_'+suffix+'.nc')):
        os.remove(os.path.join(ens_dir,'ensemble_std_'+suffix+'.nc'))
    
    files = glob.glob(os.path.join(ens_dir,'*.nc'))
    ncoutmean = os.path.join(ens_dir,'ensemble_mean_'+suffix+'.nc')
    ncoutstd = os.path.join(ens_dir,'ensemble_std_'+suffix+'.nc')
    for vv in np.arange(len(vars)):
        for month in np.arange(1,13):
            data_mean = np.zeros(mapsize,dtype=np.single)*np.NaN
            data_std = np.zeros(mapsize,dtype=np.single)*np.NaN
            for substart in np.arange(0,mapsize[1],subsetsize):
                subend = substart+subsetsize
                data_mems = np.zeros((mapsize[0],subsetsize,len(files)),dtype=np.single)*np.NaN
                for ii in np.arange(len(files)):
                    print('Loading '+vars[vv][0]+ ' month '+str(month)+ ' substart '+str(substart)+' '+files[ii])
                    dset = Dataset(files[ii])
                    data_mems[:,:,ii] = np.array(dset.variables[vars[vv][1]][month-1,:,substart:subend],dtype=np.single)
                    dset.close()
                data_mean[:,substart:subend] = np.mean(data_mems,axis=2)
                data_std[:,substart:subend] = np.std(data_mems,axis=2)
            write_to_netcdf_3d(ncoutmean,data_mean,vars[vv][1],vars[vv][2],month,1)
            write_to_netcdf_3d(ncoutstd,data_std,vars[vv][1],vars[vv][2],month,1)


def compute_kg_maps(ens_dir,koppen_table,mapsize,subsetsize,skip_existing):
    """Computes the Koppen-Geiger map of a given ensemble 
    directory using a Koppen-Geiger table. It takes the ensemble directory, 
    Koppen-Geiger table, map size, subset size, and a flag for skipping 
    existing files as inputs. It loops through each subset of the map, 
    computes the Koppen-Geiger class and confidence of each ensemble member 
    using the koppen_geiger function, and computes the mode and confidence 
    of the ensemble. It then writes the Koppen-Geiger class and confidence 
    maps to a netCDF file.
    """

    suffix = str(180/mapsize[0]).replace('.','p')
    
    if (os.path.isfile(os.path.join(ens_dir,'koppen_geiger_'+suffix+'.nc'))) \
            & (skip_existing==True):
        return
        
    if os.path.isfile(os.path.join(ens_dir,'koppen_geiger_'+suffix+'.nc')):
        os.remove(os.path.join(ens_dir,'koppen_geiger_'+suffix+'.nc'))
    
    # List of netCDF files with P and T climatologies
    # Exclude ensemble_mean_*.nc and ensemble_std_*.nc from list
    files = glob.glob(os.path.join(ens_dir,'*.nc'))
    files = [x for x in files if not 'ensemble_mean' in x and not 'ensemble_std' in x and not 'koppen_geiger' in x]
    
    ncout = os.path.join(ens_dir,'koppen_geiger_'+suffix+'.nc')
    kg_class = np.zeros(mapsize,dtype=np.int8)
    kg_confidence = np.zeros(mapsize,dtype=np.int8)
    for substart in np.arange(0,mapsize[1],subsetsize):
        subend = substart+subsetsize
        kg_class_ens = np.zeros((mapsize[0],subsetsize,len(files)),dtype=np.int8)
        for ii in np.arange(len(files)):
            dset = Dataset(files[ii])
            data_T = np.array(dset.variables['air_temperature'][:,:,substart:subend],dtype=np.single)
            data_P = np.array(dset.variables['precipitation'][:,:,substart:subend],dtype=np.single)
            dset.close()
            kg_class_ens[:,:,ii] = koppen_geiger(data_T,data_P,koppen_table)['Class']
        kg_class[:,substart:subend] = mode(kg_class_ens,axis=2)
        kg_confidence[:,substart:subend] = confidence(kg_class_ens,axis=2).astype(np.int8)
    write_to_netcdf_2d(ncout,kg_class,'kg_class','',1)
    write_to_netcdf_2d(ncout,kg_confidence,'kg_confidence','%',1)
    

def periods_to_datetimes(periods):
    """Converts an array of Pandas periods datatype to an 
    array of datetime datatype.
    """
    datetimes = np.empty(periods.shape, dtype='datetime64[s]')
    for ii in np.arange(len(periods)):
        try:
            datetimes[ii] = pd.to_datetime(pd.Period.to_timestamp(periods[ii]))
        except:
            datetimes[ii] = np.datetime64('NaT')
    return datetimes
    
    
def mode(ndarray, axis=0):
    """Computes the mode along an axis of a multi-dimensional 
    array. It takes a multi-dimensional array and an axis as inputs, and 
    returns an array of the mode along the given axis. Replaces 
    scipy.stats.mode which is extremely slow. This solution is still five 
    times slower than Matlab.
    """
    
    return np.apply_along_axis(lambda x: np.bincount(x).argmax(), axis=axis, arr=ndarray)


def confidence(ndarray, axis=0):
    """This function computes the fraction of the array equal to the mode. 
    It takes a multi-dimensional array and an axis as inputs, and returns an 
    array of the confidence along the given axis. 
    """
    
    return np.apply_along_axis(lambda x: 100*np.bincount(x).max()/len(x), axis=axis, arr=ndarray)


def modefun(x,nanint):
    """Computes the mode for 1-dimensional arrays. nanint sets the NaN value 
    for integer arrays. 
    """
    
    x = x[x!=nanint]
    if len(x)==0:
        return nanint
    else:
        return np.bincount(x).argmax()
        
 
def mapresize(A,newshape,measure='mean',nantol=0.75,nanint=0):
    """ This function computes the mode of a one-dimensional array and 
    replaces NaN values with a given integer. It takes a one-dimensional 
    array and an integer as inputs, and returns the mode of the array with 
    NaN values replaced by the given integer.
    """
    
    dtype = A.dtype
    oldshape = A.shape
    
    assert int(oldshape[0]/newshape[0])==int(oldshape[1]/newshape[1]), \
        "Input and output shapes are not proportional"
    
    # Reshape input array such that third dimension contains data of each new grid-cell
    factor = int(oldshape[0]/newshape[0])
    B = np.reshape(A,(newshape[0],oldshape[1],factor),order='C')
    C = np.reshape(B,(newshape[0],newshape[1],factor*factor),order='F')

    # If input data type is integer, convert to single
    if (np.issubdtype(A.dtype, np.integer)) & (measure!='mode'):        
        C = C.astype(np.single)
        C[C==nanint] = np.NaN
        
    # Compute mean, median, mode, max, or min in third dimension
    if measure=='mean': 
        D = np.nanmean(C, axis=2)
    elif measure=='median':
        D = np.nanmedian(C, axis=2)
    elif measure=='mode':
        if np.issubdtype(A.dtype, np.integer)==False:
            raise Exception("Only integer input using measure='mode'")
        if np.sum(A<0)>0:
            raise Exception("No numbers below zero using measure='mode'")
        D = np.apply_along_axis(modefun, axis=2, arr=C, nanint=nanint)
        D = D.astype(A.dtype)
    elif measure=='max':
        D = np.nanmax(C, axis=2)
    elif measure=='min':
        D = np.nanmin(C, axis=2)        
    
    # Compute mask based on nantol
    mask = np.sum(np.isnan(C),axis=2)>(C.shape[2]*nantol)    
    if np.issubdtype(A.dtype, np.integer):
        D[mask] = nanint
        D = D.astype(dtype)
    else:
        D[mask] = np.NaN
        
    return D
    

def write_to_geotiff(file,data,cmap,nodata):
    """exports a numpy array to a geoTIFF file. It writes the data to the 
    file and also writes a color map to the file. If the file already 
    exists, it will be removed before the new data is written. The data type
    is np.uint8 to avoid issues with the colormap.
    """
    
    if os.path.isfile(file):
        os.remove(file)
        
    with rasterio.open(
            file, 'w', driver='GTiff', dtype=np.uint8,
            width=data.shape[1], height=data.shape[0], count=1, crs='EPSG:4326',
            transform=from_origin(-180.0, 90.0, 360/data.shape[1], 180/data.shape[0]),
            nodata=nodata, tiled=True, compress='lzw') as dataset:
        dataset.write(data, indexes=1)
        dataset.write_colormap(1, cmap)
        dataset.close()
        
        
def zip_folder(zip_file,folder,pattern,compress_type):
    """Creates a zip file of files in a folder that match a given pattern. 
    It uses the python zipfile module to create the zip file and writes 
    files in the folder to the zip file using the os.walk function to 
    recursively search through subdirectories.
    """

    with zipfile.ZipFile(zip_file, 'w') as zipObj:
        for root, dirs, files in os.walk(folder):
            for file in files:
                if fnmatch.fnmatch(file,pattern):
                    print(os.path.join(root,file))
                    zipObj.write(
                        os.path.join(root,file), 
                        os.path.relpath(os.path.join(root,file), folder),
                        compress_type=compress_type
                        )
                                
                                
def sync_data(sync_cmd,dir_local,dir_remote):
    """copies data from a local directory to a remote directory using a 
    shell command. The command and the local and remote directories are 
    passed as arguments to the function.
    """

    # Insert folders into sync_cmd
    sync_cmd = sync_cmd.replace('$dir_local',dir_local).replace('$dir_remote',dir_remote)
    
    # Execute sync command
    print("Executing: '"+sync_cmd+"'")    
    out = subprocess.Popen(sync_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    stdout,stderr = out.communicate()
    print(stdout.decode("utf-8"))
    
    
def plot_map(data,data_extent,figout,figdims,cmap,plot_extent,lims,interpolation,shp,color,show_axes):
    """Creates a map figure of a data array. Inputs include the data, the 
    data extent, the figure dimensions, the color map, the plot extent and 
    limits, the interpolation type, shapefiles, color and whether to show 
    the axes or not.
    """
    
    # Plot gridded map
    data = np.ma.array(data, mask=np.isnan(data))
    ax1 = plt.imshow(data, interpolation='none', extent=data_extent, alpha=1, cmap=cmap, vmin=lims[0], vmax=lims[1])
    plt.ylim(plot_extent[2], plot_extent[3])
    plt.xlim(plot_extent[0], plot_extent[1])
    if show_axes==False:
        plt.axis('off')
    plt.gca().yaxis.set_major_formatter(StrMethodFormatter(u"{x:.0f}°"))
    plt.gca().xaxis.set_major_formatter(StrMethodFormatter(u"{x:.0f}°")) 
    
    plt.gca().set_aspect('equal')
    
    # Overlay country shapefile
    Nshp = len(shp)
    shapes  = shp.shapes()
    for nshp in np.arange(Nshp):
        ptchs = []
        pts = np.array(shapes[nshp].points)
        prt = shapes[nshp].parts
        par = list(prt)+[pts.shape[0]]
        for pij in np.arange(len(prt)):
            pts_sub = pts[par[pij]:par[pij+1]]
            plt.plot(pts_sub[:,0],pts_sub[:,1],color,linewidth=0.1,alpha=1)
    
    # Write to png
    plt.gcf().set_size_inches(figdims)
    plt.savefig(figout,dpi=600,bbox_inches='tight',pad_inches=0.0)
    plt.close()
    

def readmatfile(filepath, var): 
    """This function loads a Matlab file of any version and returns the 
    specified variable in the file. 
    """
    try: 
        f = h5py.File(filepath, 'r') 
        data = f.get(var)[()]
        data = data.transpose() 
        f.close() 
    except:
        try: 
            mat = scipy.io.loadmat(filepath)
            data = eval("mat['"+var.replace("/","'][0,0]['")+"'][:]")
        except: 
            pass 

    return data


def eliminate_trailing_zeros(data):
    """Check for long series of erroneous zeros in GSOD station records.
    """
    
    mov_avg_365 = pd.DataFrame(data).rolling(min_periods=182,window=365,center=True).mean().values
    mov_min_365 = pd.DataFrame(data).rolling(min_periods=182,window=365,center=True).min().values
    data[np.isnan(mov_avg_365+mov_min_365)] = np.NaN
    data[mov_avg_365==0] = np.NaN # Discard values in periods with only zeros
    data[mov_min_365!=0] = np.NaN # Discard values in periods with no zeros    
    return data
    
    
def mean_valid_obs(x):
    """This function takes in a data series, x, and returns the mean of the 
    series if there are at least 67% of valid observations, otherwise 
    returns NaN.
    """

    min_obs = 0.67 * x.index[0].days_in_month
    valid_obs = x.notnull().sum()
    if valid_obs < min_obs:
        return np.nan
    return x.mean()
    

def compute_monthly_climatology(data_daily,dates_daily):
    """Given a daily time series data and corresponding dates, this function 
    computes the monthly climatology. It first computes a monthly time 
    series by averaging the daily data for each month. Then it computes the 
    monthly climatology by averaging the monthly time series values over 
    multiple years for each month. The function returns an array of 12 
    monthly climatology values.
    """
    
    # Preparation
    dates_daily_year = dates_daily.year
    dates_daily_month = dates_daily.month
    dates_monthly = pd.date_range(start=dates_daily[0], end=dates_daily[-1], freq='MS')
    dates_monthly_months = dates_monthly.month
    
    # Compute monthly time series
    data_monthly = np.zeros(dates_monthly.shape)*np.NaN
    data_daily_nan = np.isnan(data_daily).flatten()
    for ii in np.arange(len(dates_monthly)):
        sel = (dates_daily_year==dates_monthly[ii].year) & (dates_daily_month==dates_monthly[ii].month) & (data_daily_nan==False)
        nobs = np.sum(sel)
        if nobs>=21:
            data_monthly[ii] = np.mean(data_daily[sel])
    
    # Compute monthly climatology
    monthly_clim = np.zeros((12,))*np.NaN
    for month in np.arange(1,13):
        sel = dates_monthly_months==month
        nobs = np.sum(np.isnan(data_monthly[sel])==False)
        if nobs>=10:
            monthly_clim[month-1] = np.nanmean(data_monthly[sel])
     
    return monthly_clim