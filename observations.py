# 3rd party modules:
import numpy as np
import netCDF4 as nc
import random
import datetime as dt
# local modules:
import experiment_setup as es


def find_nearest_idx_tol(array, value, tol=dt.timedelta(days=1.)):
    """
    Find nearest value in an array for a given tolerance
    :param array: array of values (arr)
    :param value: value for which to find nearest element (float, obj)
    :param tol: distance tolerance (float, obj)
    :return: nearest value in array, index of nearest value
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    if abs(array[idx] - value) <= tol:
        ret_val = idx
    else:
        ret_val = np.nan
    return ret_val


def extract_twin_data(obs='/group_workspaces/jasmin2/hydro_jules/data/epinnington/smap_9km_uk/smap_9km_2016.nc'):
    """
    Function for extracting observations to be assimilated in data assimilation experiments. Here we are running a twin
    experiment using a "model truth" to draw observations from.
    :param mod_truth: location of netCDF file containing "model truth" output (str)
    :param seed_val: seed value for adding random noise to the observations (int)
    :return: dictionary containing observations and observation errors (dictionary)
    """
    # open obs netCDF file
    smap = nc.Dataset(obs, 'r')
    # set position of observations
    smap_time_var = smap.variables['Soil_Moisture_Retrieval_Data_AM_tb_time_seconds']
    smap_sm = smap.variables['Soil_Moisture_Retrieval_Data_AM_soil_moisture'][:, 64, 69]
    smap_times = nc.num2date(smap_time_var[:,64,69], 'seconds since 2000-01-01 20:00:00')
    smap_none_idx = np.where(smap_times != None)[0]
    smap_times = smap_times[smap_none_idx]
    smap_sm = smap_sm[smap_none_idx]
    sm_obs = smap_sm[smap_sm.nonzero()]
    smap_times = smap_times[smap_sm.nonzero()]

    sm_err = np.ones(len(sm_obs)) * np.mean(sm_obs[sm_obs > 0]) * 0.1  # 10% error in mod obs
    # close netCDF file
    smap.close()
    return {'obs': sm_obs, 'obs_err': sm_err, 'sm_obs': sm_obs, 'sm_err': sm_err, 'sm_times': smap_times}


def extract_jules_hx(nc_file):
    """
    Function extracting the modelled observations from JULES netCDF files
    :param nc_file: netCDF file to extract observations from (str)
    :return: dictionary containing observations (dict)
    """
    nc_dat = nc.Dataset(nc_file, 'r')
    obs = '/group_workspaces/jasmin2/hydro_jules/data/epinnington/smap_9km_uk/smap_9km_2016.nc'
    smap = nc.Dataset(obs, 'r')
    smap_time_var = smap.variables['Soil_Moisture_Retrieval_Data_AM_tb_time_seconds']
    smap_sm = smap.variables['Soil_Moisture_Retrieval_Data_AM_soil_moisture'][:, 64, 69]
    smap_times = nc.num2date(smap_time_var[:,64,69], 'seconds since 2000-01-01 20:00:00')
    smap_none_idx = np.where(smap_times != None)[0]
    smap_times = smap_times[smap_none_idx]
    smap_sm = smap_sm[smap_none_idx]
    sm_obs = smap_sm[smap_sm.nonzero()]
    smap_times = smap_times[smap_sm.nonzero()]
    # set position of observations to correspond with assimilated observations
    jules_times = nc.num2date(nc_dat.variables['time'][:], nc_dat.variables['time'].units)
    sm_ob_position = np.array([find_nearest_idx_tol(jules_times, smap_t) for smap_t in smap_times])
    # extract observations from netCDF file
    sm_hx = nc_dat.variables['smcl'][sm_ob_position, 0, 0, 0]/50.
    # close netCDF file
    nc_dat.close()
    smap.close()
    return {'obs': sm_hx}


def extract_jules_hxb():
    """
    Function to extract modelled observations for prior JULES run
    :return: dictionary containing observations (dict)
    """
    return extract_jules_hx(es.output_directory+'/background/xb0.daily.nc')


def extract_jules_hxb_ens():
    """
    Function to extract ensemble of modelled observations from prior model ensemble
    :return: ensemble of modelled observations (lst)
    """
    hm_xbs = []
    for xb_fname in xrange(0, es.ensemble_size):
        hxbi_dic = extract_jules_hx(es.output_directory + '/ensemble0/ens' + str(xb_fname) + '.daily.nc')
        hxbi = hxbi_dic['obs']
        hm_xbs.append(hxbi)
    return hm_xbs
