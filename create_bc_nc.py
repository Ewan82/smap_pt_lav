import numpy as np
import netCDF4 as nc
import pickle
import itertools as itt
import shutil as sh
#third party modules
import cosby_pedo_tran as cpt



def create_bc_nc():
    soil_rec = pickle.load(open('hwsd_mu_texture.p', 'rb'))
    #soil_rec = np.recfromcsv('emma_hwsd_bc/HWSD_DATA.txt', usecols=(1,6,24,25,26,41,42,43))
    hwsd = nc.Dataset('emma_hwsd_bc/hwsd_uk_1km_grid.nc', 'r')
    bc_params = nc.Dataset('bc_params.nc', 'a')
    j=0
    for mu in np.unique(hwsd.variables['mu'][:])[:-1]:
        j+=1
        print j, ' out of ', len(np.unique(hwsd.variables['mu'][:])[:-1])
    #for mu in [10770]:
        t_sand = soil_rec['t_sand'][soil_rec['mu_global'] == mu][0]
        t_silt = soil_rec['t_silt'][soil_rec['mu_global'] == mu][0]
        t_clay = soil_rec['t_clay'][soil_rec['mu_global'] == mu][0]
        s_sand = soil_rec['s_sand'][soil_rec['mu_global'] == mu][0]
        s_silt = soil_rec['s_silt'][soil_rec['mu_global'] == mu][0]
        s_clay = soil_rec['s_clay'][soil_rec['mu_global'] == mu][0]
        t_bc = cpt.cosby(t_clay/100., t_sand/100., t_silt/100.)
        s_bc = cpt.cosby(s_clay/100., s_sand/100., s_silt/100.)
        loc = np.where(hwsd.variables['mu'][:] == float(mu))
        loc_itt = itt.izip(loc[1], loc[2])
        keys = ['b', 'vsat', 'sathh', 'satcon', 'vwilt', 'vcrit', 'hcap', 'hcon']
        for k in enumerate(keys):
            for i in loc_itt:
                #print k[1]
                #print bc_params.variables[k[1]][:2, i[0], i[1]]
                #print t_bc[k[0]]
                #print bc_params.variables[k[1]][2:, i[0], i[1]]
                #print s_bc[k[0]]
                bc_params.variables[k[1]][:2, i[0], i[1]] = t_bc[k[0]]
                bc_params.variables[k[1]][2:, i[0], i[1]] = s_bc[k[0]]
    bc_params.close()
    hwsd.close()
    return 'soil params updated!'


def create_bc_nc_params(bc_fname, param_lst=[15.70, 0.3, 0.037, 0.142, 0.63, 1.58, 0.64, 1.26]):
    soil_rec = pickle.load(open('hwsd_mu_texture.p', 'rb'))
    #soil_rec = np.recfromcsv('emma_hwsd_bc/HWSD_DATA.txt', usecols=(1,6,24,25,26,41,42,43))
    hwsd = nc.Dataset('emma_hwsd_bc/hwsd_uk_1km_grid.nc', 'r')
    sh.copy('bc_params.nc', bc_fname)
    bc_params = nc.Dataset(bc_fname, 'a')
    j=0
    for mu in np.unique(hwsd.variables['mu'][:])[:-1]:
        j+=1
        print j, ' out of ', len(np.unique(hwsd.variables['mu'][:])[:-1])
    #for mu in [10770]:
        t_sand = soil_rec['t_sand'][soil_rec['mu_global'] == mu][0]
        t_silt = soil_rec['t_silt'][soil_rec['mu_global'] == mu][0]
        t_clay = soil_rec['t_clay'][soil_rec['mu_global'] == mu][0]
        s_sand = soil_rec['s_sand'][soil_rec['mu_global'] == mu][0]
        s_silt = soil_rec['s_silt'][soil_rec['mu_global'] == mu][0]
        s_clay = soil_rec['s_clay'][soil_rec['mu_global'] == mu][0]
        t_bc = cpt.cosby(t_clay/100., t_sand/100., t_silt/100., b=param_lst[0], c=param_lst[1], e=param_lst[2],
                         f=param_lst[3], h=param_lst[4], i=param_lst[5], k=param_lst[6], l=param_lst[7])
        s_bc = cpt.cosby(s_clay/100., s_sand/100., s_silt/100., b=param_lst[0], c=param_lst[1], e=param_lst[2],
                         f=param_lst[3], h=param_lst[4], i=param_lst[5], k=param_lst[6], l=param_lst[7])
        loc = np.where(hwsd.variables['mu'][:] == float(mu))
        loc_itt = itt.izip(loc[1], loc[2])
        keys = ['b', 'vsat', 'sathh', 'satcon', 'vwilt', 'vcrit', 'hcap', 'hcon']
        for k in enumerate(keys):
            for i in loc_itt:
                #print k[1]
                #print bc_params.variables[k[1]][:2, i[0], i[1]]
                #print t_bc[k[0]]
                #print bc_params.variables[k[1]][2:, i[0], i[1]]
                #print s_bc[k[0]]
                bc_params.variables[k[1]][:2, i[0], i[1]] = t_bc[k[0]]
                bc_params.variables[k[1]][2:, i[0], i[1]] = s_bc[k[0]]
    bc_params.close()
    hwsd.close()
    return 'soil params updated!'