# 3rd party python modules:
import numpy as np
import datetime as dt
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns
import glob
import pickle
import scipy.stats as spst
# local modules:
import fourdenvar
import experiment_setup as es


def extract_ens_obs(nc_dir, var, lvl_idx=0):
    """
    Function that extracts modelled observations from an ensemble of JULES outputs
    :param nc_dir: directory where JULES ensemble output is saved (str)
    :param var: output variable to extract for each ensemble member (str)
    :param lvl_idx: level index for variable, if no level index exists this is ignored (int)
    :return: JULES output variable for each ensemble member (arr)
    """
    xbs = []
    for xb_fname in glob.glob(nc_dir+'/*.nc'):
        xbi_nc = nc.Dataset(xb_fname, 'r')
        if len(xbi_nc.variables[var].shape) == 3:
            xi = xbi_nc.variables[var][:, 0, 0]
        elif var == 'gpp':
            xi = 1000 * 60 * 60 * 24 * xbi_nc.variables[var][:, lvl_idx, 0, 0]
        elif var == 'smcl':
            xi = xbi_nc.variables[var][:, lvl_idx, 0, 0]/50.
        elif var == 'cropyield':
            xi = np.max(xbi_nc.variables[var][:, lvl_idx, 0, 0])
        elif var == 'cropstemc':
           xi = xbi_nc.variables['cropstemc'][:, lvl_idx, 0, 0] + xbi_nc.variables['cropreservec'][:, lvl_idx, 0, 0]
        else:
            xi = xbi_nc.variables[var][:, lvl_idx, 0, 0]
        xbs.append(xi)
        xbi_nc.close()
    return np.array(xbs)


def calc_mean_upp_low(x_ens):
    """
    Given an ensemble of model output will return the mean, mean + 1 stdev and mean - 1 stdev
    :param x_ens: ensemble of model output (arr)
    :return: mean, mean+1stdev, mean-1stdev
    """
    x_mean = np.mean(x_ens, axis=0)
    x_std = np.std(x_ens, axis=0)
    up_x = x_mean + x_std
    low_x = x_mean - x_std
    up_x[up_x < 0.0] = 0.0
    low_x[low_x < 0.0] = 0.0
    return x_mean, up_x, low_x


def find_nearest_idx_tol(array, value, tol=dt.timedelta(days=1.)):
    """
    Find nearest value in an array
    :param array: array of values
    :param value: value for which to find nearest element
    :return: nearest value in array, index of nearest value
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    if abs(array[idx] - value) <= tol:
        ret_val = idx
    else:
        ret_val = np.nan
    return ret_val


def plot_twin_spread(times, var, xb_dir, xa_dir, xt_var=None, obs=None, ob_times=None, err=None, ylab=None,
                     lvl_idx=0, axes=None):
    """
    Function to plot data assimilation prior and posterior experiment output for given variable
    :param times: array of times for x-axis (arr)
    :param var: JULES model variable to plot (str)
    :param xb_dir: directory of prior ensemble output (str)
    :param xa_dir: directory of posterior ensemble output (str)
    :param xt_var: model truth model output for plotted variable !optional! (arr)
    :param obs: assimilated observations !optional! (arr)
    :param ob_pos: position of assimilated observations !optional! (arr)
    :param err: assimilated observation error !optional! (arr)
    :param ylab: label for y-axis (str)
    :param lvl_idx: level index for JULES model output (int)
    :param axes: axes to plot on !optional! (obj)
    :return: figure and axis objects or just axis object, dependent on if axes arg is specified
    """
    sns.set_context('poster')
    sns.set_style('whitegrid')
    palette = sns.color_palette("colorblind", 11)
    if axes is None:
        fig, ax = plt.subplots(nrows=1, ncols=1,)
        ret_val = fig, ax
    elif axes is not None:
        ax = axes
        ret_val = ax
    xb_ens = extract_ens_obs(xb_dir, var, lvl_idx=lvl_idx)
    xa_ens = extract_ens_obs(xa_dir, var, lvl_idx=lvl_idx)
    xb_mean, up_xb, low_xb = calc_mean_upp_low(xb_ens)
    xa_mean, up_xa, low_xa = calc_mean_upp_low(xa_ens)
    ax.plot(times, xb_mean, '-', color=palette[0], label='prior')
    ax.fill_between(times, low_xb, up_xb, facecolor=palette[0],
                       alpha=0.3, linewidth=0.0)
    ax.plot(times, xa_mean, '-', color=palette[2], label='posterior')
    ax.fill_between(times, low_xa, up_xa, facecolor=palette[2],
                       alpha=0.3, linewidth=0.0)

    if xt_var is not None:
        ax.plot(times, xt_var, ':', label='model truth', color='k', linewidth=0.9)
    if obs is not None:
        ax.errorbar(ob_times, obs, yerr=err, fmt='o', alpha=0.7, color=palette[9], label='observations',
                    markeredgecolor='k', markeredgewidth=1.0, ms=8)
    ax.legend(loc=2)
    ax.set_xlabel('Date')
    ax.set_ylabel(ylab)
    ax.set_xlim([times[0],times[-1]])
    if axes is None:
        fig.autofmt_xdate()
    return ret_val


def calc_cosmos_stats(times, xb_dir, xa_dir,
                      cosmos_nc='/group_workspaces/jasmin2/hydro_jules/data/epinnington/smap_9km_uk/cosmos_daily_2016.nc'):
    sns.set_context('poster')
    sns.set_style('whitegrid')
    palette = sns.color_palette("colorblind", 11)
    site_arr = np.array([u'ALIC1', u'BALRD', u'BICKL', u'BUNNY', u'CARDT', u'CHIMN', u'CHOBH', u'COCLP',
                         u'CRICH', u'EASTB', u'ELMST', u'EUSTN', u'GISBN', u'GLENS', u'GLENW', u'HADLW', u'HARTW',
                         u'HARWD',
                         u'HENFS', u'HILLB', u'HOLLN', u'LIZRD', u'LODTN', u'LULLN', u'MOORH', u'MORLY', u'NWYKE',
                         u'PLYNL',
                         u'PORTN', u'RDMER', u'REDHL', u'RISEH', u'ROTHD', u'SHEEP', u'SOURH', u'SPENF', u'STGHT',
                         u'STIPS',
                         u'TADHM', u'WADDN', u'WYTH1'])
    xb_ens = extract_ens_obs(xb_dir, var='smcl', lvl_idx=0)
    xa_ens = extract_ens_obs(xa_dir, var='smcl', lvl_idx=0)
    xb_mean, up_xb, low_xb = calc_mean_upp_low(xb_ens)
    xa_mean, up_xa, low_xa = calc_mean_upp_low(xa_ens)

    cosmos = nc.Dataset(cosmos_nc, 'r')
    probe_idx = np.where(site_arr == 'CARDT')[0]
    cosmos_times = nc.num2date(cosmos.variables['time'][:], cosmos.variables['time'].units)
    cosmos_sm = cosmos.variables['smc'][:, probe_idx]
    fig, ax = plt.subplots(nrows=1, ncols=1)
    jules_stat_idx = np.array([find_nearest_idx_tol(times, cosmos_t) for cosmos_t in cosmos_times[cosmos_sm.nonzero()[0]]])
    jules_cosmos_sm_xb = xb_mean[jules_stat_idx]
    jules_cosmos_sm_xa = xa_mean[jules_stat_idx]
    cosmos_jules_sm = cosmos_sm[cosmos_sm.nonzero()[0]].flatten()
    xb_r = spst.stats.linregress(cosmos_jules_sm, jules_cosmos_sm_xb)[2]
    xa_r = spst.stats.linregress(cosmos_jules_sm, jules_cosmos_sm_xa)[2]
    xb_innov = [((cosmos_jules_sm[i] - np.mean(cosmos_jules_sm)) - (jules_cosmos_sm_xb[i] - np.mean(jules_cosmos_sm_xb)))**2 for i in xrange(len(jules_cosmos_sm_xb))]
    xa_innov = [
        ((cosmos_jules_sm[i] - np.mean(cosmos_jules_sm)) - (jules_cosmos_sm_xa[i] - np.mean(jules_cosmos_sm_xa))) ** 2
        for i in xrange(len(jules_cosmos_sm_xa))]
    xb_ubrmse = np.sqrt(np.sum(xb_innov) / len(jules_cosmos_sm_xb))
    xa_ubrmse = np.sqrt(np.sum(xa_innov) / len(jules_cosmos_sm_xa))
    xb_innov = [(cosmos_jules_sm[i] - jules_cosmos_sm_xb[i] )**2 for i in xrange(len(jules_cosmos_sm_xb))]
    xa_innov = [(cosmos_jules_sm[i] - jules_cosmos_sm_xa[i]) ** 2 for i in xrange(len(jules_cosmos_sm_xa))]
    xb_rmse = np.sqrt(np.sum(xb_innov) / len(jules_cosmos_sm_xb))
    xa_rmse = np.sqrt(np.sum(xa_innov) / len(jules_cosmos_sm_xa))
    ax.plot(times, xb_mean, '-', color=palette[0], label='Prior, r='+str(np.around(xb_r, decimals=2))+', ubrmse='+str(np.around(xb_ubrmse, decimals=2))+', rmse='+str(np.around(xb_rmse, decimals=2)))
    #ax.fill_between(times, low_xb, up_xb, facecolor=palette[0],
    #                   alpha=0.3, linewidth=0.0)
    ax.plot(times, xa_mean, '-', color=palette[2], label='Posterior, r=' + str(np.around(xa_r, decimals=2)) + ', ubrmse=' + str(
        np.around(xa_ubrmse, decimals=2)) + ', rmse=' + str(np.around(xa_rmse, decimals=2)))
    #ax.fill_between(times, low_xa, up_xa, facecolor=palette[2],
    #                   alpha=0.3, linewidth=0.0)
    ax.plot(cosmos_times, cosmos_sm, 'o-', color=palette[9], label='Cosmos probe', markeredgecolor='k',
            markeredgewidth=1.0, ms=8)
    ax.legend(loc=0, frameon=True, fancybox=True, framealpha=0.5)
    ax.set_title('COSMOS comparison to JULES before/after DA at CARDT')
    ax.set_xlabel('Date')
    ax.set_ylabel('Soil moisture (m$^{3}$ m$^{-3}$)')
    ax.set_xlim([times[0],times[-1]])
    fig.autofmt_xdate()
    return fig, ax


def plot_dist(xa_ens, xb_ens, idx=0, x_true=None, axes=None, title=None):
    """
    Plots prior and posterior distributions given ensembles
    :param xa_ens: posterior ensemble (arr)
    :param xb_ens: prior ensemble (arr)
    :param idx: index of parameter (int)
    :param x_true: "true" parameter vector if performing a twin experiment !optional! (arr)
    :param axes: axis to plot on !optional! (obj)
    :param title: title to give plot !optional! (str)
    :return: figure and axis objects or just axis object, dependent on if axes arg is specified
    """
    if axes is not None:
        ax = axes
        ret_val = ax
    else:
        fig, ax = plt.subplots(nrows=1, ncols=1,)
        ret_val = fig, ax
    sns.set_context('poster')
    sns.set_style('whitegrid')
    palette = sns.color_palette('Greys', 9)
    if x_true is not None:
        ax.axvline(x_true[title], color='k', linestyle='--')
    weights = np.ones_like(xb_ens[:,idx]) / float(len(xb_ens[:,idx]))
    sns.distplot(xb_ens[:,idx], kde=True, color=palette[3], ax=ax, hist_kws={'weights': weights})
    sns.distplot(xa_ens[:,idx], kde=True, color=palette[6], ax=ax, hist_kws={'weights': weights})
    if title is not None:
        ax.set_title(title)
    return ret_val


def plot_mult_dist(xa_ens, xb_ens, p_keys, x_true=None):
    """
    Creates a subplot of all the parameter distributions for the example included on Github
    :param xa_ens: posterior ensemble (arr)
    :param xb_ens: prior ensemble (arr)
    :param p_keys: list of order of parameters in prior and posterior ensemble (lst)
    :param x_true: "true" parameter dictionary !optional! (dict)
    :return: figure and axis objects
    """
    fig, ax = plt.subplots()
    ax1 = plt.subplot2grid(shape=(2,8), loc=(0,0), colspan=2)
    ax2 = plt.subplot2grid((2,8), (0,2), colspan=2)
    ax3 = plt.subplot2grid((2,8), (0,4), colspan=2)
    ax4 = plt.subplot2grid((2,8), (0,6), colspan=2)
    ax5 = plt.subplot2grid((2,8), (1,0), colspan=2)
    ax6 = plt.subplot2grid((2,8), (1, 2), colspan=2)
    ax7 = plt.subplot2grid((2,8), (1, 4), colspan=2)
    ax8 = plt.subplot2grid((2,8), (1, 6), colspan=2)
    ax1 = plot_dist(xa_ens, xb_ens, 0, x_true, axes=ax1, title=p_keys[0])
    ax2 = plot_dist(xa_ens, xb_ens, 1, x_true, axes=ax2, title=p_keys[1])
    ax3 = plot_dist(xa_ens, xb_ens, 2, x_true, axes=ax3, title=p_keys[2])
    ax4 = plot_dist(xa_ens, xb_ens, 3, x_true, axes=ax4, title=p_keys[3])
    ax5 = plot_dist(xa_ens, xb_ens, 4, x_true, axes=ax5, title=p_keys[4])
    ax6 = plot_dist(xa_ens, xb_ens, 5, x_true, axes=ax6, title=p_keys[5])
    ax7 = plot_dist(xa_ens, xb_ens, 6, x_true, axes=ax7, title=p_keys[6])
    ax8 = plot_dist(xa_ens, xb_ens, 7, x_true, axes=ax8, title=p_keys[7])
    fig.subplots_adjust(wspace=1.7, hspace=0.3)
    return fig, ax


def save_plots(xa_ens_pickle, out_dir):
    """
    Function saving plots from data assimilatoin experiment output
    :param xa_ens_pickle: location of pickled posterior ensemble arr (str)
    :param out_dir: directory to save plots in (str)
    :return: string confirming plots have been saved (str)
    """
    dat_xt = nc.Dataset('output/background/xb0.daily.nc', 'r')
    jda = fourdenvar.FourDEnVar()
    obs = es.obs_fn()
    date = nc.num2date(dat_xt.variables['time'][:], dat_xt.variables['time'].units)
    xb_dir = es.output_directory + '/ensemble' + str(es.seed_value)
    xa_dir = es.output_directory + '/ensemble_xa_' + str(es.seed_value)
    #plot SM
    fig, ax = plot_twin_spread(date[:], 'smcl', xb_dir, xa_dir,
                             ob_times=obs['sm_times'], obs=obs['sm_obs'], err=obs['sm_err'],
                             ylab=r'Soil Moisture (m$^{3}$ m$^{-3}$)')
    fig.savefig(out_dir+'/sm.png', bbox_inches='tight')
    xb_dir2 = es.output_directory + '/background_tst'
    xa_dir2 = es.output_directory + '/analysis_med'
    fig, ax = calc_cosmos_stats(date[:], xb_dir2, xa_dir2)
    fig.savefig(out_dir+'/cosmos_median.png', bbox_inches='tight')
    #plot distribution
    #true_params = {'alpha_io': 5.5e-02, 'neff_io': 5.7e-04, 'fd_io': 9.6e-03, 'mu_io': 2.0e-02, 'nu_io': 4.0e+00,
    #               'gamma_io': 1.76e+01, 'delta_io':-3.3e-01}
    xa_ens = pickle.load(open(xa_ens_pickle, 'rb'))
    p_keys = ['oneovernminusone','oneoveralpha','satcon','vsat','vcrit','vwilt','hcap','hcon']
    fig, ax = plot_mult_dist(xa_ens, jda.xbs, p_keys)
    fig.savefig(out_dir + '/distributions.png', bbox_inches='tight')
    return 'plots saved!'