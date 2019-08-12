#######################################
# File to specify setup of experiment #
#######################################
# core py modules:
import os
# local modules:
import observations
import plot


# set directory for JULES output
output_directory = os.getcwd()+'/output'
# set directory containing JULES NML files
nml_directory = os.getcwd()+'/chess_da_cardt'
# set model executable
model_exe = '/home/users/ewanp82/models/jules5.3/build/bin/jules.exe'
# set function to extract JULES modelled observations for prior JULES
jules_hxb = observations.extract_jules_hxb
jules_hxi = observations.extract_jules_hx
# set function to extract prior ensemble of modelled observations
jules_hxb_ens = observations.extract_jules_hxb_ens
# set function to extract observations to be assimilated
obs_fn = observations.extract_twin_data
# set JULES parameters to optimised during data assimilation
opt_params = {'ancillaries': {
                  'jules_soil_props': {
                      'const_val': [[0, 9.881425, (0., 20.)],
                                    [1, 0.27247956, (0.0, 0.99)],
                                    [2, 0.0287, (0.0, 0.099)],
                                    [3, 0.51, (0.0, 0.99)],
                                    [4, 0.39324558, (0.0, 0.99)],
                                    [5, 0.26872617, (0.0, 0.99)],
                                    [6, 951580., (450000., 1500000.)],
                                    [7, 0.16550983, (0.0, 0.99)]],
                }}}
# set error on prior parameter estimates
prior_err = 0.25
# set size of ensemble to be used in data assimilation experiments
ensemble_size = 100
# set number of processors to use in parallel runs of JULES ensemble
num_processes = 25
# set seed value for any random number generation within experiments
seed_value = 0
# plotting save function
save_plots = plot.save_plots
# plotting output director
plot_output_dir = os.getcwd()+'/output/plot'

