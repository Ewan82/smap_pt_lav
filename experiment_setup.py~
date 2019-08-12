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
# set function to extract prior ensemble of modelled observations
jules_hxb_ens = observations.extract_jules_hxb_ens
# set function to extract observations to be assimilated
obs_fn = observations.extract_twin_data
# set JULES parameters to optimised during data assimilation
opt_params = {'ancillaries': {
                  'jules_soil_props': {
                      'const_val': [[0, 9.881425, (3., 15.)],
                                    [1, 0.27247956, (0.05, 0.9)],
                                    [2, 0.0287, (0.005, 0.09)],
                                    [3, 0.51, (0.05, 0.9)],
                                    [4, 0.39324558, (0.05, 0.9)],
                                    [5, 0.26872617, (0.05, 0.9)],
                                    [6, 951580., (750000., 1200000.)],
                                    [7, 0.16550983, (0.005, 0.9)]],
                }}}
# set error on prior parameter estimates
prior_err = 0.25
# set size of ensemble to be used in data assimilation experiments
ensemble_size = 100
# set number of processors to use in parallel runs of JULES ensemble
num_processes = 100
# set seed value for any random number generation within experiments
seed_value = 0
# plotting save function
save_plots = plot.save_plots
# plotting output director
plot_output_dir = os.getcwd()+'/output/plot'

