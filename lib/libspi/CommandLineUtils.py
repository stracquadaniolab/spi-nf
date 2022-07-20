################################################################################
###
###     Module to generate parameters grids and parse config files.
###
################################################################################

import json
import numpy as np


##############################################################
###### load a JSON configuration file
##############################################################
def load_configuration_file(filename):
    return json.load(open(filename))

##############################################################
###### generate `step` sized vectors from `low` to `high`
##############################################################
def param_grid(low, high, step):
    return np.linspace(low, high, np.round(((high-low)/step))+1.0)


##############################################################
###### generate parameter grid
##############################################################
def build_parameters_grid(params):
    # grid definition
    lam_vec = param_grid(params['lambda']['min'],params['lambda']['max'],params['lambda']['step'])
    nu_vec = param_grid(params['nu']['min'],params['nu']['max'],params['nu']['step'])
    b_vec = param_grid(params['b']['min'],params['b']['max'],params['b']['step'])
    grid = [(l, n, b) for l in lam_vec for n in nu_vec for b in b_vec ]
    return grid
