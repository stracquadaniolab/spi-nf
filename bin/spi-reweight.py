#!/usr/bin/env python
import os
import argparse
import json
import time

import numpy as np

from libspi.Genome import *
from libspi.MonteCarloSampler import *
import libspi.IO as IO
import libspi.CommandLineUtils as cli

##############################################################
# CMD line options parser
##############################################################


def parse_command_line_options():
    parser = argparse.ArgumentParser(
        description='Monte Carlo reweighting for Cre-Lox Recombination studies.')
    parser.add_argument('--input-file', '-i', type=str,
                        help="Chromosome file structure. Chromosome is assumed to be circular.")
    parser.add_argument('--output-file', '-o', type=str,
                        required=True, help="Output file")
    parser.add_argument('--trajectory-file', '-t', type=str,
                        required=True, help="Trajectory file")
    # parser.add_argument('--config', '-c', type=str, required=True, help="Simulation configuration file")
    parser.add_argument('--lb', type=float, required=True,
                        help="Lambda parameter")
    parser.add_argument('--nu', type=float, required=True, help="Nu parameter")
    parser.add_argument('--b', type=float, required=True, help="B parameter")
    parser.add_argument('--rl', type=float, required=True,
                        help="Reweighting radius for lambda parameter")
    parser.add_argument('--rnu', type=float, required=True,
                        help="Reweighting radius for nu parameter")
    parser.add_argument('--rb', type=float, required=True,
                        help="Reweighting radius for b parameter")
    options = parser.parse_args()
    return options


##############################################################
# MAIN method
##############################################################
if __name__ == "__main__":
    #
    options = parse_command_line_options()
    # config = cli.load_configuration_file(options.config)

    # building the parameters grid
    # param_grid = cli.build_parameters_grid(config)
    param_grid = [(options.lb, options.nu, options.b)]

    print("Using %d parameters settings." % len(param_grid))

    # creating and loading genome structure file
    genome = Genome()
    genome.load_genome_from_file(options.input_file)

    # loading trajectory pool
    trajectory_pool = IO.load_trajectories_from_file(options.trajectory_file)
    filtered_pool = dict()

    for curr_param, curr_traj in trajectory_pool.items():
        c_l, c_nu, c_b = curr_param
        if (np.abs(c_l - options.lb) <= options.rl) and \
           (np.abs(c_nu - options.nu) <= options.rnu) and \
           (np.abs(c_b - options.b) <= options.rb):
            filtered_pool[curr_param] = curr_traj

    print("Reweighting using %d simulations." % len(filtered_pool.keys()))
    # creating the MC sampler
    mc = ReweightGridMonteCarloSampler(genome)

    # running and timing the sampler
    # t_start = time.time()
    trajectory_profile = mc.run(filtered_pool, param_grid)
    # t_elapsed = time.time() - t_start

    # saving objects to file
    IO.save_profiles_to_file(
        options.output_file, trajectory_profile, vars(options))
