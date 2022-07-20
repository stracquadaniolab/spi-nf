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
###### CMD line options parser
##############################################################
def parse_command_line_options():
    parser = argparse.ArgumentParser(description='The SCRaMbLE Polymer Interaction model (SPI) rejection sampling algorithm')
    parser.add_argument('--input-file', '-i', type=str, help="Chromosome file structure. Chromosome is assumed to be circular.")
    parser.add_argument('--output-file', '-o', type=str, required=True, help="Output file")
    parser.add_argument('--trajectory-file', '-t', type=str, required=True, help="Trajectory file")
    parser.add_argument('--lb', type=float, required=True, help="Lambda parameter")
    parser.add_argument('--nu', type=float, required=True, help="Nu parameter")
    parser.add_argument('--b', type=float, required=True, help="B parameter")
    parser.add_argument('--run-info', type=str, required=True, help="Run information file")
    parser.add_argument('--max-config', "-m", type=int, required=True, help="Maximum number of configurations")
    
    options = parser.parse_args()
    return options


##############################################################
###### MAIN method
##############################################################
if __name__ == "__main__":
    #
    options = parse_command_line_options()

    # creating and loading genome structure file
    genome = Genome()
    genome.load_genome_from_file(options.input_file)

    # creating the MC sampler
    mc = GridMonteCarloSampler(genome)

    # running and timing the sampler
    # t_start = time.time()
    model_parameters = (options.lb, options.nu, options.b)
    trajectory_pool, trajectory_profile, run_info = mc.run(options.max_config, [model_parameters])
    # t_elapsed = time.time() - t_start

    # saving objects to file
    IO.save_trajectories_to_file(options.trajectory_file, trajectory_pool, vars(options))
    IO.save_profiles_to_file(options.output_file, trajectory_profile, vars(options))
    IO.save_run_information_to_file(options.run_info, run_info, vars(options))