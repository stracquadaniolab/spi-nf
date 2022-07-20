import os
import sys
import time
import numpy as np
import ctypes
from collections import defaultdict

from .MonteCarloConfiguration import *

matrix_utils = ctypes.cdll.LoadLibrary("/opt/libspi/lib/spi_matrix_utils.so")

################################################################################
###
### MonteCarloSampler class
###     Basic superclass for the MC algorithm
###
################################################################################
class MonteCarloSampler:
    #
    # default constructor which takes a genome to sample from.
    #
    def __init__(self, genome):
        self.genome = genome
        self.rng_seed = int(round(time.time()))
        self.rng = np.random.default_rng(self.rng_seed)


    ##############################################################
    ###### simple function to warmup the RNG
    ##############################################################
    def burnin_rng(self):
        for burnin in range(100000):
            self.rng.random()

    #
    # compute stochastic matrix using the C library for the genotype.
    #   genotype: the genotype to analyze
    #   nu: the 3D constant
    #   b: the persistence length
    #
    def compute_stochastic_matrix(self, genotype, nu, b):
        # pre-allocating the matrix of n_segments * n_segments
        pmat = np.zeros((genotype.size, genotype.size), dtype=np.float64)

        # getting pointers for all the data to pass in input
        pmat_ptr = pmat.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        g_ptr = genotype.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        s_ptr = self.genome.segment_size.ctypes.data_as(ctypes.POINTER(ctypes.c_int))

        # compute the stochastic matrix using the C library.
        # REMARK: you have to pass the size of each element as in pure C-style
        matrix_utils.compute_stochastic_matrix(pmat_ptr, pmat.shape[0],
                                           g_ptr, genotype.size,
                                           s_ptr, self.genome.segment_size.size,
                                           ctypes.c_double(nu),
                                           ctypes.c_double(b))

        # return the matrix
        return pmat

################################################################################
###
### GridMonteCarloSampler class
###     Perform Montecarlo Samping on a grid of parameters
###
################################################################################
class GridMonteCarloSampler(MonteCarloSampler):

    #
    # simply calling the superclass constructor
    #
    def __init__(self, genome):
        MonteCarloSampler.__init__(self, genome)

    #
    # perform a Montecarlo simulations for each parameter in the paramter grid
    #   max_configuration: number of sample to sample
    #   param_set: set of parameters to use for MC simulations.
    #
    def run(self, max_configuration, param_set):
        self.burnin_rng()

        #
        start_time = time.time()
        #
        # keeping track of the trajectories that leads to feasible configs.
        # It is a simple dict indexed by the param config and has a list of list.
        #
        trajectory_pool = dict()

        #
        # keeps track of the oberved segment count normalized [0,1]
        # It is a simple dict indexed by the param config pointing to np.array
        #
        trajectory_profile = dict()

        # iterate through the parameter set.
        for param_i in param_set:
            # get three parameters form the param_i
            lam, nu, b = param_i

            # logging the parameters tested
            print("\nSimulating lambda=%f nu=%f b=%f" % (lam, nu, b))

            # initialize the pool for param_i
            trajectory_pool[param_i] = []

            # initialize the trajectory_profile for param_i
            trajectory_profile[param_i] = np.zeros(len(self.genome))

            # keeping track of configurations generated
            config_generated = 0
            config_accepted = 0

            while config_accepted < max_configuration:
                config_generated += 1
                # creating a new config
                mc_config = GenomeConfiguration(self.genome, self.rng)

                # sampling from Poisson the number of recombinations
                n_recomb = self.rng.poisson(lam)

                # performing the n_recomb
                for n_r in range(n_recomb):

                    # computing the stochastic matrix for recombination
                    prob = self.compute_stochastic_matrix(mc_config.genotype, nu, b)

                    # performing the recombination
                    mc_config.perform_recombination(prob)

                    # keep performing the recombinations if the config is still viable.
                    if not mc_config.is_viable():
                        break

                # update statistics only if the config is valid
                if mc_config.is_valid():
                    # updating the count of valid configs
                    config_accepted += 1

                    # keeping track of the valid trajectories
                    trajectory_pool[param_i].append(mc_config)

                    # updating the frequencies
                    for seg in mc_config.genotype:
                        trajectory_profile[param_i][np.abs(seg)-1] += 1

                    if config_generated % 100 == 0:
                        sys.stdout.write("> %d of %d configurations generated...\r" % (config_generated, max_configuration))
                        sys.stdout.flush()

            # computing frequencies for param_i
            trajectory_profile[param_i] = trajectory_profile[param_i] / float(max_configuration)

        return (trajectory_pool, trajectory_profile, (config_accepted, config_generated, time.time()-start_time, self.rng_seed))


################################################################################
###
### ReweightGridMonteCarloSampler class
###     peform MC reweighting using the trajectories generated by
###     GridMonteCarloSampler.
###
################################################################################
class ReweightGridMonteCarloSampler(MonteCarloSampler):

    #
    # Default constructor that takes in input the genome structure
    #
    def __init__(self, genome, rng = None):
        MonteCarloSampler.__init__(self, genome)
    #
    # peform reweighting on param_grid using trajectory_pool
    #   trajectory is a dict with a list of trajectories indexed by param.
    #
    def run(self, config_pool, reweighting_grid):
        #
        self.burnin_rng()

        #
        # we have to keep track of the weights per segment and the total sum
        # which should be the partition function over all the possible configs.
        #
        weight_grid = defaultdict(float)
        weight_segment = defaultdict(lambda:np.zeros(len(self.genome)))

        # looping through the coase grained grid and
        # reweighting over the finer mesh
        for params, configs in config_pool.items():
            # parameters from the coarse grain model
            lam, nu, b = params

            print("Reweighting using trajectories generated with %d, %f, %f" % (lam, nu, b))
            sys.stdout.flush()

            # keeping track of the start time
            start_time = time.time()

            for cnt, cgt in enumerate(configs):

                genotype, trajectory = cgt

                if cnt % 100 == 0: 
                    print("\t Processed %d samples" % cnt)

                # keep track of the weight for each trajectory
                ln_wt = defaultdict(float)

                # current trajectory
                curr_config = GenomeConfiguration(self.genome,self.rng)

                # adding weight for the number of recomb
                n_recomb = len(trajectory)

                # adding weight for the difference in numb of recomb
                for r_params in reweighting_grid:
                    rlam, rnu, rb = r_params
                    ln_wt[r_params] = n_recomb * np.log(rlam / lam) - rlam + lam
                    
                # for each move in the trajectory compute reweighting weights
                for t_e, t_i, t_j, t_prob, g_len, dist in trajectory:
                    for r_params in reweighting_grid:
                        rlam, rnu, rb = r_params
                        
                        # print "reweighting for %f %f %f" %(rlam, rnu, rb)
                        
                        # compute stochastic matrix
                        rw_mat = self.compute_stochastic_matrix(curr_config.genotype, rnu, rb)
                        
                        # reweight prob for (i,j)
                        rw_prob = rw_mat[t_i, t_j]

                        # reweight ratio
                        rw_wt = np.log(rw_prob) - np.log(t_prob)

                        # print(rw_prob, t_prob, rw_wt)
                        
                        # updating the weight for the r_params based on r_wt
                        ln_wt[r_params] += rw_wt
                                            
                    curr_config.perform_deterministic_recombination((t_e, t_i, t_j, rw_prob))

                # update the global counts based on the trajectory processed
                for r_params in reweighting_grid:
                    
                    # updating global table
                    wt = np.exp(ln_wt[r_params])
                    weight_grid[r_params] += wt

                    # assigning weights to segments
                    for seg in genotype:
                        weight_segment[r_params][np.abs(seg)-1] += wt

            print("\tElapse time: %.3f sec." % (time.time() - start_time))

        # computing the frequencies based on the reweight vectors
        results = dict()
        for k, v in weight_segment.items():
            print("effective number of reweigthed samples for %f,%f,%f: %g" % (k[0], k[1], k[2], weight_grid[k]))
            results[k] = v / weight_grid[k]

        return results
