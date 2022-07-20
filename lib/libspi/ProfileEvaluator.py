import numpy as np
import scipy.stats as stats

################################################################################
###
### Model Evaluator class
###     Basic superclass to implement model fitting evaluation
###
################################################################################
class ProfileEvaluator:

    #
    #   constructor that takes a genome in input and set the other variables to
    #   None or zero
    #
    def __init__(self, genome):
        self.genome = genome
        self.traning_set = None
        self.traning_set_len = 0

    def load_training_dataset_from_file(self, filename):
        # reading the file with training dataset.
        self.traning_set = np.zeros(len(self.genome))
        self.traning_set_len = 0

        # file is assumed in reconstruction format
        # format:
        #   1 genome per line
        #   segments separated by commas
        for line in open(filename):
            segs = line.strip().split(",")
            segs = list(map(int, segs))
            segs = np.array(segs, np.intc)
            segs = np.unique(np.abs(segs))
            segs = segs[segs <= len(self.genome)]
            segs = segs - 1

            self.traning_set[segs] += 1
            self.traning_set_len += 1

################################################################################
###
### ModelLogLik class
###     Compute the loglik of profiles given a traning set
###
################################################################################
class BinomialLogLikEvaluator(ProfileEvaluator):
    #
    #   evaluates binomial loglik for a dict of profiles
    #   returns a dict where keys are the params set and value is the loglik
    #
    def evaluate(self, profile_set):
        results = dict()

        # computing loglik for all paramters settings
        for params, profile in profile_set.items():
            logLik = 0.0
            for si, sp in enumerate(profile):
                prob = stats.binom.pmf(self.traning_set[si], self.traning_set_len, sp)
                if prob > 0:
                    logLik += np.log(prob)
            results[params] = logLik

        return results
