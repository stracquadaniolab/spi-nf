import numpy as np

################################################################################
###
### Configuration class
###     Basic superclass for MC configuration
###
################################################################################
class MonteCarloConfiguration:

    def __init__(self, rng = None):
        self.trajectory = []
        self.rng = rng

################################################################################
###
### GenomeConfiguration class
### It represents a genome structure as follows:
### 1. segment_size: length in bp of each segment in 0-based notation
### 2. segment_type: type of segment in 0-based notation:
###     0 = Non essential
###     1 = Essential
###     2 = Marker
###
################################################################################
class GenomeConfiguration(MonteCarloConfiguration):
    #
    #   default constructor
    #       it takes in input a genome object and initizialize all the
    #       instance variables.
    #
    #       genome: a Genome object
    #       genotype: the list of segments in the current
    #       trajectory: the series of steps that leads to the current config
    #       essential_deleted: number of essential deleted
    #       marker_deleted: number of markers lost
    #           Note:  they are flag to keep track of lost markers or
    #                  essentials genes.
    #                  This avoids to check the genotype each time we have
    #                  to see if the genome is valid.
    #
    def __init__(self, genome, rng = None):
        self.rng = rng
        self.genome = genome
        self.genotype = np.array(range(1,len(genome)+1), dtype=np.intc)
        self.trajectory = []
        self.essential_deleted = False
        self.marker_deleted = False

    #
    #   returns a string representation of the genotype
    #
    def __repr__(self):
        return str(self.genotype)

    #
    #   method to check if the configuration is valid or not.
    #   it just checks the flags.
    #
    def is_valid(self):
        return (not self.essential_deleted) and (self.marker_deleted)

    def genotype_length(self):
        gt_len = np.sum([self.genome.segment_size[np.abs(s)-1] for s in self.genotype])
        return gt_len

    def site_distance(self, start, stop):
        dist = 0
        itr = start
        while itr != stop:
            curr_seg = np.abs(self.genotype[itr])-1
            dist += self.genome.segment_size[curr_seg]
            # updating the itr var, keep in mind that the genome is circular.
            itr = (itr + 1) % self.genotype.size
        return dist

    #
    #   method to check if the configuration is valid or not.
    #   it just checks the flags.
    #
    def is_viable(self):
        return not self.essential_deleted

    #
    #   perform an inversion or deletion between two segments
    #   based on the stochastic matrix passed in input
    #
    def perform_recombination(self, prob_matrix):
        # generating a coin
        coin = self.rng.random()
        # norm_prob_matrix = prob_matrix.copy() 
        # norm_prob_matrix = norm_prob_matrix / np.sum(norm_prob_matrix)

        # probabilisticly picking a recombination using the coin random value
        i, j = divmod(prob_matrix.cumsum().searchsorted(coin), prob_matrix.shape[1])

        # recording the type of move taken
        event = None
        site_dist = self.site_distance(i,j)
        
        # picking uniform at random between inversion and deletion
        coin = self.rng.random()
        if coin < 0.5:
            # picking inversion
            event = 'i'
            self.__perform_inversion(i,j)
        else:
            # picking deletion
            event = 'd'
            self.__perform_deletion(i,j)
        
        # computing genome length after the recombination happened
        gt_len = self.genotype_length()

        # saving the trajectory that leads to the current state
        self.trajectory.append((event, i,j, prob_matrix[i,j], gt_len, site_dist))

    #
    #   perform a specific recombination based on a trajectory tuple
    #   1. event type
    #   2. first breakpoint
    #   3. second breakpoint
    #   4. probability of recomb (just recorded for compatibility)
    #
    def perform_deterministic_recombination(self, recomb):
        event, i, j, prob = recomb
        if event == 'i':
            self.__perform_inversion(i,j)
        else:
            self.__perform_deletion(i,j)

    #
    #   perform a deletion between `start` and `stop`
    #
    def __perform_deletion(self, start, stop):
        # keeping track of the elements to delete
        del_id = []

        # flags to keep track whether an essential or marker was deleted
        essential_deleted = False
        marker_deleted = False

        # iterating from start to stop
        itr = start

        while itr != stop:
            # tracking the indexes to delete
            del_id.append(itr)

            # getting the segment id in 0-based notation regardless of inv.
            curr_segment = np.abs(self.genotype[itr])-1

            # checking if curr_segment is an essential gene
            if self.genome.segment_type[curr_segment] == 1:
                essential_deleted = True

            # checking if curr_segment is a marker gene
            if self.genome.segment_type[curr_segment] == 2:
                marker_deleted = True

            # updating the itr var, keep in mind that the genome is circular.
            itr = (itr + 1) % self.genotype.size

        # using numpy delete to remove the elements form the current genotype.
        self.genotype = np.delete(self.genotype, del_id)

        # updating the flags for deleted marker and essential genes.
        self.essential_deleted = self.essential_deleted or essential_deleted
        self.marker_deleted = self.marker_deleted or marker_deleted

    #
    #   perform an inversion between `start` and `stop`
    #
    def __perform_inversion(self, start, stop):
        # keeping track of the segments to invert
        inv_id = []

        # iterating from start to stop
        itr = start

        while itr != stop:
            # adding the itr-id to the inv_id
            inv_id.append(itr)

            # updating the itr var, keep in mind that the genome is circular.
            itr = (itr + 1) % self.genotype.size

        # gettin the subset of segments to invert
        inv_seg = self.genotype[inv_id]

        # the inv_id segments will be inverted
        self.genotype[inv_id] = (inv_seg[::-1] * -1)
