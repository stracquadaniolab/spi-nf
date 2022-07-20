import numpy as np

################################################################################
###
### Genome class
### It represents a genome structure as follows:
### 1. segment_size: length in bp of each segment in 0-based notation
### 2. segment_type: type of segment in 0-based notation:
###     0 = Non essential
###     1 = Essential
###     2 = Marker
###
################################################################################
class Genome:
    #
    #   default constructor
    #       it just initialize everything to None or zero
    #
    def __init__(self):
        self.length = 0
        self.segment_size = None
        self.segment_type = None

    #
    #   Return the number of segments in the genome
    #
    def __len__(self):
        return self.length

    #
    #   load_genome_from_file
    #       Load a genome structure file and set the relevant instance variables.
    #       Each record has three tab-separated fields:
    #       1. segment number
    #       2. segment length in BP
    #       3. segment type
    def load_genome_from_file(self,filename):
        #   local variables to store size, essential and marker segments.
        size_l = []
        ess_l = []
        mark_l = []

        #
        # reading the structure file line by line
        # WARN: exeception to parsing should be added
        #
        for line in open(filename):

            # skipping comments
            if line.startswith("#"):
                continue

            # parsing the three fields
            s_id, s_len, s_type = line.strip().split("\t")

            # saving segment length
            size_l.append(float(s_len))

            # checking the segment type
            if s_type == "ESSENTIAL":
                ess_l.append(int(s_id)-1)

            # checking the segment type
            if s_type == "MARKER":
                mark_l.append(int(s_id)-1)

        #
        # saving the size of the genome and the number of segments
        self.segment_size = np.array(size_l, dtype=np.float64)
        self.length  = self.segment_size.size

        # the array segment type is used as a mask to check segment status
        # it is faster than doing `in`.
        self.segment_type = np.zeros(self.length)
        self.segment_type[ess_l] = 1
        self.segment_type[mark_l] = 2
