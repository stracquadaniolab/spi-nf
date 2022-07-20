################################################################################
### Module to generate output files for the different types of data returned
### by the sampler
################################################################################
import time
import datetime
import numpy as np

def __get_timestamp_header():
    timestamp_str = datetime.datetime.now().strftime("%H:%M:%S %d/%m/%Y")
    timestamp_str = "#\n#\tBuildtime: %s\n#\n" % timestamp_str
    return timestamp_str

def __get_arguments_header(args):
    args_str = "#\n"
    for k,v in args.items():
        args_str =  args_str + "#\t %s: %s\n" % (k, str(v))
    args_str = args_str + "#\n"
    return args_str



################################################################################
####
#### write genome segment structure to file
####    filename: name of the file
####    distance matrix: the distance matrix between loxp sites
####
####    File format:
####        segment_id TAB length TAB type
####        segment_id: integer
####        length: integer
####        type: string ESSENTIAL|MARKER|NONESSENTIAL
####
################################################################################
def save_genome_structure_to_file(filename, distance_matrix, essential_list,
                                  marker_list, args):
    fh = open(filename, 'w')
    # timestamping the file
    fh.write(__get_timestamp_header())
    # adding the list of params used to generate the file
    fh.write(__get_arguments_header(args))

    n_site = distance_matrix.shape[0]
    for i in range(n_site):
        segment_type = None
        if (i+1) in essential_list:
            segment_type = "ESSENTIAL"
        elif (i+1) in marker_list:
            segment_type = "MARKER"
        else:
            segment_type = "NONESSENTIAL"

        fh.write("%d\t%d\t%s\n" % (i+1, distance_matrix[i, (i+1)%43], segment_type))

    fh.close()

################################################################################
####
#### write genome segment profiles to file
####    filename: name of the file
####    profile: a dict with parameters grid as keys and np.array for the freq
####
####    File format:
####        params TAB freqs
####        params: comma separated
####        freqs : comma separated
####
################################################################################
def save_profiles_to_file(filename, profile, args):
    fh = open(filename, 'w')

    # timestamping the file
    fh.write(__get_timestamp_header())
    # adding the list of params used to generate the file
    fh.write(__get_arguments_header(args))

    for k,v in profile.items():
        fh.write("%f,%f,%f\t%s\n" % (k[0], k[1], k[2], ",".join(map(str,v))))
    fh.close()

################################################################################
#### load genome segment profiles to file
################################################################################
def load_profiles_from_file(filename):
    profiles = dict()

    for line in open(filename):
        if line.startswith("#"):
            continue

        params, profile = line.strip().split("\t")

        # converting everything to float
        lam, nu, b = list(map(float, params.split(",")))

        # converting everything to float
        profile = list(map(float, profile.split(",")))

        # updating profiles
        profiles[(lam, nu, b)] = np.array(profile)

    return profiles


################################################################################
####
#### write config trajectories to file
####    filename: name of the file
####    profile: a dict with parameters grid as keys and list containing the
####             list of each move for each config
####
####    File format:
####        params TAB trajectories
####        trajectories: traj_1 PIPE traj_2 PIPE ... PIPE traj_N
####        traj: move_1 COMMA move_2 COMMA ... COMMA move_m
####        move: comma separated
####
################################################################################
def save_trajectories_to_file(filename, pool, args):
    # opening file for writing
    fh = open(filename, 'w')

    # timestamping the file
    fh.write(__get_timestamp_header())

    # adding the list of params used to generate the file
    fh.write(__get_arguments_header(args))

    # looping through trajectories per config
    for params, configs in pool.items():
        # keeping a list of sample
        # REMARK: I could have done it more compact but I preferred to keep it
        #           clean
        config_list = []

        for t in configs:
            t_str = ";".join([",".join(map(str, m)) for m in t.trajectory])
            gt_str = ",".join(map(str,t.genotype)) + ";" + t_str
            config_list.append(gt_str)

        # joining samples together
        s_str = "|".join(config_list)

        # joining params together
        p_str = ",".join(map(str, params))

        # writing record to file
        fh.write("%s\t%s\n" % (p_str, s_str))
    fh.close()

################################################################################
#### write trajectories to file
################################################################################
def load_trajectories_from_file(filename):
    pool = dict()
    for line in open(filename):
        if line.startswith("#"):
            continue

        params, config_list = line.strip().split("\t")

        lam, nu, b = list(map(float, params.split(",")))
        param_tuple = (lam, nu, b)
        pool[param_tuple] = []

        for t in config_list.split("|"):
            fields = t.split(";")
            genotype = list(map(int,fields[0].split(",")))
            trajectory = []
            for m in fields[1:]:
                tfields = m.split(",")
                m_p = (tfields[0], int(tfields[1]), int(tfields[2]), float(tfields[3]), float(tfields[4]), float(tfields[5]))
                trajectory.append(m_p)
            pool[param_tuple].append((genotype,trajectory))

    return pool

def save_trajectories_stats_to_file(in_filename, out_filename):

    # open stats file
    fh = open(out_filename, 'w')
    header = "lambda,nu,b,genome_length, num_events, num_segments,num_unique_segments\n"
    fh.write(header)

    for line in open(in_filename):
        if line.startswith("#"):
            continue

        params, config_list = line.strip().split("\t")

        lam, nu, b = list(map(float, params.split(",")))
        param_tuple = (lam, nu, b)

        for t in config_list.split("|"):
            fields = t.split(";")
            genotype = list(map(int,fields[0].split(",")))
            trajectory = []
            for m in fields[1:]:
                tfields = m.split(",")
                m_p = (tfields[0], int(tfields[1]), int(tfields[2]), float(tfields[3]), float(tfields[4]), float(tfields[5]))
                trajectory.append(m_p)

            genome_length = trajectory[-1][-2]
            num_events = len(trajectory)
            num_segments = len(genotype)
            num_unique_segments = len(np.unique(np.abs(genotype)))

            record = "%f,%f,%f,%d,%d,%d,%d\n" %(lam, nu, b, genome_length, num_events, num_segments, num_unique_segments)
            fh.write(record)
    fh.close()

################################################################################
#### save profile evaluation to file
################################################################################
def save_profile_evaluation_to_file(filename, results, args):
    fh = open(filename, 'w')

    # timestamping the file
    fh.write(__get_timestamp_header())
    # adding the list of params used to generate the file
    fh.write(__get_arguments_header(args))

    fh.write("lambda,nu,b,loglik\n")
    params = sorted(results.keys())
    for k in params:
        fh.write("%f,%f,%f,%f\n" % (k[0], k[1], k[2], results[k]))
    fh.close()

################################################################################
#### save run information to file
################################################################################
def save_run_information_to_file(filename, results, args):
    fh = open(filename, 'w')

    # timestamping the file
    fh.write(__get_timestamp_header())
    # adding the list of params used to generate the file
    fh.write(__get_arguments_header(args))

    fh.write("# lambda,nu,b,accepted,generated,elapsed_time,seed\n")
    fh.write("%f,%f,%f,%d,%d,%.3f,%d\n" % (args['lb'], args['nu'], args['b'], results[0], results[1], results[2], results[3]))
    fh.close()
