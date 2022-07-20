#!/usr/bin/env python
import argparse
import re
import numpy as np
from Bio import SeqIO

import spi.IO as IO

def map_sites(reference, site_seq="ATAACTTCGTATAATGTACATTATACGAAGTTAT"):
    site_map = []
    site_regex = re.compile(site_seq)
    for m in site_regex.finditer(str(reference.seq)):
        site_map.append(m)
        #print m.start(), " ", m.end()
    return site_map

def compute_distance_matrix(reference, site_map):
    distance_matrix = np.zeros((len(site_map),len(site_map)))
    len_genome = len(reference)

    for i in range(len(site_map)):
        for j in range(i+1,len(site_map)):
            distance_matrix[i,j] = np.abs(site_map[i].start()-site_map[j].start())
            distance_matrix[j,i] = len_genome - distance_matrix[i,j]

    return distance_matrix


##############################################################
###### CMD line options parser
##############################################################
def parse_command_line_options():
    parser = argparse.ArgumentParser(description='Compute chromosome structure file, which provides the length and the type of each segment bounded by loxp-sites.')
    parser.add_argument('--input-file', '-i', type=str, help="FASTA file. Chromosome is assumed to be circular.")
    parser.add_argument('--output-file', '-o', type=str, required=True, help="Output structure file")
    parser.add_argument('--marker', '-m', type=int, nargs='+', help="Segments to be reported as markers")
    parser.add_argument('--essential', '-e', type=int, nargs='+', help="Segments to be reported as essential")
    options = parser.parse_args()
    return options

if __name__ == "__main__":
    options = parse_command_line_options()
    reference = SeqIO.read(options.input_file, "fasta")
    site_map = map_sites(reference)

    distance_matrix = compute_distance_matrix(reference, site_map)

    IO.save_genome_structure_to_file(options.output_file, distance_matrix, options.essential, options.marker, vars(options))
