#! /usr/bin/env python

import sys, subprocess, os
import utils

def main(args):

    hout = open(args.output_file + ".tmp.mutation.sorted.txt", 'w')
    subprocess.check_call(["sort", "-k1,1", "-k2,2n", "-k4,4r", args.mutation_file], stdout = hout)
    hout.close()

    utils.get_multi_mutation_region(args.output_file + ".tmp.mutation.sorted.txt", args.output_file + ".tmp.multi_mutation_region.txt")

    utils.generate_template_seq(args.output_file + ".tmp.multi_mutation_region.txt", args.output_file + ".tmp.template.seq.fa", args.reference_genome)




