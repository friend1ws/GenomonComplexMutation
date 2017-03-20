#! /usr/bin/env python

import sys, subprocess, os
import utils
import my_utils.pyssw

def main(args):

    hout = open(args.output_file + ".tmp.mutation.sorted.txt", 'w')
    subprocess.check_call(["sort", "-k1,1", "-k2,2n", "-k4,4r", args.mutation_file], stdout = hout)
    hout.close()

    utils.get_multi_mutation_region(args.output_file + ".tmp.mutation.sorted.txt", args.output_file + ".tmp.multi_mutation_region.txt")

    with open(args.output_file + ".tmp.multi_mutation_region.txt", 'r') as hin:
        for line in hin:

            region_chr, region_start, region_end, region_mutations = line.rstrip('\n').split('\t')

            # generate template sequences with and without mutations
            utils.generate_template_seq(region_chr, region_start, region_end, args.reference_genome, 
                                        region_mutations, args.output_file + ".tmp.template.seq.fa")

            # extract short reads around the mutations
            utils.extract_short_read(args.bam_file, args.output_file + ".tmp.short_read.seq.fa", 
                                     region_chr, region_start, region_end)

            # perform smith-waterman alignment check
            type2count = my_utils.pyssw.main2(args.output_file + ".tmp.short_read.seq.fa", args.output_file + ".tmp.template.seq.fa", 0)
   
 
            tnum = 0
            ref_count = 0
            ref_conf = ""
            is_joint = "independent"
            mut_class = "simple"
            first_count = 0
            second_count = 0
            print_conf = []
            print_count = []
            for ttype, tcount in sorted(type2count.items(), key = lambda x: x[1], reverse=True):
                treg, tconf = ttype.split(';')

                # treatment for reference alignment 
                if tconf.find("1") == -1: 
                    ref_count = tcount
                    ref_conf = tconf
                    continue
        
                # ignore template with no counts
                if tcount == 0: continue

                if tnum == 0:
                    if tconf.count("1") >= 2:
                        is_joint = "joint"
                        mut_class = utils.classify_complex_mutation(region_mutations, tconf)

                print_conf.append("(" + tconf + ")")
                print_count.append(str(tcount))

                tnum = tnum + 1

            print_conf.append("(" + ref_conf + ")")
            print_count.append(str(ref_count))

            print region_chr + '\t' + region_start + '\t' + region_end + '\t' + region_mutations + '\t' + \
                    is_joint + '\t' + mut_class + '\t' + ';'.join(print_conf) + '\t' + ';'.join(print_count)



