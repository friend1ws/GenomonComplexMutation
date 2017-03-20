#! /usr/bin/env python

import sys
import pysam
import my_utils.seq

def get_multi_mutation_region(input_file, output_file, dist_thres = 20):

    temp_chr = ""
    temp_start = 0
    temp_end = 0
    temp_mut = []

    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[0].startswith("Chr") or F[0].startswith("#"): continue

            if F[0] != temp_chr or int(F[1]) - temp_end > dist_thres:

                if len(temp_mut) >= 2:

                    # filter SNVs within deletions
                    del_regions = []
                    for tmut in temp_mut:
                        FF = tmut.split(',')
                        if FF[4] == "-": del_regions.append([int(FF[1]), int(FF[2])])

                    temp_mut2 = []
                    for tmut in temp_mut:
                        FF = tmut.split(',')
                        del_flag = 0
                        for del_region in del_regions:
                            if int(FF[1]) >= del_region[0] and int(FF[2]) <= del_region[1]: del_flag = 1
                        
                        if del_flag == 0:
                            temp_mut2.append(tmut)
                        
                    if len(temp_mut2) >= 2:
                        print >> hout, temp_chr + '\t' + str(temp_start) + '\t' + str(temp_end) + '\t' + ';'.join(temp_mut2)
  
                temp_chr = F[0]
                temp_start = int(F[1])
                temp_end = int(F[2])
                temp_mut = [','.join(F[0:5])]

            else:
                temp_mut.append(','.join(F[0:5]))
                temp_end = int(F[2])


    hout.close()



def generate_configurations(dim):

    conf = [[0], [1]]

    for k in range(dim - 1):

        new_conf = []
        for elm in conf:
            new_conf.append(elm + [0])
            new_conf.append(elm + [1])
            conf = new_conf

    return conf


def generate_template_seq(region_chr, region_start, region_end, reference, region_mutations, output_file, template_margin = 200):

    hout = open(output_file, 'w')

    tmp_region = region_chr + ':' + str(int(region_start) - template_margin + 1) + '-' + str(int(region_end) + template_margin - 1)
    org_seq = my_utils.seq.get_seq(reference, tmp_region)

    mutations_in_region = region_mutations.split(';')
    for conf in generate_configurations(len(mutations_in_region)):

        mut_seq = org_seq
        off_set = 0
        for i in range(len(mutations_in_region)):

            if conf[i] == 0: continue

            mut_chr, mut_start, mut_end, mut_ref, mut_alt = mutations_in_region[i].split(',')

            # in this case, we remove nucleotides from concatenated sequences
            mut_start_rel = int(mut_start) - int(region_start) + template_margin - 1 + off_set
            mut_end_rel = int(mut_end) - int(region_start) + template_margin - 1 + off_set

            # SNV
            if mut_ref != '-' and mut_alt != '-': 

                # for debug
                if mut_seq[mut_start_rel] != mut_ref:
                    print >> sys.stderr, '\t'.join([mut_chr, str(mut_start), str(mut_end), mut_ref, mut_alt])
                    print >> sys.stderr, "mutation inconsistent!!!"
                    sys.exit(1)
                mut_seq = mut_seq[:mut_start_rel] + mut_alt + mut_seq[(mut_end_rel + 1):]


            elif mut_alt == '-': # deletion

                # for debug
                if mut_seq[mut_start_rel:(mut_end_rel + 1)] != mut_ref != '-':
                    print >> sys.stderr, '\t'.join([mut_chr, str(mut_start), str(mut_end), mut_ref, mut_alt])
                    print >> sys.stderr, "mutation inconsistent!!!"
                    sys.exit(1)
                mut_seq = mut_seq[:mut_start_rel] + mut_seq[(mut_end_rel + 1):]
                off_set = off_set - len(mut_ref)

            elif mut_ref == '-': #insertion
    
                mut_seq = mut_seq[:(mut_start_rel + 1)] + mut_alt + mut_seq[(mut_start_rel + 1):]
                off_set = off_set + len(mut_alt)

        print >> hout, '>' + region_chr + ':' + region_start + '-' + region_end + ';' + ','.join([str(x) for x in conf])
        print >> hout, mut_seq


    hout.close()



def extract_short_read(bam_file, output_file, region_chr, region_start, region_end, short_read_margin = 5):

    bamfile = pysam.Samfile(bam_file, 'rb')
    hout = open(output_file, 'w') 
    for read in bamfile.fetch(region_chr, int(region_start), int(region_end)):

        # get the flag information
        flags = format(int(read.flag), "#014b")[:1:-1]

        # skip unmapped read 
        if flags[2] == "1" or flags[3] == "1": continue 

        # skip supplementary alignment
        if flags[8] == "1" or flags[11] == "1": continue

        # skip duplicated reads
        if flags[10] == "1": continue

        reference_start_pos = read.reference_start + 1
        reference_end_pos = read.reference_end

        # remove short reads that does not cover the entire mutation regions 
        if reference_start_pos > int(region_start) - short_read_margin: continue
        if reference_end_pos < int(region_end) + short_read_margin: continue

        read_id = read.qname + '/1' if flags[6] == '1' else read.qname + '/2'
        
        
        print >> hout, '>' + read_id + '\n' + read.seq

    bamfile.close()
    hout.close()



def classify_complex_mutation(region_mutations, configuration):

    mutations_in_region = region_mutations.split(';')
    active_vector = configuration.split(',')

    positions = []
    is_cindel = 0
    for i in range(len(mutations_in_region)):
        if active_vector[i] == "0": continue
        mut_chr, mut_start, mut_end, mut_ref, mut_alt = mutations_in_region[i].split(',')
        if mut_ref == "-" or mut_alt == "-": is_cindel = 1
        positions.append(int(mut_start))

    positions = sorted(positions)
    diff_positions = []
    is_continuous = 1
    for i in range(len(positions) - 1):
        if positions[i + 1] - positions[i] != 1: is_continuous = 0 


    if is_cindel == 1:
        return "CINDEL"
    elif is_continuous == 1:
        if len(positions) == 2:
            return "DNV"
        elif len(positions) == 3:
            return "TNV"
        else:
            return "MNV"
    else:
        return "INV"

            


