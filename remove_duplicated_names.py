#!/usr/bin/env python3
#
#
# Identifies name duplications in bed annotation; output the LONGER one
# e.g usage: remove_duplicated_names.py -f HLparMaj1.ncbi.bed > uniq.HLparMaj1.ncbi.bed
#

import argparse
from collections import defaultdict
from operator import itemgetter


__author__ = "Ekaterina Osipova, 2021."


def read_annotation(file):
    ## Reads bed12 annotation file into a dict

    anno_dict = defaultdict(list)
    with open(file, 'r') as inf:
        for line in inf.readlines():
            trans_name = line.split()[3]
            exons = line.split()[10].rstrip(',').split(',')
            exons_cov = sum([int(i) for i in exons])
            anno_dict[trans_name].append((exons_cov, line.rstrip()))
    return anno_dict


def output_uniq_transcripts(anno_dict):
    ## Outputs unique tramscripts; id case of duplications, outputs the longer one

    for trans in anno_dict:
        if len(anno_dict[trans]) == 1:
            trans_info = anno_dict[trans][0][1]
            print(trans_info)
        else:
            trans_list = anno_dict[trans]
            longest_trans = max(trans_list, key=itemgetter(0))[1]
            print(longest_trans)


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--filebed', type=str, help='bed12 file')
    args = parser.parse_args()

    ## Make a dictionary of transcripts
    anno_dict = read_annotation(args.filebed)

    ## Print unique transcripts
    output_uniq_transcripts(anno_dict)


if __name__ == "__main__":
    main()