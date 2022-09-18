#!/usr/bin/env python3
#
#
#

import argparse
from collections import defaultdict
import numpy as np
import sys


__author__ = "Ekaterina Osipova, 2022."


def read_into_dict(file):
    ## Read isoformes into dict

    iso_dict = defaultdict(list)
    with open(file, 'r') as inf:
        for line in inf.readlines():
           gene = line.rstrip().split()[0]
           trans = line.rstrip().split()[1]
           iso_dict[gene].append(trans)
    return iso_dict


def read_anno_values(anno, field, col):
    ## Reads annotation file into a dictionary; consideres IDs to be in field position

    anno_dict = {}
    with open(anno, 'r') as inf:
        for line in inf.readlines():
            id = line.split()[field - 1]
            value = line.split()[col - 1]
            anno_dict[id] = value
    return anno_dict


def assign_values_to_genes(iso_dict, anno_dict):
    ## Assigns all values from each isoform to a gene

    value_gene_dict = {}
    for gene in iso_dict:
        iso_ids = iso_dict[gene]
        iso_values = []
        for i in iso_ids:
            if i in anno_dict:
                iso_values.append(anno_dict[i])
        if iso_values:
            value_gene_dict[gene] = iso_values
    return value_gene_dict


def output_gene_stats(value_gene_dict, stats):
    ## Calculates requested stats; prints

    for gene in value_gene_dict:
        all_values = [float(i) for i in list(value_gene_dict[gene])]
        if stats == 'all':
            results = ','.join([str(i) for i in all_values])
        elif stats == 'mean':
            results = mean(all_values)
        elif stats == 'std':
            results = np.var(all_values)
        elif stats == 'max':
            results = max(all_values)
        else:
            print('{} is not an appropriate value for stats!'.format(stats))
            sys.exit(1)
        print('{}\t{}'.format(gene, results))


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--anno', type=str, help='annotation file')
    parser.add_argument('-i', '--iso', type=str, help='isoformes file: gene \t transcript')
    parser.add_argument('-c', '--col', type=int, help='column with values of interest')
    parser.add_argument('-s', '--stats', type=str, default='max', help='all/mean/std/max; default=max(or longest)')
    args = parser.parse_args()

    ## Make a dictionary of transcripts
    iso_dict = read_into_dict(args.iso)

    ## Read annotation file into a dict: ID : value
    field = 1
    anno_dict = read_anno_values(args.anno, field, args.col)

    value_gene_dict = assign_values_to_genes(iso_dict, anno_dict)

    ## Output calculated stistics per gene
    output_gene_stats(value_gene_dict, args.stats)


if __name__ == "__main__":
    main()