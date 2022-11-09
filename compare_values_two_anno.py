#!/usr/bin/env python3
#
'''
This script takes two files with annotations of the same list of IDs (e.g genes);
and checks if values corresponding to the same ID are similar enough.
'''

import argparse
import itertools


__author__ = "Ekaterina Osipova, 2022."


def read_into_dict(file):
    ## Read genes and corresponding values into dict

    anno_dict = {}
    with open(file, 'r') as inf:
        for line in inf.readlines():
           gene = line.rstrip().split()[0]
           values = line.rstrip().split()[1]
           anno_dict[gene] = [float(v) for v in values.split(',')]
    return anno_dict


def get_pairwise_div_two_lists(l1, l2):
    ## Returns results of division of all pairs of two provided lists

    all_pairs = list(itertools.product(l1, l2))
    all_ratio = [p[0]/p[1] for p in all_pairs] + [p[1]/p[0] for p in all_pairs]
    return all_ratio


def filter_genes_by_values_diff(anno1, anno2, diff, exclude, pprint):
    ## Output genes with difference in values smaller than the threshold

    for gene in anno1:
        if gene in anno2:
            v1 = anno1[gene]
            v2 = anno2[gene]

            if pprint:
                print('{}\t{}\t{}'.format(gene, v1, v2))
            else:
                # ratio = max(v1 / v2, v2 / v1)
                all_ratio = get_pairwise_div_two_lists(v1, v2)

                if not exclude:
                    if any(ratio <= 1 + diff / 100 for ratio in all_ratio):
                        print(gene)
                else:
                    if any(ratio > 1 + diff / 100 for ratio in all_ratio):
                        print(gene)


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-a1', '--anno1', type=str, help='annotation file 1; gene name \t float value')
    parser.add_argument('-a2', '--anno2', type=str, help='annotation file 2; gene name \t float value')
    parser.add_argument('-d', '--diff', type=int, help='max difference in percent of two values for gene to stay')
    parser.add_argument('-e', '--excluded', action='store_true', help='specify if you want EXCLUDED genes instead')
    parser.add_argument('-p', '--pprint', action='store_true', help='specify if you want to just print both values side by side')
    args = parser.parse_args()

    ## Make a dictionary of genes with corresponding values from both annotations
    anno1 = read_into_dict(args.anno1)
    anno2 = read_into_dict(args.anno2)

    ## Output genes with allowed difference in values
    filter_genes_by_values_diff(anno1, anno2, args.diff, args.excluded, args.pprint)


if __name__ == "__main__":
    main()
