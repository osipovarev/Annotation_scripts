#!/usr/bin/env python3
#

import argparse
import pyfastx
import collections


__author__ = "Ekaterina Osipova, 2020."


def read_rename_dict(file):
    ## Reads correspondence table into dictionary

    rename_table = {}
    with open(file, 'r') as inf:
        for line in inf.readlines():
            rename_table[line.split(',')[0]] = line.rstrip().split(',')[1]
    return rename_table


def read_fasta(fasta_file):
    ## Reads fasta into a dictionary

    fasta_dict = {}
    for name, seq in pyfastx.Fasta(fasta_file, uppercase=False,  build_index=False):
        fasta_dict[name] = seq
    return collections.OrderedDict(sorted(fasta_dict.items()))


def rename_headers_from_dict(fasta_dict, rename_table):
    ## Replaces headers of fasta with names from rename_table

    for name in fasta_dict:
        seq = fasta_dict[name]
        print('>{}'.format(rename_table[name]))
        print(seq)


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', type=str, help='fasta file')
    parser.add_argument('-d', '--renamedict', type=str, help='table of name correspondence, usually renaming_dictionary.csv')
    args = parser.parse_args()

    ## Read fasta file into a dictionary
    fasta_dict = read_fasta(args.fasta)

    ## Read renaming dictionary
    rename_table = read_rename_dict(args.renamedict)

    ## Output fasta with headers renamed
    rename_headers_from_dict(fasta_dict, rename_table)


if __name__ == "__main__":
    main()