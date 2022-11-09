#!/usr/bin/env python3
#
import argparse
from collections import defaultdict


__author__ = "Ekaterina Osipova, 2022."


def read_into_dict(file):
    ## Read isoformes into dict

    iso_dict = defaultdict(list)
    with open(file, 'r') as inf:
        for line in inf.readlines():
            elements = line.rstrip().split()
            if len(elements) == 2:
                gene = line.rstrip().split()[0]
                trans = line.rstrip().split()[1]
                iso_dict[gene].append(trans)
    return iso_dict


def reverse_dict(d):
    ## Flips keys and values is a dict; makes set of new keys

    rev_d = {}
    for k, v in d.items():
        new_k = tuple(set(v))
        new_v = k
        rev_d[new_k] = new_v
    return rev_d


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i1', '--iso1', type=str, help='isoformes file 1: gene \t transcript')
    parser.add_argument('-i2', '--iso2', type=str, help='isoformes file 2: gene \t transcript')
    args = parser.parse_args()

    ## Read both isoform files inot dictionaries
    iso_dict1 = read_into_dict(args.iso1)
    iso_dict2 = read_into_dict(args.iso2)

    ## Reverse both dictionaries so they represent set(transcripts) : gene
    rev_iso_dict1 = reverse_dict(iso_dict1) 
    rev_iso_dict2 = reverse_dict(iso_dict2)

    ## Try to find for each entry in isofomrs 1 an equivalent from isoforms 2
    for t in rev_iso_dict1:
        if t in rev_iso_dict2:
            print("{} == {}".format(rev_iso_dict1[t], rev_iso_dict2[t]))
            pass
        else:
            #print("Found no equivalent for {}. {}. {}".format(rev_iso_dict1[t], t, iso_dict2[rev_iso_dict1[t]]))
            pass


if __name__ == "__main__":
    main()
