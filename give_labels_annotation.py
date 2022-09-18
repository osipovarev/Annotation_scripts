#!/usr/bin/env python3

import argparse
import sys


__author__ = "Ekaterina Osipova, 2021."



def read_anno_into_dict(anno, field):
    ## Reads bed12 annotation file into a dictionary

    anno_dict = {}
    with open(anno, 'r') as inf:
        for line in inf.readlines():
            if len(line.split()) < field - 1:
                print('There is no field {} in the annotation! Abort'.format(field))
                sys.exit(1)
            else:
                id = line.split()[field - 1]
                transc_info = line.rstrip()
                anno_dict[id] = transc_info
    return anno_dict


def give_labels(anno_dict, list_to_label, label):
    ##

    for t in anno_dict:
        if t in list_to_label:
            new_t = t + label
            print(anno_dict[t].replace(t, new_t))
        else:
            print(anno_dict[t])
            pass


def read_list_from_file(file):
    ##

    with open(file, 'r') as inf:
        list_from_file = [t.rstrip() for t in inf.readlines()]
    return list_from_file


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--trans_list', type=str, help='list of transcripts to label')
    parser.add_argument('-a', '--anno', type=str, help='annotation file')
    parser.add_argument('-l', '--label', type=str, help='label to give; e.g: _potentialNMD')
    parser.add_argument('-f', '--field', type=int, default=4, help='field number to label; default=4(transciprt_ID in bed12)')
    args = parser.parse_args()

    ## bed12 format:
    ## chrom[0] start[1] end[2] name[3] score[4] strand[5] cds_start[6] cds_end[7] rgb[8] count[9]\
    ##  block_sizes[10] block_starts[11]

    ## Read annotation into a dictionary: {transc_id : transc_info}
    anno_dict = read_anno_into_dict(args.anno, args.field)

    ## Read file with labels
    list_to_label = read_list_from_file(args.trans_list)

    ## Give labels to the transcripts in list_to_label
    give_labels(anno_dict, list_to_label, args.label)


if __name__ == '__main__':
    main()
