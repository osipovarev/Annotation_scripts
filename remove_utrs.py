#!/usr/bin/env python
#


##### NOT FUNCTIONAL SO FAR!!! #####

"""
This script takes annotation in bed12 format with UTRs and removes non-coding regions leaving only CDS
"""

import argparse
from collections import defaultdict
import sys

__author__ = "Ekaterina Osipova, 2020."


def remove_utrs(block_sizes, block_starts, utr5_size, utr3_size, exon_number):
    ## Gets a transcripts line in bed12 format with UTRs, returns bed12 line with CDS only

    block_sizes_list = map(int, [i for i in block_sizes.rstrip(',').split(',')])
    block_starts_list = map(int, [i for i in block_starts.rstrip(',').split(',')])

    if exon_number == 1:
        new_block_starts = str(block_starts_list[0]) + ','
        new_block_sizes = str(block_sizes_list[0] - utr5_size - utr3_size) + ','
    elif exon_number > 1:
        new_block_starts_list = [block_starts_list[0]] + [i - utr5_size for i in block_starts_list[1:]]
        new_block_sizes_list = [block_sizes_list[0] - utr5_size] + block_sizes_list[1:-1] + [
            block_sizes_list[-1] - utr3_size]
        new_block_starts = ','.join(map(str, new_block_starts_list)) + ','
        new_block_sizes = ','.join(map(str, new_block_sizes_list)) + ','
    else:
        sys.stderr.write('ERROR: number of exons should be at least 1! transcript has: {} exons'.format(exon_number))
        sys.exit(1)

    return new_block_sizes, new_block_starts


def process_anno_bed(bed_file):
    ## Reads annotation file and removes UTRs from transcripts

    with open(bed_file, 'r') as inf:
        for bed_line in inf.readlines():
            chrom = bed_line.split()[0]
            start = int(bed_line.split()[1])
            end = int(bed_line.split()[2])
            exon_number = int(bed_line.split()[9])
            cds_start = int(bed_line.split()[6])
            cds_end = int(bed_line.split()[7])
            block_sizes = bed_line.split()[10]
            block_starts = bed_line.split()[11]
            utr5_size = cds_start - start
            utr3_size = end - cds_end

            # get other params of annotation line
            ## chrom[0] start[1] end[2] name[3] score[4] strand[5] cds_start[6] cds_end[7] rgb[8] count[9]\
            ##  block_sizes[10] block_starts[11]
            name, score, strand, rgb = bed_line.split()[3], bed_line.split()[4], bed_line.split()[5], bed_line.split()[8]

            if (utr5_size < 0) or (utr3_size < 0):
                sys.stderr.write('ERROR: bed line is corrupted: {}'.format(bed_line))
                sys.exit(1)
            elif (utr5_size > 0) or (utr3_size > 0):
                new_block_sizes, new_block_starts = remove_utrs(block_sizes, block_starts, utr5_size, utr3_size, exon_number)
                new_exon_number = str(len(new_block_sizes.rstrip(',').split(',')))
            else:
                new_block_sizes, new_block_starts = block_sizes, block_starts
                new_exon_number = str(exon_number)

            # prepare a new annotation line
            cds_start = str(cds_start)
            cds_end = str(cds_end)
            new_bed_line = '\t'.join([chrom, cds_start, cds_end, name, score, strand, cds_start, cds_end,\
                                      rgb, new_exon_number, new_block_sizes, new_block_starts])
            print(new_bed_line)


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--anno', type=str, help='bed12 annotation file with UTRs')
    args = parser.parse_args()

    ## Read annotation file. Remove UTRs. output to stdout
    process_anno_bed(args.anno)


if __name__ == "__main__":
    main()

