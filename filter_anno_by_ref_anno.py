#!/usr/bin/env python3


"""
This script checks every transcript of a given annotation if it has a close
enough (by size) equivalent in the reference annotation
"""

import argparse
import subprocess
import numpy as np
import tempfile
import sys


__author__ = "Ekaterina Osipova, 2021."


def read_anno_into_dict(anno):
    ## Reads bed12 annotation file into a dictionary

    anno_dict = {}
    with open(anno, 'r') as inf:
        for line in inf.readlines():
            id = line.split()[3]
            transc_info = line.rstrip()
            anno_dict[id] = transc_info
    return anno_dict


# def get_overlapping_transcripts(bed12_line, bed12_file, fraction):
#     ## Used bedtools to get all transcripts from bed12_file overlapping at least $fraction of bed12_line
#
#
#     # write bed12_line to a named tmp file
#     tmp = tempfile.NamedTemporaryFile('w+t')
#     tmp.write(bed12_line + '\n')
#     tmp.seek(0)
#
#     # bedtools intersect -split -a $anno_bed -b <(echo "$toga_loci") -wo -f $fraction
#     bedtools_cmd = 'bedtools intersect -split -a {} -b {} -wo -f {}'.format(bed12_file, tmp.name, fraction)
#     run_cmd = subprocess.Popen(bedtools_cmd, shell=True, stdout=subprocess.PIPE)
#     stdout_cmd, err = run_cmd.communicate()
#     transc_list = stdout_cmd.decode().rstrip().split('\n')
#
#     # remove temp file
#     tmp.close()
#     return transc_list


def run_bedtools_intersect(file_a, file_b, fraction):
    ## Runs bedtools intersect for two provided files to get pairs of overlapping transcripts

    bedtools_cmd = 'bedtools intersect -split -a {} -b {} -wo -f {}'.format(file_a, file_b, fraction)
    run_cmd = subprocess.Popen(bedtools_cmd, shell=True, stdout=subprocess.PIPE)
    stdout_cmd, err = run_cmd.communicate()
    a_b_transcript_pairs = stdout_cmd.decode().rstrip().split('\n')
    return a_b_transcript_pairs


def get_bases(starts, sizes, abs_start):
    ## Converts provided intervals into lists of bases

    list_2d = [list(range(abs_start + starts[i], abs_start + starts[i] + sizes[i])) for i in range(len(starts))]
    return list(np.concatenate(list_2d))


def check_a_b_overlap(transc_list, minr, maxr):
    ## Checks a list of pairs of transcripts from bedtools intersect

    good_transcripts = []
    dropped_transcripts = []
    for transc_pair in transc_list:
        a_bed = transc_pair.split()[:12]
        b_bed = transc_pair.split()[12:25]
        trans_name = b_bed[3]

        a_abs_start = int(a_bed[6])
        b_abs_start = int(b_bed[6])
        a_starts = [int(i) for i in a_bed[11].rstrip(',').split(',')]
        b_starts = [int(i) for i in b_bed[11].rstrip(',').split(',')]
        a_block_sizes = [int(i) for i in a_bed[10].rstrip(',').split(',')]
        b_block_sizes = [int(i) for i in b_bed[10].rstrip(',').split(',')]

        a_bases = get_bases(a_starts, a_block_sizes, a_abs_start)
        b_bases = get_bases(b_starts, b_block_sizes, b_abs_start)
        bases_overlap = len(list(set(a_bases) & set(b_bases)))

        big_delta = float(len(b_bases)) / float(bases_overlap)
        small_delta = float(bases_overlap) / float(len(a_bases))
        # print('a_bases min, max = {}, {}'.format(min(a_bases), max(a_bases)))
        # print('b_bases min, max = {}, {}'.format(min(b_bases), max(b_bases)))
        # print('bases_overlap = {}'.format(bases_overlap))
        # print('big_delta = {}'.format(big_delta))
        # print('small_delta = {}'.format(small_delta))

        if (big_delta <= maxr) and (small_delta >= minr):
            good_transcripts.append(trans_name)
        else:
            dropped_transcripts.append(trans_name)
    return good_transcripts, dropped_transcripts


def output_transcripts(transc_list, anno_dict, stdout=True):
    ## Outputs full bed12 lines of transcripts

    for trans in set(transc_list):
        if stdout:
            print(anno_dict[trans])
        else:
            # sys.stderr(anno_dict[trans])
            print(anno_dict[trans], file=sys.stderr)


def main():
    ## Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--refanno', type=str, help='reference annotation')
    parser.add_argument('-a', '--anno', type=str, help='annotation to check with reference')
    parser.add_argument('-minr', '--minratio', type=float, default=0.3, help='min length of a transcript in fraction-ref to keep; default=0.3')
    parser.add_argument('-maxr', '--maxratio', type=float, default=1.5, help='min length of a transcript in fraction-ref to keep; default=1.5')
    parser.add_argument('-d', '--drop', action='store_true', help='output dropped transcripts into stderr')
    args = parser.parse_args()

    ## bed12 format:
    ## chrom[0] start[1] end[2] name[3] score[4] strand[5] cds_start[6] cds_end[7] rgb[8] count[9]\
    ##  block_sizes[10] block_starts[11]

    ## Read annotation into a dictionary: {transc_id : transc_info}
    query_anno_dict = read_anno_into_dict(args.anno)

    ## Run bedtools intersect
    a_b_transcript_pairs = run_bedtools_intersect(args.refanno, args.anno, 0.5)

    ## Filter transcripts based on a_b pairs overlap
    good_transcripts, dropped_transcripts = check_a_b_overlap(a_b_transcript_pairs, args.minratio, args.maxratio)
    all_transcripts = [k for k in query_anno_dict]

    ## Output filtered annotation
    good_and_nooverlap_transcripts = list(set(all_transcripts) - set(dropped_transcripts))
    output_transcripts(good_and_nooverlap_transcripts, query_anno_dict)

    ## If requested, output dropped transcripts to stderr
    if args.drop:
        output_transcripts(dropped_transcripts, query_anno_dict, stdout=False)



if __name__ == "__main__":
    main()