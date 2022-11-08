#!/usr/bin/env python3

'''
This script takes an annotation gff3 file
and outputs table: transcript\t geneID\t geneName
'''

import argparse


__author__ = "Ekaterina Osipova, 2022."


def extract_ids_from_exon_lines(file):

	transc_dict = {}
	with open(file, 'r')as inf:
		for line in inf.readlines():
			line_elements = line.split('\t')

			## check if it a gff annotation line
			if len(line_elements) == 9:
				if line_elements[2] == 'exon':
					exon_info = line_elements[8]
					exon_elements = exon_info.replace(',', ';').split(';')
					for el in exon_elements:
						if el.startswith('Parent='):
							transc = el.replace('Parent=', '')
						elif el.startswith('Dbxref=GeneID:'):
							ncbigene = el.replace('Dbxref=GeneID:', '')
						elif el.startswith('gene='):
							genename = el.replace('gene=', '')
					transc_dict[transc] = (ncbigene, genename)
	return transc_dict


def output_transcripts(transc_dict):

	for transc in transc_dict:
		ncbigene = transc_dict[transc][0]
		genename = transc_dict[transc][1]
		print('{}\t{}\t{}'.format(transc, ncbigene, genename))


def main():
	## Parse arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', '--annogff', type=str, help='annotation file in gff format')
	args = parser.parse_args()

	## Parse gff file into a transcript dictionary
	transc_dict = extract_ids_from_exon_lines(args.annogff)

	## Output the requested table
	output_transcripts(transc_dict)


if __name__ == '__main__':
	main()