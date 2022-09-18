#!/usr/bin/env python3
#

import argparse
from collections import defaultdict


__author__ = "Ekaterina Osipova, 2022."


def read_anno_to_dict(file, column):
	## Reads annotation file into a dict; checks IDs for uniqueness

	anno_dict = defaultdict(list)
	with open(file , 'r') as inf:
		for line in inf.readlines():
			info = line.rstrip()
			id_curr = info.split()[column - 1]
			anno_dict[id_curr].append(info)

	return anno_dict


def enumerate_nonuniq(anno_dict, suffix):
	## Goes through dictionary, if a name corresponds \
	## to multiple transcripts, gives them uniq lables

	for id_curr in anno_dict:
		all_transcripts = anno_dict[id_curr]
		
		if len(all_transcripts) == 1:
			print(all_transcripts[0])
		else:
			for i in range(len(all_transcripts)):
				new_id = '{}{}{}'.format(id_curr, suffix, i)
				new_info = all_transcripts[i].replace(id_curr, new_id)
				print(new_info)


def main():
	## Parse arguments
	parser = argparse.ArgumentParser()
	
	parser.add_argument(
		'-a',
		'--annotation',
		type=str,
		help='annotation file'
		)
	parser.add_argument(
		'-c',
		'--column',
		type=int,
		default=4,
		help='column number with IDs; default=4 (bed12)'
		)
	parser.add_argument(
		'-s',
		'--suffix', 
		type=str,
		default='_', 
		help='suffix to give in addition to enumeration')
	
	args = parser.parse_args()

	## Read annotation file into a dictionary
	anno_dict = read_anno_to_dict(args.annotation, args.column)

	## Give unique label to duplicated IDs and output updated annotation
	enumerate_nonuniq(anno_dict, args.suffix)



if __name__ == "__main__":
	main()
