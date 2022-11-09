#!/usr/bin/env python3

'''
This script takes an isoforms file: ABC1\tABC1_rna-XMxxxx
and outputs table: ABC1\trna-XMxxxx
'''

import argparse


__author__ = "Ekaterina Osipova, 2022."


## Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--isoforms', type=str, help='isoforms file: ABC1\tABC1_rna-XMxxxx')
args = parser.parse_args()

with open(args.isoforms, 'r') as inf:
	for line in inf.readlines():
		g = line.split()[0]
		t = line.split()[1]
		new_t = g.replace(g + '_', '')
		new_line = '{}\t{}'.format(g, new_t)
		print(new_line)