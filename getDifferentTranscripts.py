#!/usr/bin/env python
#
#
"""
This script takes two files in bed12 format and outputs transcripts that are present in only one file;
ignores color and score fields.
"""

import argparse
import sys

__author__ = "Ekaterina Osipova, 2019."


parser = argparse.ArgumentParser()
parser.add_argument('-f1', '--file1', type=str, help='file1 in bed12 format')
parser.add_argument('-f2', '--file2', type=str, help='file2 in bed12 format')
args = parser.parse_args()

transcDict1 = {}
transcDict2 = {}
fileCount = 0

for file in [args.file1, args.file2]:
    fileCount += 1
    with open(file, 'r') as inf:
        for line in inf.readlines():
            transcrName, score, color = line.split()[3], line.split()[4], line.split()[8]
            transcrInfo = '\t'.join(line.split()[:3] + line.split()[5:8] + line.split()[9:])
            if fileCount == 1:
                transcDict1[transcrInfo] = (transcrName, score, color)
            else:
                transcDict2[transcrInfo] = (transcrName, score, color)

#transcNotInDict1 = { k: transcDict2[k] for k in set(transcDict2) - set(transcDict1)}
#transcNotInDict2 = { k: transcDict1[k] for k in set(transcDict1) - set(transcDict2)}

dictCount = 0
for dict in [transcDict1, transcDict2]:
    dictCount += 1
    if dictCount == 1:
        print('NOT in the first file: ')
        otherDict = transcDict2
    else:
        print('NOT in the second file: ')
        otherDict = transcDict1
    for transcrInfo in dict:
        if transcrInfo not in otherDict:
            transcrNames = dict[transcrInfo]
            newBedLine = '\t'.join(transcrInfo.split()[:3]) + '\t' + transcrNames[0] + '\t' + \
                         transcrNames[1] + '\t' + '\t'.join(transcrInfo.split()[3:6]) + '\t' + \
                         transcrNames[2] + '\t' + '\t'.join(transcrInfo.split()[6:])
            print(newBedLine)

