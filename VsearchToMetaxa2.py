#!/usr/bin/env python

import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--VsearchOutput', required = True, help = "\n ttt\n")
parser.add_argument('-t', '--Taxonomies', required = True, help = "\n ttt\n")
parser.add_argument('-o', '--OutPutFile', required = True, help = "\n ttt\n")
args = parser.parse_args()

TaxDct = {}
with open(args.Taxonomies, 'r') as TaxFile:
	for line in TaxFile:
		TaxDct[line.split('\t')[0]] = line.strip().split('\t')[1] 

vDct = {}
ID = {}
Lngth = {}
with open(args.VsearchOutput, 'r') as Vsearch:
	for line in Vsearch:
		vDct[line.split('\t')[0]] = line.split('\t')[1]
		ID[line.split('\t')[0]] = line.split('\t')[2]
		Lngth[line.split('\t')[0]] = line.split('\t')[3]

with open(args.OutPutFile, 'w') as OutFile:
	for key in vDct:
		OutFile.write(str(key) + '\t' + str(TaxDct[vDct[key]]) + '\t' + str(ID[key]) + '\t' + str(Lngth[key]) + '\n')
