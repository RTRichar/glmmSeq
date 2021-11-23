#!/usr/bin/env python

# sys.argv[1]: CV_VsearchFiles, sys.argv[2]: OutfileName

import sys

# after combining all vsearch files, open combined file
# Scroll through and pick top hit for each query, add it to list
# add first Q to list
# when Q changes, add that line to list
# print list as new vsearch outfile mimic

# get initial values from line 0
#with open(sys.argv[1], 'r') as VsearchOut:	
#	Q = VsearchOut.readlines()[0].split('\t')[0]

############ Qs must be sorted and in sequenctial order (i.e. best hist first). Cannot repeat throughout file, second instance will be ignored
###### Assumes that Vsrch out has Qs ordered from best to worst hit
Qs_Found = []
TopHitList = []
with open(sys.argv[1], 'r') as VsearchOut:
	for line in VsearchOut:
		Q = line.split('\t')[0]
		if Q not in Qs_Found:
			TopHitList.append(line.strip())
			Qs_Found.append(Q)

with open(sys.argv[2], 'w') as OutFile:
	for i in TopHitList:
		OutFile.write(str(i) + '\n')
