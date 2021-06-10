#!/usr/bin/env python

#--$ CombineCVs.py CTEMPDIR (Parts-1) suffix

import sys

tmpLst = []
for i in range(0,int(sys.argv[2])):
	with open(str(sys.argv[1]+'/'+str(i)+'_'+sys.argv[3])) as CVfile:
		for line in CVfile:
			tmpLst.append(line)

with open(str(sys.argv[1]+'/'+sys.argv[3]), 'w') as OutFile:  # '/CV.vsrch.txt'
	for line in tmpLst:
		OutFile.write(line)
