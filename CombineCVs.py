#!/usr/bin/env python

#--$ CombineCVs.py CTEMPDIR (Parts-1) suffix

import sys

tmpLst = []
for i in range(0,int(sys.argv[2])):
	with open(str(sys.argv[1]+'/'+str(sys.argv[2])+'Fold'+str(i)+sys.argv[3])) as CVfile:
		for line in CVfile:
			tmpLst.append(line)

with open(str(sys.argv[1]+'/'+str(sys.argv[2])+'Fold'+sys.argv[3]), 'w') as OutFile: 
	for line in tmpLst:
		OutFile.write(line)
