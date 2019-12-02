#!/usr/bin/env python

#--$ CombineCVs.py CTEMPDIR

import sys

tmpLst = []
for i in range(0,10):
	with open(str(sys.argv[1]+'/'+str(i)+'_CV.vsrch.txt')) as CVfile:
		for line in CVfile:
			tmpLst.append(line)

with open(str(sys.argv[1]+'/CV.vsrch.txt'), 'w') as OutFile:
	for line in tmpLst:
		OutFile.write(line)
