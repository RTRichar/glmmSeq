#!/usr/bin/env python

#--$ CombineCVs.py CTEMPDIR kSeries

import sys

tmpLst = []
kSeries = str(sys.argv[2]).split(',')
for i in kSeries:
	with open(str(sys.argv[1]+'/'+str(i)+'Fold_CV.LogReg.csv')) as file:
		for line in file:
			tmpLst.append(line)

with open(str(sys.argv[1]+'/FinalLogReg.csv'), 'w') as OutFile: 
	for line in tmpLst:
		OutFile.write(line)
