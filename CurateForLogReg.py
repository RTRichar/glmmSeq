#!/usr/bin/env python

# sys.argv[1]: input predicted taxonomies, sys.argv[2]: actual taxonomies, sys.argv[3]: Covariates table , sys.argv[4]: outfile name, argv[5]: k

import sys

TAX = open(sys.argv[2], 'r')
CLASS = open(sys.argv[1], 'r')
OUT = open(sys.argv[4], 'w')

# make dictionary of actual GIs and taxonomies
ACT_TAX = {}
for line in TAX:
	line = line.strip().split('\t')
	GI = line[0]
	LINEAGE = line[1].split(';')
	for i in range(0,8):
		if len(LINEAGE) < 8:
			LINEAGE.append('')
	ACT_TAX[GI] = ','.join(LINEAGE)
# make dictionary of GIs and predicted Taxonomies
PRDCT_TAX = {}
for line in CLASS:
        line = line.strip().split('\t')
        GI = line[0]
        LINEAGE = line[1].split(';')
	for i in range(0,8):
                if len(LINEAGE) < 8:
                        LINEAGE.append('')
	SCORE = line[-1]
	PI = line[-2]
	X = str(','.join(LINEAGE) + ',' + SCORE + ',' + PI)
        PRDCT_TAX[GI] = X

# make dict of GIs and distance to second tax hit and number of top tax hits
COVs = {}
with open(sys.argv[3], 'r') as Covariates:
	for line in Covariates:
		GI = line.split(',')[0]
		Info = ','.join(line.strip('\n').split(',')[1:])
		COVs[GI] = Info

for key in PRDCT_TAX:
	OUT.write(str(key) + ',' + ACT_TAX[key] + PRDCT_TAX[key] + ',' + COVs[key] + ',' + str(sys.argv[5]) + '\n')

OUT.close()
TAX.close()
CLASS.close()
