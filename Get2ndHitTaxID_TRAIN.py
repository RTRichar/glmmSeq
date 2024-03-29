#!/usr/bin/env python

# sys.argv[1]: Taxonomies, sys.argv[2]: CV_2nd_VsearchFiles, sys.argv[3]: OutfileName

import sys

# make dictionary of tax info 
TaxDct = {}
with open(sys.argv[1], 'r') as TaxFile:
	for line in TaxFile:
		TaxDct[line.split('\t')[0]] = line.strip().split('\t')[1].split(';')

# run through 2nd tax alignment file once for each taxonomic rank, record second tax percent ID
# get initial values from line 0
ScndTaxIDs = {} # ScndTaxIDs[Q] = [Ord2ndID,Fam2ndID,Gen2ndID,Sp2ndID]
FrstTaxQs = {} # number of hits to first tax
with open(sys.argv[2], 'r') as VsearchOut:	
	L = VsearchOut.readlines()[0].split('\t')[0:2]
	Q = L[0]
	T = L[1]
ScndIDstat = 'notfound'
Ql = 0

for i in [3,4,5,6]:
	with open(sys.argv[2], 'r') as VsearchOut:
		for line in VsearchOut:
			if line.split('\t')[0] == Q: # if Q same			
				if ScndIDstat == 'notfound': # and Second ID not found
					if TaxDct[line.split('\t')[0]][i] != '': # Q taxa not blank # Q Taxa only during training, not classification
						if TaxDct[line.split('\t')[1]][i] != '': # T not blank
							Ql += 1
							if TaxDct[line.split('\t')[1]][i] != TaxDct[T][i]: # and T changed
								ScndID = line.split('\t')[2]
								ScndIDstat = 'found'
								if i == 3:
									ScndTaxIDs[Q] = [str(ScndID)]
									FrstTaxQs[Q] = [str(Ql)]
								else:
									ScndTaxIDs[Q].append(str(ScndID))
									FrstTaxQs[Q].append(str(Ql))
			else: # if Q changed
				if ScndIDstat == 'notfound': # but alnmt ID of second tax not yet found
					if i == 3:
						ScndTaxIDs[Q] = ['']
						FrstTaxQs[Q] = [str(Ql)]
					elif Q not in ScndTaxIDs.keys():# somehow we aren't on 3 but Q hasn't been initiated in Dct {prob fixed}
						ScndTaxIDs[Q] = [''] # need to double check {prob fixed}
						FrstTaxQs[Q] = [str(Ql)]
					else:
						ScndTaxIDs[Q].append('')
						FrstTaxQs[Q].append(str(Ql))
				Q = line.split('\t')[0] 
				T = line.split('\t')[1]
				ScndID = ''
				ScndIDstat = 'notfound'
				Ql = 0

with open(sys.argv[3], 'w') as OutFile:
	OutFile.write(','.join(['Query','o2ndID','f2ndID','g2ndID','s2ndID','oTpTxLngth','fTpTxLngth','gTpTxLngth','sTpTxLngth','\n']))
	for key in ScndTaxIDs:
		OutFile.write(str(key + ',' + ','.join(ScndTaxIDs[key]) + ',' + ','.join(FrstTaxQs[key]) + '\n'))
