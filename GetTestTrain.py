#!/usr/bin/env python 

#--$ python GetTestTrain.py input.fasta outbase ProportonForCV

import sys
from collections import defaultdict
import random

# populate dictionary with fasta
fasta = defaultdict(list)
with open(sys.argv[1]) as file_one:
    for line in file_one:
        if line.startswith(">"):
            fasta[line.strip(">\n")].append(next(file_one).rstrip())

# set number of seqs to retrieve
count = {'count': 0}
for key in fasta:
	count['count'] += 1
SampleSize = int(count['count'] * float(sys.argv[3]))

# subsampe fasta into test dictionary
tmpFasta = fasta.copy()
TestFastaLst = range(0,10)
for a in range(0,10):
	TestFastaLst[a] = {}
	for i in random.sample(list(tmpFasta), SampleSize):
		TestFastaLst[a][i] = tmpFasta[i]
		del tmpFasta[i]

# subsample fasta into train dictionary
TrainFastaLst = range(0,10)
for i in range(0,10):
	TrainFastaLst[i] = {}
	for key, value in fasta.items():
		if key not in TestFastaLst[i].keys():
			TrainFastaLst[i][key] = value

# write dictionaries to file
for i in range(0,10):
	with open(str(str(sys.argv[2])+str(i)+'_Test.fasta'), 'w') as File:
		for key in TestFastaLst[i]:
		        fkey = '>' + key
		        seq = str(TestFastaLst[i][key])
		        File.write("%s\n%s\n" % (fkey, seq[2:-2]))
	with open(str(str(sys.argv[2])+str(i)+'_Train.fasta'), 'w') as File2:
	        for key in TrainFastaLst[i]:
	                fkey = '>' + key
	                seq = str(TrainFastaLst[i][key])
	                File2.write("%s\n%s\n" % (fkey, seq[2:-2]))
