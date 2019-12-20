#!/usr/bin/env python 

#--$ python GetTestTrain.py input.fasta outbase ProportonForCV

import sys
from collections import defaultdict
import random
import time

sys.stderr.write('\n\n### ' + time.ctime(time.time()) + '\n\n')

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

# set number of partitions
Parts = int(1/float(sys.argv[3]))

sys.stderr.write('\n\n### performing ' + str(Parts) + '-fold partitioning of the reference data' + '\n\n')

# subsampe fasta into test dictionary
#tmpFasta = fasta.copy()
#TestFastaLst = range(0,Parts)
for a in range(0,Parts):
	tmpTestLst = []
	with open(str(str(sys.argv[2])+str(a)+'_Test.fasta'), 'w') as tstFile:
		with open(str(str(sys.argv[2])+str(a)+'_Train.fasta'), 'w') as trnFile:
			for i in random.sample(list(fasta), SampleSize):
				fkey = '>' + i
				seq = str(fasta[i])
				#print(seq)
				tstFile.write("%s\n%s\n" % (fkey, seq[2:-2]))
				tmpTestLst.append(i)
			for key, value in fasta.items():
				if key not in tmpTestLst:
					fkey = str('>' + key)
					seq = str(value)
					trnFile.write(fkey+'\n'+seq+'\n')

sys.stderr.write('\n\n### ' + time.ctime(time.time()) + '\n\n')

# subsample fasta into train dictionary
#TrainFastaLst = range(0,Parts)
#for i in range(0,Parts):
#	TrainFastaLst[i] = {}
#	for key, value in fasta.items():
#		if key not in TestFastaLst[i].keys():
#			TrainFastaLst[i][key] = value

# write dictionaries to file
#for i in range(0,Parts):
#	with open(str(str(sys.argv[2])+str(i)+'_Test.fasta'), 'w') as File:
#		for key in TestFastaLst[i]:
#		        fkey = '>' + key
#		        seq = str(TestFastaLst[i][key])
#		        File.write("%s\n%s\n" % (fkey, seq[2:-2]))
#	with open(str(str(sys.argv[2])+str(i)+'_Train.fasta'), 'w') as File2:
#	        for key in TrainFastaLst[i]:
#	                fkey = '>' + key
#	                seq = str(TrainFastaLst[i][key])
#	                File2.write("%s\n%s\n" % (fkey, seq[2:-2]))
