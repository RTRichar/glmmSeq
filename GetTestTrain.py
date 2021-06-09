#!/usr/bin/env python 

#--$ python GetTestTrain.py input.fasta outbase kFolds

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
SampleSize = int(count['count'] / float(sys.argv[3]))

# set number of k-fold partitions
Parts = int(sys.argv[3])

sys.stderr.write('\n\n### performing ' + str(Parts) + '-fold partitioning of the reference data' + '\n\n')

# NOT Really 10-fold!!!!!
# currently bootstrapped, need to track test set additively so they aren't resampled
for a in range(0,Parts):
	tmpTestLst = []
	with open(str(str(sys.argv[2])+str(a)+'_Test.fasta'), 'w') as tstFile:
		with open(str(str(sys.argv[2])+str(a)+'_Train.fasta'), 'w') as trnFile:
			for i in random.sample(list(fasta), SampleSize):
				fkey = '>' + i
				seq = str(fasta[i])
				tstFile.write("%s\n%s\n" % (fkey, seq[2:-2]))
				tmpTestLst.append(i)
			for key, value in fasta.items():
				if key not in tmpTestLst:
					fkey = str('>' + key)
					seq = str(value)
					trnFile.write(fkey+'\n'+seq+'\n')

sys.stderr.write('\n\n### ' + time.ctime(time.time()) + '\n\n')

