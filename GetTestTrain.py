#!/usr/bin/env python 

#--$ python GetTestTrain.py input.fasta outbase k-folds  # (used to be ProportonForCV)

import sys
import random
import time

sys.stderr.write('\n\n### ' + time.ctime(time.time()) + '\n\n')

def diff(list1, list2):
	return list(set(list1).symmetric_difference(set(list2)))

fasta = {}
FastStatus = str('\n') # printed after fasta is parsed, tells user if fasta had duplicates
with open(sys.argv[1], 'r') as Fasta:
	for line in Fasta:
		if not line.strip():
			continue
		if line.startswith('>'):
			if str(line.strip())[1:] not in fasta:
				header = str(line.strip())[1:]
				continue
			else:
				FastStatus = str('\n# WARNING: duplicate sequence headers skipped')
		fasta[header] = line.strip()
sys.stderr.write(FastStatus)

# set number of partitions
Parts = int(sys.argv[3])
print('k-folds during test/train partitioning: ' + str(Parts))
Prop = float(1)/int(sys.argv[3])

# set number of seqs to retrieve
count = {'count': 0}
for key in fasta:
	count['count'] += 1
SampleSize = int(count['count'] * Prop)
print('Test sequence case sample size: ' + str(SampleSize))

# set number of partitions
#Parts = int(1/float(sys.argv[3]))

sys.stderr.write('\n\n### performing ' + str(Parts) + '-fold partitioning of the reference data' + '\n\n')

# make dict of lists, one list for each part
UsedSeqsLst = []
SeqsNotUsed = fasta.keys()
PartLstsDct = {}
for i in range(0,Parts):
	PartLstsDct[i] = random.sample(SeqsNotUsed, SampleSize)
	for seq in PartLstsDct[i]:
		UsedSeqsLst.append(seq)
	SeqsNotUsed = []
	for seq in diff(fasta.keys(),UsedSeqsLst):
		SeqsNotUsed.append(seq)

# for each part, write file
for a in range(0,Parts):
	with open(str(str(sys.argv[2])+str(a)+'_Test.fasta'), 'w') as tstFile:
		with open(str(str(sys.argv[2])+str(a)+'_Train.fasta'), 'w') as trnFile:
			for i in fasta:
				fkey = '>' + i
				seq = str(fasta[i])
				if i in PartLstsDct[a]:
					tstFile.write("%s\n%s\n" % (fkey, seq[2:-2]))
				else:
					trnFile.write("%s\n%s\n" % (fkey, seq[2:-2]))

sys.stderr.write('\n\n### ' + time.ctime(time.time()) + '\n\n')
