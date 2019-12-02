#!/usr/bin/env python

#--$ FrmtLineages.py input output

import sys

with open(sys.argv[2], 'w') as Outfile:
	with open(sys.argv[1], 'r') as Infile:
		for line in Infile:
			line = line.replace('\t',',')
			line = line.replace(';',',').split(',')
			Nums = line[-2:]
			Tax = line[:-2]
	#		print('#######NEW LINE')
	#		print(Nums)
	#		print(Tax)
	#		print(str(len(Tax)))
			for i in range(0,9):
				if len(Tax) < 9:
					Tax.append('')
	#		print(Tax)
	#		print(str(len(Tax)))
			line = str(','.join(Tax)+','+','.join(Nums))
			Outfile.write(line)
