#!/usr/bin/env python

#--$ FrmtLineages.py MtxaFrmInput 2ndTaxInput Output

import sys

# process 2ndTaxInput
COVs = {}
with open(sys.argv[2], 'r') as Covariates:
        for line in Covariates:
                GI = line.split(',')[0]
                Info = ','.join(line.strip('\n').split(',')[1:])
                COVs[GI] = Info

# process mtxa formated hits and write output
with open(sys.argv[3], 'w') as Outfile:
	with open(sys.argv[1], 'r') as Infile:
		for line in Infile:
			line = line.strip('\n').replace('\t',',')
			line = line.replace(';',',').split(',')
			GI = line[0]
			Nums = line[-2:]
			Tax = line[:-2]
			for i in range(0,9):
				if len(Tax) < 9:
					Tax.append('')
			line = str(','.join(Tax)+','+','.join(Nums)+','+COVs[GI]+'\n')
			Outfile.write(line)
