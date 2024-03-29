#!/usr/bin/env python

import subprocess, argparse, sys, time, os

parser = argparse.ArgumentParser(description="")
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
# required
required.add_argument('-i', '--InputFasta', required = True, help = "\n Input file of reference sequences to be extracted (sequence headers should include a unique sequence identifier, NCBI accessions are recomended, preceded only by '>')\n")
required.add_argument('-db', '--Database', required = True, help = "\n Name of database to be used for classification. Database must be ... \n")
required.add_argument('-o', '--Output', required = True, help = "\n Na \n")
#optional
optional.add_argument('-t', '--Threads', required = False, default = 1, help = "\n Number of processors for Vsearch alignment\n")
optional.add_argument('-id', '--idCutoffs', required = False, default = 1, help = "\n Option for")
optional.add_argument('--SaveTemp', default=False, type=lambda x: (str(x).lower() == 'true'), help = "\nOption for \n")
args = parser.parse_args()

### create temp directory
CTEMPDIR = str('Classify_tmp_' + '_'.join(time.ctime(time.time()).replace(':','_').split()[1:]))
subprocess.call(['mkdir', CTEMPDIR])

# retreive db directory
DBDIR = str(os.path.abspath(os.path.dirname(sys.argv[0]))+'/DBs/')

# Set RE level and highest rank to be analyzed
with open(str(DBDIR+args.Database+'/'+'InfoFile.txt'), 'r') as File:
	line = File.readlines()[0].split(';')
	print(line)
	reStructure = line[0] 
	HighestRank = line[1]
	ModelStructure = line[2]
	idCutoffs = line[3]
	gRichness = line[4]
	sqrt = line[5]

sys.stderr.write('\n### '+time.ctime(time.time())+': Classifying with '+HighestRank+' as highest rank and '+reStructure+' as lowest rank used for random intercept specification \n')

# Run Vsearch algnmnt
#subprocess.call(['vsearch', '--usearch_global', str(args.InputFasta), '--db', str(DBDIR+args.Database+'/DB.fa'), '--id', '0.6', \
#	'--maxaccepts', '100', '--maxrejects', '50', '--maxhits', '1', '--gapopen', '0TE', '--gapext', '0TE', '--userout', str(CTEMPDIR+'/Alnmt.txt'), \
#	'--userfields', 'query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits+qcov', '--query_cov', '0.95', '--threads', str(args.Threads)])
subprocess.call(['vsearch', '--usearch_global', str(args.InputFasta), '--db', str(DBDIR+args.Database+'/DB.fa'), '--id', '0.6', \
        '--maxaccepts', '100', '--maxrejects', '50', '--maxhits', '50', '--gapopen', '0TE', '--gapext', '0TE', '--userout', str(CTEMPDIR+'/Alnmt_2nd.txt'), \
        '--userfields', 'query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits+qcov', '--query_cov', '0.95', '--threads', str(args.Threads)])

# Get Vsearch Top hits
subprocess.call(['GetVsrchTopHit.py', str(CTEMPDIR+'/Alnmt_2nd.txt'), str(CTEMPDIR+'/Alnmt.txt')])

# Get Mtxa2 formatted output
subprocess.call(['VsearchToMetaxa2.py', '-v', str(CTEMPDIR+'/Alnmt.txt'), '-t', str(DBDIR+args.Database+'/DB.tax'), '-o', str(CTEMPDIR+'/tmp.tax')])

# reformat lineages  # not used in training
subprocess.call(['Get2ndHitTaxID_CLASSIFY.py', str(DBDIR+args.Database+'/DB.tax'), str(CTEMPDIR+'/Alnmt_2nd.txt'), str(CTEMPDIR+'/Alnmt_2nd.csv')])
subprocess.call(['FrmtLineages.py', str(CTEMPDIR+'/tmp.tax'), str(CTEMPDIR+'/Alnmt_2nd.csv'), str(CTEMPDIR+'/tmp2.tax')]) 

# Run GLMM analysis and output calls and probabilities    # replace str(CTEMPDIR+'/tmp2.tax') with str(CTEMPDIR+'/tmp3.tax')
subprocess.call(['subClassifyGLMM.r', str(CTEMPDIR+'/tmp2.tax'), str(DBDIR+args.Database+'/'), str(args.Output), HighestRank, reStructure, \
	ModelStructure, idCutoffs, gRichness, sqrt])

# Clean up tmp
if bool(args.SaveTemp) == False:
	subprocess.call(['rm', '-rf', str(CTEMPDIR)]) 
