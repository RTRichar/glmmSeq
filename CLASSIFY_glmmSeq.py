#!/usr/bin/env python

import subprocess, argparse, sys, time, os

parser = argparse.ArgumentParser(description="")
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('universal optional arguments')
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
DBDIR = str(os.path.abspath(os.path.dirname(sys.argv[0]))+'/')

# Run Vsearch algnmnt
subprocess.call(['vsearch', '--usearch_global', str(args.InputFasta), '--db', str(DBDIR+args.Database+'/DB.fa'), '--id', '0.75', \
	'--maxaccepts', '100', '--maxrejects', '50', '--maxhits', '5', '--gapopen', '0TE', '--gapext', '0TE', '--userout', str(CTEMPDIR+'/Alnmt.txt'), \
	'--userfields', 'query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits+qcov', '--query_cov', '0.95', '--threads', str(args.Threads)])

# Get Mtxa2 formatted output
subprocess.call(['VsearchToMetaxa2.py', '-v', str(CTEMPDIR+'/Alnmt.txt'), '-t', str(DBDIR+args.Database+'/DB.tax'), '-o', str(CTEMPDIR+'/tmp.tax')])

# reformat lineages
subprocess.call(['FrmtLineages.py', str(CTEMPDIR+'/tmp.tax'), str(CTEMPDIR+'/tmp2.tax')])

# Run GLMM analysis and output calls and probabilities
subprocess.call(['ClassifyGLMM.r', str(CTEMPDIR+'/tmp2.tax'), str(DBDIR+args.Database+'/'), str(args.Output)])

# Clean up tmp
if bool(args.SaveTemp) == False:
	subprocess.call(['rm', '-rf', str(CTEMPDIR)]) 
