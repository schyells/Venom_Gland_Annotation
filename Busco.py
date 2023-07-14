#! /usr/bin/env python # tells the command line to use python
import sys, os, shutil
import argparse
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from random import *
from Bio.Blast import NCBIXML
from Bio import SeqIO
# Set up and extract the command line options.
parser = argparse.ArgumentParser(description='Generates a .fasta file from the result table given from BUSCO. This file can then be used to search the a blast database made from the ancestral sequences given from BUSCO. This file can then be used in BUSCO2.py')
parser.add_argument("-i","--input",
                    type=str,
                    help="full_table.tsv")
parser.add_argument("-f","--fasta",
                    nargs='?',
                    type=str,
                    help="Name of trinity fasta file")
parser.add_argument("-a","--animal",
                    nargs='?',
                    type=str,
                    help="Name of trinity fasta file")
args = parser.parse_args()
Dic = {}
Infile = open(args.input, 'rU')
for line in Infile:
	if line[0] == '#':
		Fuck = "you"
	else:	
		line = line.split('\t')
		if line[1] == "Complete":
			Dic[line[2]] = line[8].strip()
record_dict = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))
for key, value in Dic.items():
	Name = key
	Prot = value
	seq = record_dict[Name]
	n= open(args.animal+"_busco.fasta",'a+')
	n.write('>'+Name+'\n'+str(seq.seq)+'\n')
	n.close()