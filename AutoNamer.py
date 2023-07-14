#! /usr/bin/env python # tells the command line to use python
import sys #imports the module that allows for user input into the program
import re
import os
import argparse
import datetime
import textwrap
from Bio import SeqIO
from Bio.Seq import Seq
from random import *
from Bio.Blast import NCBIXML
from Bio import SeqIO
Usage = '''This program requires a fasta file and name
A blast file .xml
command should be AutoNamer.py file.fasta file.xml 
July 14 2023'''
now = datetime.datetime.now()


# Set up and extract the command line options.
parser = argparse.ArgumentParser(description='''Automated naming of fasta file.
Version July 14 2023.
takes a fasta file from Blast_Namer and an xml file from the same contigs to generate a .seq file that includes the cds, protein and signal peptide annotated the file name will include part of the name of the blast hit along with animal name and contig name''')

parser.add_argument("-i","--input",
                    type=str,
                    help="Fasta file of contigs with names the exact same as xml file")
parser.add_argument("-x","--xml",
                    type=str,
                    help="XML file of BLAST results")
parser.add_argument("-o","--output",
                    nargs='?',
                    type=str,
                    default="Annotations",
                    help="Name of file to hold all of the .seq files generated here")
parser.add_argument("-a","--animal",
                    nargs='?',
                    type=str,
                    default="Animal",
                    help="Name of file to hold all of the .seq files generated here")                    
parser.add_argument("-s","--sig",
                    nargs='?',
                    type=str,
                    default="Y",
                    help="Use a Y or a N to indicate if you want to use signalp searches or not")
args = parser.parse_args()


#Translates Sequence
def translate(seq): 

	ambiguities=['Y','R','W','S','K','M','N','D','V','H','B']
	table = { 
		'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
		'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
		'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
		'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',				 
		'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
		'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
		'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
		'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
		'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
		'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
		'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
		'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
		'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
		'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
		'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
		'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
	} 
	protein ="" 
	if len(seq)%3 == 0: 
		for i in range(0, len(seq), 3): 
			codon = seq[i:i + 3] 
			if any(base in codon for base in ambiguities):
				protein+='X'
			else:
				protein+= table[codon] 
	return protein 



FASTA = args.input # reads txt
xml = args.xml #reads fasta
animal = args.animal
blast = open(xml,"r")
blast_records = NCBIXML.parse(blast)
List=[]
for rec in blast_records: 
	for seq_record in SeqIO.parse(FASTA, "fasta"):
		filese = r'(.+)_str(\d+)_stp(\d+)_sp(\d+)'
		filere = re.search(filese, str(seq_record.id))
		file = filere.group(1)
		start = int(filere.group(2))
		stop = int(filere.group(3))
		sp = int(filere.group(4))
		query = file
		if rec.query.split()[0] == query:
			aligncount = 1
			for alignment in rec.alignments:
				if aligncount == 1:
					title = alignment.title	
					aligncount +=1
					titlere = r'.+ TSA: (.+)'
					namere = re.search(titlere,title) #searches through the length string according to the search query above
					name = title
					sequence = str(seq_record.seq) 
					cds = str(sequence[start-1:(stop-3)])
					translation = (translate(str(cds)))
					translation1 = translation[0:43]
					translation2 = textwrap.wrap(translation[43:],58)
					print(file)
					n=open(animal+'_'+name[18:50].replace(' ','')+'_'+file.replace('_','-')+'.seq','w+')
					if args.sig == "Y" or args.sig =="y":
						n.write('LOCUS       '+animal+'_'+name+'_'+file.replace('_','-')+'        '+str(len(str(seq_record.seq)))+' bp    DNA     linear       '+now.strftime("%d-%b-%y")+'\n'+file+'\n'+'FEATURES             Location/Qualifiers'+'\n'+'     CDS             '+str(start)+'..'+str(stop)+'\n'+'                     /note="'+name+'"\n'+'                     /translation="'+translation1 +'\n                     '+'\n                     '.join(translation2)+'"\n     sig_peptide     '+str(start)+'..'+str(start-1+(sp*3))+'\n                     /note="Signal Peptide"'+'\n     source          1..'+str(len(str(seq_record.seq)))+'\n                     /dnas_title="'+animal+'_'+name+'_'+file.replace('_','-')+'"'+'\n^^\n'+str(seq_record.seq))
					else:
						n.write('LOCUS       '+animal+'_'+name+'_'+file.replace('_','-')+'        '+str(len(str(seq_record.seq)))+' bp    DNA     linear       '+now.strftime("%d-%b-%y")+'\n'+file+'\n'+'FEATURES             Location/Qualifiers'+'\n'+'     CDS             '+str(start)+'..'+str(stop)+'\n'+'                     /note="'+name+'"\n'+'                     /translation="'+translation1 +'\n                     '+'\n                     '.join(translation2)+'\n     source          1..'+str(len(str(seq_record.seq)))+'\n                     /dnas_title="'+animal+'_'+name+'_'+file.replace('_','-')+'"'+'\n^^\n'+str(seq_record.seq))
					n.close()
					List.append(query)
for seq_record in SeqIO.parse(FASTA, "fasta"):
	filese = r'(.+)_str(\d+)_stp(\d+)_sp(\d+)'
	filere = re.search(filese, str(seq_record.id))
	file = filere.group(1)
	start = int(filere.group(2))
	stop = int(filere.group(3))
	sp = int(filere.group(4))
	query = file
	if not query in List:
			sequence = str(seq_record.seq) 
			cds = str(sequence[start-1:(stop-3)])
			translation = (translate(str(cds)))
			translation1 = translation[0:43]
			translation2 = textwrap.wrap(translation[43:],58)
			print(file)
			n=open(animal+'_VP_'+file.replace('_','-')+'.seq','w+')
			if args.sig == "Y" or args.sig =="y":
				n.write('LOCUS       '+animal+'_VP_'+file.replace('_','-')+'        '+str(len(str(seq_record.seq)))+' bp    DNA     linear       '+now.strftime("%d-%b-%y")+'\n'+file+'\n'+'FEATURES             Location/Qualifiers'+'\n'+'     CDS             '+str(start)+'..'+str(stop)+'\n'+'                     /note="'+'VP'+'"\n'+'                     /translation="'+translation1 +'\n                     '+'\n                     '.join(translation2)+'"\n     sig_peptide     '+str(start)+'..'+str(start-1+(sp*3))+'\n                     /note="Signal Peptide"'+'\n     source          1..'+str(len(str(seq_record.seq)))+'\n                     /dnas_title="'+animal+'_'+'VP'+'_'+file.replace('_','-')+'"'+'\n^^\n'+str(seq_record.seq))
			else:
				n.write('LOCUS       '+animal+'_VP_'+file.replace('_','-')+'        '+str(len(str(seq_record.seq)))+' bp    DNA     linear       '+now.strftime("%d-%b-%y")+'\n'+file+'\n'+'FEATURES             Location/Qualifiers'+'\n'+'     CDS             '+str(start)+'..'+str(stop)+'\n'+'                     /note="'+'VP'+'"\n'+'                     /translation="'+translation1 +'\n                     '+'\n                     '.join(translation2)+'\n     source          1..'+str(len(str(seq_record.seq)))+'\n                     /dnas_title="'+animal+'_'+'VP'+'_'+file.replace('_','-')+'"'+'\n^^\n'+str(seq_record.seq))
			n.close()
fin = 'mkdir '+args.output
os.system(fin)
fin = 'mv *.seq '+args.output
os.system(fin)	