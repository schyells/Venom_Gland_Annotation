#! /usr/bin/env python # tells the command line to use python
#imports the module that allows for user input into the program
import re
import datetime
import textwrap
import argparse
import sys, os, shutil
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from random import *
from Bio.Blast import NCBIXML
from Bio import SeqIO
now = datetime.datetime.now()

# Set up and extract the command line options.
parser = argparse.ArgumentParser(description='''Automated annotation of transcripts.
Version Apr 14 2020.
Will return a fasta file of toxin sequences that can be used in the AutoNamer script and a protein fasta file of toxin sequences that can be used in protein searches. Output on the screen is the contig names with the start stop and signal peptide locations''')

parser.add_argument("-i","--input",
                    type=str,
                    help="Fasta file of contigs with names the exact same as xml file")
parser.add_argument("-x","--xml",
                    type=str,
                    help="XML file of BLAST results from a protein based search")
parser.add_argument("-o","--output",
                    nargs='?',
                    type=str,
                    default="Annotations",
                    help="Name of fasta output file without the extention")
parser.add_argument("-s","--sig",
                    nargs='?',
                    type=str,
                    default="Y",
                    help="Use a Y or a N to indicate if you want to use signalp searches or not")
args = parser.parse_args()


#Reverse transcribes the sequence	
def reverse_sequence(s):
	"""Return a reverse sequence of `s`"""
	s = s.replace('G','B')
	s = s.replace('C','G')
	s = s.replace('B','C')
	s = s.replace('T','U')
	s = s.replace('A','T')
	s = s.replace('U','A')
	return(s)
	
	
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

#Finds the signal Peptide
def signalp(p):
	f = open("signalp.fasta",'w+')
	f.write('>P'+'\n'+p[:-1])
	f.close()
	sp = 'signalp -u 0.34 -U 0.34 signalp.fasta  > signalp.txt'
	os.system(sp)
	cut = 0
	file = open('signalp.txt', 'rU')
	for line in file:
		if line[0] == '#':
			line = line
		else:
			line = line.strip('\n')
			signalre = r'([NY])'
			signal = re.search(signalre,line)
			sig = signal.group(1)
			if sig == "Y":
				line = line.replace(' ','l')
				cut = line[34:36]
				cut = cut.replace('l','')
	return(sig,cut)

# Finds potential starts
def annotate(a):
	bestre = r'(M\w+)'
	best = re.search(bestre,a)
	best = best.group(1)
	return(best)
	
# Finds correct start postion for the output file	
def startpos(P, n):
    x = P.find('M')
    while x >= 0 and n > 1:
        x = P.find('M', x+1)
        n -= 1
    return x


#Looks for the reading frame that matches the blast result and returns the longest protein or the protein with a signal peptide from that reading frame and the correct position. 
blast = open(args.xml,"r")
blast_records = NCBIXML.parse(blast)
List=[]
for rec in blast_records:
	for seq_record in SeqIO.parse(args.input, "fasta"):
		if rec.query.split()[0] == seq_record.id:
			aligncount = 1
			for alignment in rec.alignments:
				if aligncount == 1:
					aligncount +=1
					hspcount = 1
					for hsp in alignment.hsps:
						if hspcount == 1:
							hspcount +=1
							count=0
							if hsp.frame[0] < 0 :
								start = hsp.query_end
								end = hsp.query_start
								Tstart = start 
								Tend= end
							else :
								start = hsp.query_start
								end = hsp.query_end
								Tstart = start 
								Tend= end
							if start > end:
								#print ('Backwards',Tstart,Tend,start,end)
								start = (len(str(seq_record.seq))-start+1)
								start = (start - (int(start/3)*3))
								if start == 0:
									start = 3
								end = (len(str(seq_record.seq))-end+4)
								seq = str(reverse_sequence(str(seq_record.seq[::-1])))
							else:
								start = (start - (int(start/3)*3))
								seq = seq_record.seq
								if start == 0:
									start = 3
							if int(Tend) < int(Tstart):
								Tstart = (len(str(seq_record.seq))-Tstart+1)
								Tend = (len(str(seq_record.seq))-Tend+4)
							if len(seq[start-1:]) % 3 == 0:
								CDS = seq[start-1:]
							elif len(seq[start-1:-1]) % 3 == 0:
								CDS = seq[start-1:-1]
							elif len(seq[start-1:-2]) % 3 == 0:
								CDS = seq[start-1:-2]
							prot = translate(str(CDS))
							prot2 = prot
							while prot2.count('_') >= 1:
								num = prot2.find('_')
								prot3 = prot2[:num+1]
								if prot3.count('M') == 0:
									prot2 = prot2[num+1:]
								else:
									while prot3.count('M') >=1:
										mcount= 0
										while count == 0:
											count = 0
											if prot3.count('M') == 0:
												break
											if prot3.count('M') ==1 and prot3[-2] =='M':
												break
											prot3 = annotate(prot3)
											mcount+=1
											mstart = startpos(prot2,mcount)
											new=prot2[mstart:] 
											b = (len(prot[:])-len(new))*3
											#print(Tstart, len(prot3), b, Tend,start,len(new),len(prot))
											if Tstart <= ((len(prot3)*3)+start+b-1) and Tstart >= (start+b):
												if (((Tend-(start+b)) - (len(prot3)*3))/(Tend-(start+b))) <= 0.1:
													mstart = startpos(prot2,mcount)
													new=prot2[mstart:]
													if args.sig == "Y" or args.sig =="y":
														sig = signalp(prot3)
														if sig[0] == 'Y':
															n= open(args.output+".fasta",'a+')
															n.write('>'+rec.query.split()[0]+'_str'+str(start+b)+'_stp'+str((len(prot3)*3)+start+b-1)+'_sp'+str(int(sig[1])-1)+'\n'+str(seq)+'\n')
															n.close()
															m = open(args.output+'_prot.fasta','a+')
															m.write('>'+rec.query.split()[0]+'\n'+prot3+'\n')
															m.close()
															print('>'+rec.query.split()[0]+'_str'+str(start+b)+'_stp'+str((len(prot3)*3)+start+b-1)+'_sp'+str(int(sig[1])-1))
															count = 1
														else:
															prot3 = prot3[1:]
													elif args.sig =='N' or args.sig =='n':
														n= open(args.output+".fasta",'a+')
														n.write('>'+rec.query.split()[0]+'_str'+str(start+b)+'_stp'+str((len(prot3)*3)+start+b-1)+'_sp0'+'\n'+str(seq)+'\n')
														n.close()
														m = open(args.output+'_prot.fasta','a+')
														m.write('>'+rec.query.split()[0]+'\n'+prot3+'\n')
														m.close()
														print('>'+rec.query.split()[0]+'_str'+str(start+b)+'_stp'+str((len(prot3)*3)+start+b-1)+'_sp0')
														count = 1
												else:
													mstart = startpos(prot2,mcount)
													new=prot2[mstart:]
													if args.sig == "Y" or args.sig =="y":
														sig = signalp(prot3)
														if sig[0] == 'Y':
															n= open(args.output+"_partial.fasta",'a+')
															n.write('>'+rec.query.split()[0]+'_str'+str(start+b)+'_stp'+str((len(prot3)*3)+start+b-1)+'_sp'+str(int(sig[1])-1)+'\n'+str(seq)+'\n')
															n.close()
															m = open(args.output+'_partial_prot.fasta','a+')
															m.write('>'+rec.query.split()[0]+'\n'+prot3+'\n')
															m.close()
															print('>partial'+rec.query.split()[0]+'_str'+str(start+b)+'_stp'+str((len(prot3)*3)+start+b-1)+'_sp'+str(int(sig[1])-1))
															count = 1
														else:
															prot3 = prot3[1:]
													elif args.sig =='N' or args.sig =='n':
														n= open(args.output+"_partial.fasta",'a+')
														n.write('>'+rec.query.split()[0]+'_str'+str(start+b)+'_stp'+str((len(prot3)*3)+start+b-1)+'_sp0'+'\n'+str(seq)+'\n')
														n.close()
														m = open(args.output+'_partial_prot.fasta','a+')
														m.write('>'+rec.query.split()[0]+'\n'+prot3+'\n')
														m.close()
														print('>partial'+rec.query.split()[0]+'_str'+str(start+b)+'_stp'+str((len(prot3)*3)+start+b-1)+'_sp0')
														count = 1		
											else:
												prot3 = prot3[1:]
										if prot3.count('M') ==1 and prot3[-2] =='M':
											break
										if count == 1:
											break
									prot2 = prot2[num+1:]
if args.sig == "Y" or args.sig =="y":
	fin = 'rm signalp.*'
	os.system(fin)
