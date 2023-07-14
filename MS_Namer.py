#! /usr/bin/env python # tells the command line to use python

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


# Set up and extract the command line options.
parser = argparse.ArgumentParser(description='''Automated annotation of transcripts from scaffold.
Version March 16 2021.
Will return a fasta file of toxin sequences that can be used in the AutoNamer script and a protein fasta file of toxin sequences that can be used in protein searches. Output on the screen is the contig names with the start stop and signal peptide locations''')

parser.add_argument("-i","--input",
                    type=str,
                    help="Fasta file of contigs with names the exact same as scaffold file")
parser.add_argument("-scaff","--scaffold",
                    type=str,
                    help="Scaffold Protein Accession Number report")
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

# Runs Signal P and returns whether the protein has a signal peptide or not
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

# rewrites the fasta and the txt file so that it only contains the information needed to process the files
def geneious():
	Infile = open(csv, 'rU')
	Infile2 = open(fastaL, 'rU')
	count = 0
	fasta = Infile2.read()
	genList=[]
	for line in Infile:
		if count < 2:
			count +=1
		elif line[0:3] != 'END' and line[0:3] != 'AMP' and line[0:3] != 'FIM' and line[0:3] != 'HCH':
			x = line.split('\t')[0]
			y = line.split('\t')[1]
			count2 = 1
			for t in x[::-1]:
				if t == '_':
					break
				else:
					count2 +=1
			contig = x[:(0-count2)]
			if contig not in genList:
				genList.append(contig)
				for seq_record in SeqIO.parse(fastaL, "fasta"):
					if seq_record.id == contig:
						filtered = open(("Filtered.fasta"),'a+')
						filtered.write('>'+seq_record.id+'\n'+str(seq_record.seq)+'\n')
						filtered.close()
						break
				location = open(("location.txt"),'a+')
				location.write(x+'\t'+y+'\n')
				location.close()
			else:
				print(contig,y, 'Multiple CDS in same Contig Need to Check by hand')
	Infile.close()
	Infile2.close()



csv = args.scaffold # reads txt
fastaL = args.input #reads fasta
geneious()
CSV = "location.txt"
FASTA = "Filtered.fasta"
InCSV = open(CSV,'rU')
CSVList= []
Dict = {}
for line in InCSV:
	line = line.strip('\n')
	MyList=line.split('\t')
	Newline = MyList[0]
	count = 1
	for t in Newline[::-1]:
		if t == '_':
			break
		else:
			count +=1
	count2 = 0
	length2 = ""
	length = MyList[1]
	length = length.replace('"', '')
	length = length.strip()
	for t in length:
		if t == "]":
			break
		else:
			count2 +=1
			length2 +=t			
	Dict[(Newline[:(0-count)])] = length2[1:]
	
for seq_record in SeqIO.parse(FASTA, "fasta"):
	length = Dict[seq_record.id]
	sequence = str(seq_record.seq)		
	length = length.strip() #removes spaces from string
	numre = r'(\d+)( - )(\d+)' #define search using grep into distinct groups
	num = re.search(numre,length) #searches through the length string according to the search query above
	first = int(num.group(1))
	second = int(num.group(3))
	if second < 4 or second >= len(seq_record.seq)-3:
		f = open("MSlog.txt",'a+')
		f.write('bad '+seq_record.id+' Truncated Sequence'+'\n')
		f.close()
	else:
		if first > second:
			cds = sequence[(second-4):first]
			cds = cds[::-1]
			cds = str(reverse_sequence(str(cds)))
		else:
			cds = str(sequence[first-1:(second+3)])
		prot = (translate(str(cds)))
		if len(cds)%3 != 0:
			f = open("MSlog.txt",'a+')
			f.write('bad '+seq_record.id+' Not a multiple of 3'+'\n')
			f.close()
		else:	
			if prot[-1] == '_' and prot[0] =='M': #selects for stop and start sequences
				if args.sig == "Y" or args.sig =="y":
					sig = signalp(prot)			
					if sig[0] == 'Y':
						print(seq_record.id+'\t'+Dict[seq_record.id]+'\t'+str(int(sig[1])-1))
						j = open("MS_Namer.fasta",'a+')
						j.write('>'+seq_record.id+'\n'+cds+'\n')
						j.close()
						t = open("MS_Namer_prot.fasta",'a+')
						t.write('>'+seq_record.id+'\n'+prot[:-1]+'\n')
						t.close() 
						n= open("Annotations.fasta",'a+') 
						if first > second:
							n.write('>' + seq_record.id+'_str'+str(len(str(seq_record.seq))-first+1)+'_stp'+str(len(str(seq_record.seq))-second+4)+'_sp'+str(int(sig[1])-1)+ '_m0\n' + str(reverse_sequence(str(seq_record.seq[::-1]))) + '\n')
						else:	
							n.write('>' + seq_record.id+'_str'+str(first)+'_stp'+str(second+3)+'_sp'+str(int(sig[1])-1)+ '_m0\n' + str(seq_record.seq) + '\n') 
						n.close()
					else:
						if prot[1:-2].count('M')==0:
							f = open("MSlog.txt",'a+')
							f.write('bad '+seq_record.id+' Missing signal peptide\n')
							f.close()
						else:
							startre = r'(M\w+)'
							start = re.search(startre,prot[1:-1])
							start = start.group(1)
							sig2 = signalp(start)	
							if sig2[0] == 'Y':
								print(seq_record.id+'\t'+Dict[seq_record.id]+'\t'+str(int(sig2[1])-1)+'\t'+'moved one')
								a = (len(prot[:-1])-len(start))*3
								j = open("MS_Namer.fasta",'a+')
								j.write('>'+seq_record.id+' moved 1'+'\n'+cds[a:]+'\n')
								j.close()
								t = open("MS_Namer_prot.fasta",'a+')
								t.write('>'+seq_record.id+'\n'+start+'\n')
								t.close() 
								n= open("Annotations.fasta",'a+')
								if first > second:
									n.write('>' + seq_record.id+'_str'+str(len(str(seq_record.seq))-(first-a)+1)+'_stp'+str(len(str(seq_record.seq))-second+4)+'_sp'+str(int(sig2[1])-1)+ '_m1\n' + str(reverse_sequence(str(seq_record.seq[::-1]))) + '\n')
								else:	 
									n.write('>' + seq_record.id+'_str'+str(first+a)+'_stp'+str(second+3)+'_sp'+str(int(sig2[1])-1)+ '_m1\n' + str(seq_record.seq) + '\n') 
								n.close() 
							else:	
								if start[1:-1].count('M') == 0:
									f = open("MSlog.txt",'a+')
									f.write('bad '+seq_record.id+' Missing signal peptide\n')
									f.close()
								else:
									startre2 = r'(M\w+)'
									start2 = re.search(startre2,start[1:])
									start = start2.group(1)
									sig3 = signalp(start2.group(1))
									if sig3[0] == 'Y':
										print(seq_record.id+'\t'+Dict[seq_record.id]+'\t'+str(int(sig3[1])-1)+'\t'+'moved two')
										a = (len(prot[:-1])-len(start2.group(1)))*3
										j = open("MS_Namer.fasta",'a+')
										j.write('>'+seq_record.id+' moved 2'+'\n'+cds[a:]+'\n')
										j.close()
										t = open("MS_Namer_prot.fasta",'a+')
										t.write('>'+seq_record.id+'\n'+start+'\n')
										t.close() 
										n= open("Annotations.fasta",'a+') 
										if first > second:
											n.write('>' + seq_record.id+'_str'+str(len(str(seq_record.seq))-(first-a)+1)+'_stp'+str(len(str(seq_record.seq))-second+4)+'_sp'+str(int(sig3[1])-1)+ '_m2\n' + str(reverse_sequence(str(seq_record.seq[::-1]))) + '\n')
										else:
											n.write('>' + seq_record.id+'_str'+str(first+a)+'_stp'+str(second+3)+'_sp'+str(int(sig3[1])-1)+ '_m2\n' + str(seq_record.seq) + '\n') 
										n.close() 
									else:	
										if start[1:-1].count('M') == 0:
											f = open("MSlog.txt",'a+')
											f.write('bad '+seq_record.id+' Missing signal peptide\n')
											f.close()
										else:
											startre3 = r'(M\w+)'
											start3 = re.search(startre3,start[1:])
											start = start3.group(1)
											sig4 = signalp(start3.group(1))
											if sig4[0] == 'Y':
												print(seq_record.id+'\t'+Dict[seq_record.id]+'\t'+str(int(sig4[1])-1)+'\t'+'moved three')
												a = (len(prot[:-1])-len(start3.group(1)))*3
												j = open("MS_Namer.fasta",'a+')
												j.write('>'+seq_record.id+' moved 3'+'\n'+cds[a:]+'\n')
												j.close()
												t = open("MS_Namer_prot.fasta",'a+')
												t.write('>'+seq_record.id+'\n'+start+'\n')
												t.close() 
												n= open("Annotations.fasta",'a+') 
												if first > second:
													n.write('>' + seq_record.id+'_str'+str(len(str(seq_record.seq))-(first-a)+1)+'_stp'+str(len(str(seq_record.seq))-second+4)+'_sp'+str(int(sig4[1])-1)+ '_m3\n' + str(reverse_sequence(str(seq_record.seq[::-1]))) + '\n')
												else:
													n.write('>' + seq_record.id+'_str'+str(first+a)+'_stp'+str(second+3)+'_sp'+str(int(sig4[1])-1)+ '_m3\n' + str(seq_record.seq) + '\n') 
												n.close() 
											else:
												if start[1:-1].count('M') == 0:
													f = open("MSlog.txt",'a+')
													f.write('bad '+seq_record.id+' Missing signal peptide\n')
													f.close()
												else:
													startre4 = r'(M\w+)'
													start4 = re.search(startre4,start[1:])
													start = start4.group(1)
													sig5 = signalp(start4.group(1))
													if sig5[0] == 'Y':
														print(seq_record.id+'\t'+Dict[seq_record.id]+'\t'+str(int(sig5[1])-1)+'\t'+'moved four')
														a = (len(prot[:-1])-len(start4.group(1)))*3
														j = open("MS_Namer.fasta",'a+')
														j.write('>'+seq_record.id+' moved 4'+'\n'+cds[a:]+'\n')
														j.close()
														t = open("MS_Namer_prot.fasta",'a+')
														t.write('>'+seq_record.id+'\n'+start+'\n')
														t.close() 
														n= open("Annotations.fasta",'a+')
														if first > second:
															n.write('>' + seq_record.id+'_str'+str(len(str(seq_record.seq))-(first-a)+1)+'_stp'+str(len(str(seq_record.seq))-second+4)+'_sp'+str(int(sig5[1])-1)+ '_m4\n' + str(reverse_sequence(str(seq_record.seq[::-1]))) + '\n')
														else:
															n.write('>' + seq_record.id+'_str'+str(first+a)+'_stp'+str(second+3)+'_sp'+str(int(sig5[1])-1)+ '_m4\n' + str(seq_record.seq) + '\n') 
														n.close() 
													else:
														if start[1:-1].count('M') == 0:
															f = open("MSlog.txt",'a+')
															f.write('bad '+seq_record.id+' Missing signal peptide\n')
															f.close()
														else:
															startre5 = r'(M\w+)'
															start5 = re.search(startre5,start[1:])
															start = start5.group(1)
															sig6 = signalp(start5.group(1))
															if sig6[0] == 'Y':
																print(seq_record.id+'\t'+Dict[seq_record.id]+'\t'+str(int(sig6[1])-1)+'\t'+'moved five')
																a = (len(prot[:-1])-len(start5.group(1)))*3
																j = open("MS_Namer.fasta",'a+')
																j.write('>'+seq_record.id+' moved 5'+'\n'+cds[a:]+'\n')
																j.close()
																t = open("MS_Namer_prot.fasta",'a+')
																t.write('>'+seq_record.id+'\n'+start+'\n')
																t.close() 
																n= open("Annotations.fasta",'a+') 
																if first > second:
																	n.write('>' + seq_record.id+'_str'+str(len(str(seq_record.seq))-(first-a)+1)+'_stp'+str(len(str(seq_record.seq))-second+4)+'_sp'+str(int(sig6[1])-1)+ '_m5\n' + str(reverse_sequence(str(seq_record.seq[::-1]))) + '\n')
																else:
																	n.write('>' + seq_record.id+'_str'+str(first+a)+'_stp'+str(second+3)+'_sp'+str(int(sig6[1])-1)+ '_m3\n' + str(seq_record.seq) + '\n') 
																n.close() 
															else:
																if start[1:-1].count('M') == 0:
																	f = open("MSlog.txt",'a+')
																	f.write('bad '+seq_record.id+' Missing signal peptide\n')
																	f.close()
																else:
																	startre6 = r'(M\w+)'
																	start6 = re.search(startre6,start[1:])
																	start = start6.group(1)
																	sig7 = signalp(start6.group(1))
																	if sig7[0] == 'Y':
																		print(seq_record.id+'\t'+Dict[seq_record.id]+'\t'+str(int(sig7[1])-1)+'\t'+'moved six')
																		a = (len(prot[:-1])-len(start6.group(1)))*3
																		j = open("MS_Namer.fasta",'a+')
																		j.write('>'+seq_record.id+' moved 6'+'\n'+cds[a:]+'\n')
																		j.close()
																		t = open("MS_Namer_prot.fasta",'a+')
																		t.write('>'+seq_record.id+'\n'+start+'\n')
																		t.close() 
																		n= open("Annotations.fasta",'a+') 
																		if first > second:
																			n.write('>' + seq_record.id+'_str'+str(len(str(seq_record.seq))-(first-a)+1)+'_stp'+str(len(str(seq_record.seq))-second+4)+'_sp'+str(int(sig7[1])-1)+ '_m6\n' + str(reverse_sequence(str(seq_record.seq[::-1]))) + '\n')
																		else:
																			n.write('>' + seq_record.id+'_str'+str(first+a)+'_stp'+str(second+3)+'_sp'+str(int(sig7[1])-1)+ '_m3\n' + str(seq_record.seq) + '\n') 
																		n.close() 
																	else:
																		if start[1:-1].count('M') == 0:
																			f = open("MSlog.txt",'a+')
																			f.write('bad '+seq_record.id+' Missing signal peptide\n')
																			f.close()
																		else:
																			startre7 = r'(M\w+)'
																			start7 = re.search(startre7,start[1:])
																			start = start7.group(1)
																			sig8 = signalp(start7.group(1))
																			if sig8[0] == 'Y':
																				print(seq_record.id+'\t'+Dict[seq_record.id]+'\t'+str(int(sig8[1])-1)+'\t'+'moved seven')
																				a = (len(prot[:-1])-len(start7.group(1)))*3
																				j = open("MS_Namer.fasta",'a+')
																				j.write('>'+seq_record.id+' moved 7'+'\n'+cds[a:]+'\n')
																				j.close()
																				t = open("MS_Namer_prot.fasta",'a+')
																				t.write('>'+seq_record.id+'\n'+start+'\n')
																				t.close() 
																				n= open("Annotations.fasta",'a+')
																				if first > second:
																					n.write('>' + seq_record.id+'_str'+str(len(str(seq_record.seq))-(first-a)+1)+'_stp'+str(len(str(seq_record.seq))-second+4)+'_sp'+str(int(sig8[1])-1)+ '_m7\n' + str(reverse_sequence(str(seq_record.seq[::-1]))) + '\n')
																				else:
																					n.write('>' + seq_record.id+'_str'+str(first+a)+'_stp'+str(second+3)+'_sp'+str(int(sig8[1])-1)+ '_m3\n' + str(seq_record.seq) + '\n') 
																				n.close() 
																			else:
																				f = open("MSlog.txt",'a+')
																				f.write('bad '+seq_record.id+' Missing signal peptide\n')
																				f.close()
				else:
					print(seq_record.id+'\t'+Dict[seq_record.id]+'\t 0')
					j = open("MS_Namer.fasta",'a+')
					j.write('>'+seq_record.id+'\n'+cds+'\n')
					j.close()
					t = open("MS_Namer_prot.fasta",'a+')
					t.write('>'+seq_record.id+'\n'+prot[:-1]+'\n')
					t.close() 
					n= open("Annotations.fasta",'a+') 
					if first > second:
						n.write('>' + seq_record.id+'_str'+str(len(str(seq_record.seq))-first+1)+'_stp'+str(len(str(seq_record.seq))-second+4)+'_sp0'+ '_m0\n' + str(reverse_sequence(str(seq_record.seq[::-1]))) + '\n')
					else:	
						n.write('>' + seq_record.id+'_str'+str(first)+'_stp'+str(second+3)+'_sp0'+ '_m0\n' + str(seq_record.seq) + '\n') 
					n.close()																			
			else:
				f = open("MSlog.txt",'a+')
				f.write('bad '+seq_record.id+' Missing M or Stop\n')
				f.close()	
fin = 'rm signalp.*'
os.system(fin)

		
		
		
		
		
	
		
		
		
		
		
		
		
		