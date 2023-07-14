#!/usr/bin/env python3

# ChimeraCheck is designed to screen CDS of annotated transcripts for evidence of a shitty assembly. 
# It looks for coverage anomalies.

# Additional software necessary to run this:
# (1) bwa 
# (2) GATK
# (3) bedtools

# The input required includes:
# (1) a fasta file of the contigs to annotate
# (2) a fasta file of the nt sequences of previously annotated coding sequences

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import sys, os, shutil
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq

# Set up and extract the command line options.
parser = argparse.ArgumentParser(description='SE version 12/11/19')
parser.add_argument("-i","--input",
                    type=argparse.FileType('r'),
                    help="Fasta file of contigs to check. Use only CODING SEQUENCES.")
parser.add_argument("-r","--reads",
                    type=argparse.FileType('r'),
                    help="Fastq file of UNPAIRED reads.")
parser.add_argument("-b","--bwa",
                    nargs='?',
                    type=str,
                    default="bwa",
                    help="Path to bwa executable. Default assumes it is in your PATH.")
parser.add_argument("-s","--samtools",
                    nargs='?',
                    type=str,
                    default="samtools",
                    help="Path to samtools executable. Default assumes it is in your PATH.")
parser.add_argument("-bt","--bedtools",
                    nargs='?',
                    type=str,
                    default="bedtools",
                    help="Path to bedtools executable. Default assumes it is in your PATH.")
parser.add_argument("-p","--picard",
                    nargs='?',
                    type=str,
                    default="picard ",
                    help="Picard command. Something like: java -jar /PATH/TO/PICARD/picard.jar")
parser.add_argument("-g","--gatk",
                    nargs='?',
                    type=str,
                    default="java -jar /home/sellsworth/scripts/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar ",
                    help="GATK command. Some like: java -jar /PATH/TO/GATK/GenomeAnalysisTK.jar")
parser.add_argument("-m","--mismatches",
                    nargs='?',
                    type=int,
                    default=0,
                    help="Number of allowable mismatches to keep in alignment. The default is 0.")
parser.add_argument("-ts","--tooShort",
                    nargs='?',
                    type=int,
                    default=150,
                    help="Minimum length of clipped reads. Default is 150.")
parser.add_argument("-mq","--mapQuality",
                    nargs='?',
                    type=int,
                    default=0,
                    help="Minimum mapping quality. Default is 0. Note that reads with multiple mappings get assigned a quality of 0 by bwa.")
parser.add_argument("-min","--minRead",
                    nargs='?',
                    type=int,
                    default=150,
                    help="Minimum read length. Default is 150.")
parser.add_argument("-max","--maxRead",
                    nargs='?',
                    type=int,
                    default=1000,
                    help="Maximum read length. Default is 1000.")
parser.add_argument("-X","--foldCoverage",
                    nargs='?',
                    type=float,
                    default=20.0,
                    help="Maximum allowed fold coverage differential. Default is 20.")
parser.add_argument("-mer","--merSize",
                    nargs='?',
                    type=int,
                    default=100,
                    help="Mer size for sequence comparisons. Default is 100.")
parser.add_argument("-win","--window",
                    nargs='?',
                    type=int,
                    default=151,
                    help="Window size for coverage comparisons. Default is 151.")
parser.add_argument("-pdf","--pdfprint",
                    nargs='?',
                    type=str,
                    default="Y",
                    help="Use a Y or a N to indicate if you want to print pdfs of the good and bad coverage sequences. N will only print bad coverage sequences")                    
args = parser.parse_args()

columns =20


print("*"*int(columns))

# Check the input fasta file and output the number of detected sequences.
totalContigs = 0
for seq in SeqIO.parse(args.input,format="fasta") :
    totalContigs += 1
print("Total number of input contigs = " + str(totalContigs))
print("Maximum allowed mismatches per read retained in the alignment = " + str(args.mismatches))

print("*"*int(columns))

# Run the alignment
grepNM = "grep -E 'NM:i:[0-" + str(args.mismatches) + "][[:space:]]|^@'"
name = args.input.name.split(".")[0] 
print("Generating bwa alignment: " + name + ".bam")
# Build the bwa index
command = args.bwa + " index " + args.input.name
subprocess.call(command,shell=True)
# Generate the initial sam alignment
command = args.bwa + " mem -M -t 8 -R \'@RG\\tID:" + args.input.name + "\\tSM:" + args.reads.name + "' " + args.input.name + " " + args.reads.name + " | " + grepNM  + " > tmp1.sam"
subprocess.call(command,shell=True)
# Create a sorted bam file
command = args.picard + " SortSam INPUT=tmp1.sam OUTPUT=tmp2.bam SORT_ORDER=coordinate" 
subprocess.call(command,shell=True)
command = args.picard + " BuildBamIndex INPUT=tmp2.bam"
subprocess.call(command,shell=True)
# Remove overclipped reads
command = args.picard + "CreateSequenceDictionary REFERENCE=" + args.input.name + " OUTPUT=" + args.input.name.split(".")[0] + ".dict"
subprocess.call(command,shell=True)
command = args.samtools + " faidx " + args.input.name
subprocess.call(command,shell=True)
command = args.gatk + "-T PrintReads -R " + args.input.name + " -I tmp2.bam  -rf OverclippedRead --filter_is_too_short_value " + str(args.tooShort) + " --do_not_require_softclips_both_ends -rf MappingQuality -mmq " + str(args.mapQuality) + " -rf ReadLength -minRead " + str(args.minRead) + " -maxRead " + str(args.maxRead) + " -o " + name + ".bam"
subprocess.call(command,shell=True)
# Calculate coverage
# command = args.bedtools + " genomecov -d -ibam " + name + ".bam > coverage.txt"
# subprocess.call(command,shell=True)
subprocess.call("rm tmp*[bs]a[im]",shell=True)
# Generate Bed file
for seq_record in SeqIO.parse(args.input.name, "fasta"):
	count = 0
	for i in range(1,len(seq_record.seq)): 			
		if count < len(seq_record.seq):
			n= open("Test.bed",'a+')
			n.write(seq_record.id+'\t'+str(i)+'\t'+str(i+args.window)+'\n')
			n.close()
			count = i+args.window
# Calculate coverage
command = args.bedtools + ' coverage -a Test.bed -b '+name + '.bam -f 1 > coverage.txt'
subprocess.call(command,shell=True)
# 

# Read in the coverage data
print("*"*int(columns))
print("Importing coverage information.")
X = pd.read_csv("coverage.txt",sep='\t',names=['transcript','pos1','pos2','cov','basecov','window','percent'])
print("Identified coverage data for " + str(len(set(X['transcript']))) + " transcripts.")
T = list(set(X['transcript']))

reportFile = open("report.txt","w")
goodFasta = open("good.fasta","w")
badcovFasta = open("badcov.fasta","w")
zeroBad = []
# Check for any sites with no coverage
for i in T :
    zeros = list(X[X['transcript'] == i]['cov']).count(0)
    if zeros :
        zeroBad.append(i)
        reportFile.write("BADZERO: " + i + " has " + str(zeros) + " sites with no coverage.\n")
reportFile.write("\n\n")
# Remove the contigs with zero coverage sites from the master list
T = [x for x in T if x not in zeroBad]
print("Removed " + str(len(zeroBad)) + " contigs for sites with zero coverage.")

# Check for bad coverage
covBad = []
for i in T :
    x = float(max(list(X[X['transcript'] == i]['cov'])))/float(min(list(X[X['transcript'] == i]['cov'])))
    if x > args.foldCoverage :
        covBad.append(i)
        reportFile.write("BADCOV: " + i + " has a " + str(round(x)) + "-fold coverage differential.\n")
    else :
        reportFile.write("GOOD: " + i + " has a " + str(round(x)) + "-fold coverage differential.\n")
print("Identified " + str(len(covBad)) + " contigs with >" + str(round(args.foldCoverage)) + "-fold coverage differentials.")

T = [x for x in T if x not in covBad]

# Function to generate mers for a single sequence (character string)
def MerSplit(merSize, seq) :
    mers = []
    if len(seq) < merSize : 
        print("Sequence length (" + str(len(seq)) + ") is less than mer size (" + str(merSize) + ").")
        sys.exit()
    for i in range(0,len(seq)-merSize+1) :
        mers.append(seq[i:i+merSize])
    return mers

# Function to identify exactly matching regions between two sequences
# The return value will be a list of three-element lists, where the first element is a list of the starting and ending positions
# in seqence 1, the second is a list for sequence 2, and the third is the actual sequence match.
def ExactMatches(merSize,seq1,seq2) :
    if len(seq1) < merSize or len(seq2) < merSize :
        print("Sequence length is less than mer size (" + str(merSize) + ").")
        sys.exit()
    matches = []
    seq1mers = MerSplit(merSize,seq1)
    seq2mers = MerSplit(merSize,seq2)
    i = 0
    while i < len(seq1mers) :
        if seq1mers[i] in seq2mers : 
            match1 = [i,i+merSize]                # get seq1 coordinates for initial match
            seq = seq1mers[i]                     # get initial match sequence
            j = seq2mers.index(seq1mers[i])       # get seq2 coordinates for initial match
            match2 = [j,j+merSize]
            i += 1
            j += 1
            while i < len(seq1mers) and j < len(seq2mers) and seq1mers[i] == seq2mers[j] :    # extend match while 2 sequences agree
                match1[1] += 1
                match2[1] += 1
                seq += seq1mers[i][-1]
                i += 1
                j += 1
            match1[1] += -1                       # to correct for overstepping by one 
            match2[1] += -1
            matches.append([match1,match2,seq])
        else :
            i += 1
            continue
    return matches

# Create a dictionary of all of the contigs.
seqFile = open(args.input.name,"r")
S = {}
for seq in SeqIO.parse(seqFile,"fasta") :
    #create fasta file for good and bad coverage sequences
    if seq.id in T:
    	goodFasta.write('>'+str(seq.id)+'\n'+str(seq.seq)+'\n')
    if seq.id in covBad:
    	badcovFasta.write('>'+str(seq.id)+'\n'+str(seq.seq)+'\n')
    S[seq.description] = str(seq.seq)
seqFile.close()

reportFile.write("\n\n")

# Generate coverage plots and exact matches for bad-coverage contigs
print("Generating coverage plots for " + str(len(covBad)) + " contigs with >" + str(round(args.foldCoverage)) + "-fold coverage differentials.")
for i in covBad :
    matches = {}
    for j in T :
        match = ExactMatches(args.merSize,S[i],S[j])
        if match :
            matches[j] = match
    if len(matches.keys()) == 0 :
        reportFile.write("No GOOD matches identified for " + i +".\n")
        x = list(X[X['transcript'] == i]['pos1'])
        y = list(X[X['transcript'] == i]['cov'])
        plt.figure(figsize=(12,4))
        plt.title(i)
        plt.ylabel('Coverage')
        plt.xlabel('CDS position')
        plt.plot(x,y,'b-',linewidth=2)
        plt.xlim(x[0],x[-1])
        plt.savefig('badcov_'+i+".pdf")            
        plt.close()
        continue

    x = list(X[X['transcript'] == i]['pos1'])
    y = list(X[X['transcript'] == i]['cov'])

    plt.figure(figsize=(12,4+4))
    plt.subplot(2,1,1)
    plt.title(i)
    plt.ylabel('Coverage')
    plt.xlabel('CDS position')
    plt.plot(x,y,'b-',linewidth=2)
    plt.xlim(x[0],x[-1])

    if not len(matches.keys()) :
        continue
    else :
        plt.subplot(2,1,2)
        plt.xlim(x[0],x[-1])
        plt.xlabel('CDS position')
        plt.ylim(0,len(matches.keys())+1)
        pos = 1
        for m in matches.keys() :    
            plt.text(0.1*x[-1],pos+0.25,m)
            for n in range(0,len(matches[m])) :
                a = [matches[m][n][0][0],matches[m][n][0][1]]
                plt.hlines(pos,a[0],a[1],"r",linewidth=3)
                plt.vlines(a[0],pos+0.1,pos-0.1)
                plt.vlines(a[1],pos+0.1,pos-0.1)
            pos += 1
    plt.savefig('badcov_'+i+".pdf")            
    plt.close()
if args.pdfprint == 'Y':
	for i in T:
		x = list(X[X['transcript'] == i]['pos1'])
		y = list(X[X['transcript'] == i]['cov'])
		plt.figure(figsize=(12,4))
		plt.title(i)
		plt.ylabel('Coverage')
		plt.xlabel('CDS position')
		plt.plot(x,y,'b-',linewidth=2)
		plt.xlim(x[0],x[-1])
		plt.savefig('good_'+i+".pdf")            
		plt.close()
reportFile.close()
goodFasta.close()
badcovFasta.close()
print("*"*int(columns))
