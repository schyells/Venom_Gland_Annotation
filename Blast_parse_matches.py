#!/usr/bin/env python

# The arguments to be passed at the command line are the contigs file, the blast results file (xml), and the percentage match desired. 
# You must have Biopython installed to run this.

import os, sys
from Bio.Blast import NCBIXML
from Bio import SeqIO

# The name of the file containing the sequences that were blasted.
# This should be in fasta format and provided as the first command-line argument.
inputfile=sys.argv[1]
if not os.path.isfile(inputfile) : 
    print("Fasta file does not exist: " + inputfile)
    exit(1)
else :
    inputFasta = open(inputfile,"r")

# The name of the file storing the blast results in xml format. This should be provided as the second command-line argument.
blastfile=sys.argv[2]
if not os.path.isfile(blastfile) : 
    print("Blast file does not exist: " + blastfile)
    exit(1)
else :
    blastRecords = open(blastfile,"r")

# This will limit the minimum percentage of the total length of a match for a sequence to be included.
# It is meant to exclude small fragements, exons, etc. It should be provided as the third command-line argument. 
minmatch=int(sys.argv[3])
if minmatch not in range(0,100) : 
    print("Unallowable minimum match percentage: " + str(minmatch))
    exit(1)

# How many characters to print for the blast match
titlelength=125

contigs=[]
all=[]

###################################################################################################################################
# First deal with the sequences that match the set of hard-coded keywords. Note that only contigs passing the minimum match
# percentage threshold get added to the "all" list, which is why it is not necessary to check this again for outputting the
# non-keyword-matching contigs. 
count=0
blastRecs = NCBIXML.parse(blastRecords)

# Output for the toxin sequences
outputFasta = open("matches_int.fasta","w")

# Find the contigs that have matches to the keywords
excluded=0
nomatch=0

for rec in blastRecs :
    count += 1
    if len(rec.alignments)==0 : 
        nomatch += 1
        continue
    if 100*max([float(abs(hsp.sbjct_end - hsp.sbjct_start))/al.length for al in rec.alignments for hsp in al.hsps]) <  minmatch : 
        excluded += 1
        continue
    all.append(str(rec.query.split()[0]))

total=count

# Eliminate duplicates
unique = list(set(all))

# Extract the keyword contigs from the fasta file and create a new fasta file
for rec in SeqIO.parse(inputFasta,"fasta") : 
    for contig in unique : 
        if contig == rec.description.split()[0] :
            SeqIO.write(rec,outputFasta,"fasta")

# make an uninterleaved version
output=open("matches.fasta","w")

outputFasta.close()

with open("matches_int.fasta","r") as f :
    count=0
    for line in f :
        if line[0] == ">" and count == 0 :
            count = 1
            output.write(line.rstrip() + "\n")
            continue
        if line[0] == ">" and count != 0 :
            output.write("\n" + line.rstrip() + "\n")
        else : 
            output.write(line.rstrip())

output.write("\n")

# Delete the interleaved versions and close the output file
os.remove("matches_int.fasta")
output.close()

blastRecords.close()
inputFasta.close()

print("\nExcluded {0} of {1} contigs from consideration for match percentages less than {2}% and {3} had no matches.".format(excluded,total,minmatch,nomatch))


