# Venom_Gland_Annotation
Custom Python scripts that are used to annotate assembled transcriptomic data.

-Assembled transcriptome needs to be searched against a protein database using ***BLASTX*** and save the .xml file

-Use the ***Blast_parse_matches.py*** script to narrow down the .fasta file to only those sequences that have a minimum percent of the total length of the protein identified from the .xml file

-Use the ***Blast_Namer.py*** script with the *matches* file and the .xml file to annotate the transcriptome

-The .fasta file that is output from ***Blast_Namer.py*** can then be used in the ***AutoNamer.py*** script to give names based on a .xml file. 

-The protein based annotations are done using the ***MS_Namer.py*** script that takes the ***Scaffold*** protein accession number report and annotates the genes from the proteins detected through mass spectrometery 

-The output file from ***MS_Namer*** can then be ran in ***AutoNamer.py*** script to give names based on a .xml file.

-Fully annotated transcripts can then be ran through the script ***ChimeraCheck.py*** to remove potential chimeric sequences. 

-The ***Busco.py*** script is used on a fasta file and the full_table.tsv file that is output from the program ***BUSCO***. The output from the ***Busco.py*** script needs to have a generated .xml file from the ancestral sequences given by ***BUSCO*** and then it can be used to name the sequences using ***Busco2.py***
