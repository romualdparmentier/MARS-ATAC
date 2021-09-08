#!/bin/bash

# esearch comes from NCBI Sequence Read Archive (SRA) Toolkit (see : https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit )
# Also see : https://www.reneshbedre.com/blog/ncbi_sra_toolkit.html for information and another way to do it

# GEO dataset page link : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156735
# GEO accession number : GSE156735
# BioProject ID : PRJNA658921

#########
## Download all paired-end sra files (i.e ID starting with SRR) and convert it into .fastq files (fasterq-dump -S)
#########

esearch -db sra -query 	PRJNA658921  | efetch --format runinfo | cut -d ',' -f 1 | grep SRR | xargs fasterq-dump -S

#########
# Run the following to extract only the ten first sport of the 5 first files of the project to test
#########
# esearch -db sra -query PRJNA658921  | efetch --format runinfo | cut -d ',' -f 1 | grep SRR | head -5 | xargs fastq-dump -X 10 --split-files


#########
## Download all supplement data to skip genome alignment (ATAC-seq series)
#########
wget -r -np -e robots=off 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156733/suppl/' | GREP GSE

#########
## Download all supplement data to skip genome alignment (RNA-seq series)
#########
wget -r -np -e robots=off 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156734/suppl/' | GREP GSE
