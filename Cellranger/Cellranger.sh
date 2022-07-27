#!/bin/bash

#source /etc/profile.d/set.sh

if [ $# -ne 4 ]
then
    echo "Usage: ./cellranger.sh [id] [fastqpath] [sampleID] [t] "
    exit 65
fi

id=$1
fastqpath=$2
sampleID=$3
t=$4

cellranger count --id=$id \
		 --fastqs=$fastqpath \
		 --sample=$sampleID \
		 --localcores=$t \
              	 --transcriptome=/share/swap/refdata-cellranger-hg19-1.2.0/  
		 #--transcriptome=/export/data0/refdata-cellranger-hg19-3.0.0/ 

