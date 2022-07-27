### File Format : EMs1_S1_L001_R1_001.fastq.gz,EMs1_S1_L001_R2_001.fastq.gz;Eut1_S1_L001_R1_001.fastq.gz,Eut2_S1_L001_R1_001.fastq.gz
### According to the dataset name
for i in `echo EMs1 Eut1`;do
     bash Cellranger.sh $i /share/data5/SingleCellSeq/op-scRNA/ $i 40 
done

