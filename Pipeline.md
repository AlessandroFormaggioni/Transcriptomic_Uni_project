# Transcriptomic project
## SRA selection
The species in the *Strongyloides* genus are soil-transmitted gastrointestinal parasites of human and other animals. The female parasites produce genetically indentical offspring by mitotic parthenogenesis. Eggs leave the host trough faeces and develop larvare (iL3). These are able to infect a new host as well as develop into a dioecious free-living adult stage. Therefore, the *Strongyloides* genus has the almost unique feature to produce genetically identical parasite adults and free-living adults. In these analysis we retrive the SRA files from 3 samples of free-living females and 3 samples of parasitic females, in order to test the differences in expression, since the two stages are genetically identical we hypotise that the differences between the two ecological life-style are due to a different gene transcription. 

The SRA code detected for this analysis are:
| SRA code | Ecological feature | Sample |
| -------- | ------------------ | ------ |
| DRR106347 | free living | free_s1 |
| DRR106349 | free living | free_s2 |
| DRR1063?? | free living | free_s3 |
| DRR106353 | parasitic | para_s1 |
| DRR106354 | parasitic | para_s2 |
| DRR106356 | parasitic | para_s3 |

## SRA download, quality evaluation and trimming

* First, with fastqc we download the SRA, which have been previously selected. The reads are pair ends, therefore we will specift it with the flag `<--split-files>`, this comand
will divide the right reads and the left reads in two distinctive fastq. The command:
```
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRACODE
```

* Then, to be sure the downloads have been successful, we validated the files in the NCBI folder
```
vdb-validate *.sra
```

* To check the reads quality we used **fastqc** for each pair of fastq: 
```
fastqc *_1.fastq *_2.fastq -o fastqc
```
* According to the .html files, the quality of the raw reads is .... before the trimmming, as we can see from the figure below, the 


**metti le immagini e ricordati di guardare le statistiche di trimmomatic**

* Trimming
According to the statics related to the trimming, most of the reads were paired, between 88.8% and 90%. The unpaired forward ends were between 6.6%  and 3.9, whereas the unpaired reverse ends were 2%-1%. Discarded reads were around 1%.  

## Assembly and evoluetion of the assemblage

The de novo assembly allow us to obtain a assembled transcriptome from our raw reads. The assemble will be then used as a refernce to map the reads. Since we are dealing with a single spiecies, moreover the sample are genetically identical, all the trimmed files will be merged in a file that contains all the left reads and a file with all the right reads. The Trinity algorithm performs better when the size of both files is around 10G, the bash scrpit *random_subsampling_PE.sh* samples randomly a number of reads. We used the script to reduce the number of reads and the size of the files. A bash script was written to cycle the porcess for every sample
```
#!/bin/bash

cd Para
for a in 1 2 3; do
cd $a
        cp ../../random_subsampling_PE.sh .
conto=$(wc -l *_pr1 | awk '{print$1}')
conto=$(( conto/12 ))
bash random_subsampling_PE.sh *_pr1 *_pr2 $conto
cd ..
done
cd ..
cd Free
for a in 1 2 3; do
cd $a
cp ../../random_subsampling_PE.sh .
conto=$(wc -l *_pr1 | awk '{print$1}')
conto=$(( conto/12 ))
        bash random_subsampling_PE.sh *_pr1 *_pr2 $conto  
        cd ..   
done 
```

Now, with *cat* we merge all the subsamples in two files
```
cat Para/1/*1_subsampl.fastq Para/2/*1_subsampl.fastq Para/3/*1_subsampl.fastq Free/1/*1_subsampl.fastq Free/2/*1_subsampl.fastq Free/3/*1_subsampl.fastq > tot_pr1
cat Para/1/*2_subsampl.fastq Para/2/*2_subsampl.fastq Para/3/*2_subsampl.fastq Free/1/*2_subsampl.fastq Free/2/*2_subsampl.fastq Free/3/*2_subsampl.fastq > tot_pr2
```

The two file are used to assemble the transcriptome with Trinity:
```
```
....
....


## Mapping

....index ..
The pair ends of each sample now will be mapped on the indexed reference trough `<bowtie2>`. We specify the left trimmed reads and the right trimmed reads, moreover with `<--no-discordant>` we specify that only pair ends that map concordantly will be included in the final SAM file. Two pair ends map concordantly when between them there is a fixed distance. The comand launched for the first sample is as follows:
```
bowtie2 -x ../../cdhit/indexed_ref/references -1 sp1_pr1 -2 sp1_pr2 -S mapped_sp1.sam  --no-discordant  > alignment_rate
```
The SAM file is heavy, therefore we translate the SAM file in BAM format, which is a less heavy binary file. We create an original copy of the mapping file and a filtered one, in the latter we keep only paired reads (actually all the reads are paired, since trimmomatic split the paired and unpaired reads in different files), with a mapping quality above 30, moreover we exclude all the secondary alignments. 
```
for a in 1 2 3; do cid $a; samtools view -h -Sb mapped_sp"$a".sam > mapped_sp"$a".bam; samtools view -h -f 0x2 -F 256 -q 30 -Sb mapped_sp"$a".sam > mapped_sp"$a"_filtered.bam; cd ..; done
```
Then we sort and index the bam file, to finally obtain the raw counts for each transcript
```
for a in 1 2 3; do cd $a; samtools sort mapped_sp"$a"_filtered.bam > mapped_sp"$a"_sortfilt.bam; samtools index mapped_sp"$a"_sortfilt.bam; samtools idxstats mapped_sp"$a"_sortfilt.bam >  sp"$a"_rawcounts.txt; cd ..; done
```

To analyse the differential expression between conditions we will use the R package ```NOISeq2```. We first have to create a data frame where the rows are the transcpts and the columns the six samples, the intesection between rows and columns is the raw counts of the reads of the sample that map on that specific transcript:
```
#First import raw counts output, we are interested just in the first and third column (respectively the name of the trasncript and he raw count).
sf1=read.table("sf1_rawcounts.txt",col.names=c("transc","","sf1",""))[,c(1,3)] #Here we report the upload of the first output, however all the output were uploaded in six different tables (sf1, sf2, sf3, sp1, sp2, sp3)

#With the function full_join from the dplyr package we merge all the tables according to the "transc" column (the column that contains the transcripts names). Full_join is able to merge two tables at one time, therefore we use the pipe simbol (%>%) to direct the output of the prevoius function as the input ofthe following function. 
library(dplyr)
stot=full_join(sf1,sf2,by="transc") %>% full_join(.,sf3, by="transc") %>% full_join(.,sp1, by="transc") %>% full_join(.,sp2, by="transc") %>% full_join(.,sp3, by="transc")

#We modify the complete dataframe, in order to have the tanscript names as row names. 
stot=data.frame(stot[,2:7],row.names=stot[,1])

#We create a dataframe that indicates at which group each sample (or column) belongs, as requested from NOISeq2
myfactors=data.frame(LifeStyle=c("Free","Free","Free","Para","Para","Para"))
