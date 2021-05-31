# Transcriptomic project

## SRA selection




## SRA download, quality evaluation and trimming

First, with fastqc we download the SRA, which have been previously selected. The reads are pair ends, therefore we will specift it with the flag `<--split-files>`, this comand
will divide the right reads and the left reads in two distinctive fastq. The command:
```
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRACODE
```

Then, to be sure the downloads have been successful, we validated the files in the NCBI folder
```
vdb-validate *.sra
```

To check the reads quality we used **fastqc** for each pair of fastq: 
```
fastqc *_1.fastq *_2.fastq -o fastqc
```
According to the .html files, the qulity of the raw reads is .... before the trimmming, as we can see from the figure below, the 
