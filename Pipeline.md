# Transcriptomic project
## SRA selection
The species in the *Strongyloides* genus are soil-transmitted gastrointestinal parasites of human and other animals. The female parasites produce genetically indentical offspring by mitotic parthenogenesis. Eggs leave the host trough faeces and develop larvare (iL3). These are able to infect a new host as well as develop into a dioecious free-living adult stage. Therefore, the *Strongyloides* genus has the almost unique feature to produce genetically identical parasite adults and free-living adults. In these analysis we retrive the SRA files from 3 samples of free-living females and 3 samples of parasitic females, in order to test the differences in expression, since the two stages are genetically identical we hypotise that the difference is in terms of transcription. 

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
