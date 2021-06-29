# Transcriptomic project
## SRA selection
The species in the *Strongyloides* genus are soil-transmitted gastrointestinal parasites of human and other animals. The female parasites produce genetically indentical offspring by mitotic parthenogenesis. Eggs leave the host trough faeces and develop larvare (iL3). These are able to infect a new host as well as develop into a dioecious free-living adult stage. Therefore, the *Strongyloides* genus has the almost unique feature to produce genetically identical parasite adults and free-living adults. In these analysis we retrive the SRA files from 3 samples of free-living females and 3 samples of parasitic females, in order to test the differences in expression, since the two stages are genetically identical we hypotise that the differences between the two ecological life-style are due to a different gene transcription. 

The SRA code detected for this analysis are:
| SRA code | Ecological feature | Sample |
| -------- | ------------------ | ------ |
| DRR106347 | free living | free_s1 |
| DRR106349 | free living | free_s2 |
| DRR106350 | free living | free_s3 |
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

The de novo assembly allow us to obtain an assembled transcriptome from our raw reads. The assemble will be then used as a refernce to map the reads. Since we are dealing with a single spiecies, moreover the sample are genetically identical, all the trimmed files will be merged in a file that contains all the left reads and a file with all the right reads. The Trinity algorithm performs better when the size of both files is around 10G, the bash scrpit *random_subsampling_PE.sh* samples randomly a number of reads. We used the script to reduce the number of reads and the size of the files. A bash script was written to cycle the porcess for every sample
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

First the reference, which is the transcriptome, is indexed. Then, the pair ends of each sample are mapped on the indexed reference trough `<bowtie2>`. The comand launched for the first sample is as follows:
```
bowtie2-build cdhit_ouput.fasta indexed_ref
bowtie2 -x ../../cdhit/indexed_ref/references -1 sp1_pr1 -2 sp1_pr2 -S mapped_sp1.sam    > alignment_rate
```
The results of the mapping for each sample

| Sample | Al. conc. | Al. discord. | overall all. |
| ------ | --------- | ------------ | ------------ |
| free_s1 | 68% | 13% | 90% |
| free_s2 | 73% | 12% | 90% |
| free_s3 | 70% | 15% | 92% |
| para_s1 | 61% | 14% | 90% |
| para_s2 | 59% | 15% | 90% |
| para_s3 | 60% | 13% | 90% |


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

#We edit the complete dataframe, in order to have the tanscript names as row names. 
stot=data.frame(stot[,2:7],row.names=stot[,1])


#We create a dataframe that indicates at which group each sample (or column) belongs, as requested from NOISeq2
myfactors=data.frame(LifeStyle=c("Free","Free","Free","Para","Para","Para"))

#Now we pack all the information in a NOISeq object
mydata=readData(data=stot, factors=myfactors)

#We create the saturation plot
mysaturation=dat(mydata,k=0, ndepth=7,type="saturation")
explo.plot(mysaturation, toplot = 1, samples = 1:6, yleftlim = NULL, yrightlim = NULL)

#and the sensitivity plot
mycountsbio = dat(mydata, factor = NULL, type = "countsbio")
explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "barplot") #acqtually, I did not understand what the horizontal lines mean.

#Filtering the loci with low counts, 4767 out of 31727 features have been selected. 
myfilt10 = filtered.data(stot, factor = myfactors$LifeStyle, norm = FALSE, depth = NULL, method = 1, cv.cutoff = 100, cpm = 10, p.adj = "fdr")

#Normalization
mydata_TMM10 = tmm(assayData(mydata)$exprs, long = 1000, lc = 0)

#Differential analysis
mydata = readData(data = mydata_TMM10, factors=myfactors)
mynoiseqbio_para_free_t0=noiseqbio(mydata, k=0.1, norm="n", filter=0, factor="LifeStyle")


#With the function degenes we can select the features that have a high probability to be differentially expressed, firstly we select all the significantly differentially expressed features, then we select the ones more expressed in free saples, and finally tyhe ones more expressed in the parasitic samples.
mynoiseqbio_para_free_t0_degtot=degenes(mynoiseqbio_para_free_t0, q= 0.95, M = NULL) #3444 featues
mynoiseqbio_para_free_t0_degfree=degenes(mynoiseqbio_para_free_t0, q= 0.95, M = "up") #1864 features
mynoiseqbio_para_free_t0_degpara=degenes(mynoiseqbio_para_free_t0, q= 0.95, M = "down") #1580 features

#We plot resplectively the expression plot and the MD plot (D: the absolute value of the difference in expression betweeen the two conditions; M: log-fold change). 
DE.plot(mynoiseqbio_para_free_t0,q = 0.95, graphic = "expr", log.scale = TRUE)
DE.plot(mynoiseqbio_para_free_t0,q = 0.95, graphic = "MD")
```

## GO enrichment

From the output of Panzzer2 we create a file, where for each transcript are listed the GO terms assigned to the genes located on that transcript. I realised that most of the transcripts harbor more than one ORF, therefore in these cases the GO terms of one transcript refer to different genes localized on the same transcript
```
for b in `awk '{print$1}' GO_filtered.out | sort | uniq`; do stringa=""; for a in `grep -w $b GO_filtered.out | awk '{print$2}'`; do stringa="${stringa}GO:${a}, "; done; echo -e "$b\t${stringa%, }"; done > GO_per_transc
```
Then we use the topGO R package to perform the GO enrichment
```
geneID2GO=readMappings(file="all_genes_annotation_ok")
geneNames = names(geneID2GO)

#First we load the genes that are significantively differentially expressed
tab=read.table("diff_expr_para-free_tot.txt")
genes_int_list=as.vector(rownames(tab))

#We check how many transcripts do not have a any GO term, 963 on 3444
length(genes_int_list[!(genes_int_list %in% geneNames)])

#We create the vector geneList, topGO requires a vector of gene with a p-value referred to the probability to be differentially transcribed (1 is significative), so we create a vector of 0 and 1, 1 for the transcript that we have already seen are diff. transc. (the ones in genes_int_list), the 0 are the transcript that are not in genes_int_list, thus not diff. transc.. Each element of the vector will be named with the transcript name
geneList <- factor(as.integer(geneNames %in% genes_int_list))
names(geneList) = geneNames

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GOdata, classicFisher = resultFis,ranksOf = "classicFisher", topNodes = 447) #477 is the number of significative terms. 
#I did not quite understand how they calculate the expected number for each GO term
```
Then, the GO enrichment has been performed for the other two categories (CC and MF). Moreover, we analysed the GO terms of the transcripts uptranscribed in the parasitic condition and in the free living codition (the pipeline in R is the same as above apart from: `tab=read.table("diff_expr_para-free_para.txt") and tab=read.table("diff_expr_para-free_free.txt")`. The results has been loaded in the folder ....

!!!! IDEA PER LA CONCLUSIONE: i para hanno in generale una minore trascrizione e nei GO terms i primi riguardano RNA processing and ncRNA metabolic process, una hp potrebbe essere che durante il parassitismo c'è una attivazione di ncRNA che vanno a inattivare la trascrizione di molti geni, andando a ridurre in generale la trascrizione. Infatti i GO terms per i campioni a vita libera sono più collegati alla organizzazione e allo sviluppo cellulare.
## Annotation

First we will annotate the nucleotide sequences of the transcripts on the uniprot database. The database has been already built with `make db`. The output is in the TSV format, each column respectively rapresents: 1) Query id 2)Name of the target 3) evalue 4)bitscore 5)percentage of identical matches 6)description of the target
```
diamond blastx --db /var/local/uniprot/uniprot_sprot.fasta-norep_per_diamond.dmnd --query ../cdhit/cdhit_ouput.fasta -p 16 -o output --outfmt 6 qseqid sseqid evalue bitscore pident stitle --max-target-seqs 5 --evalue 0.005
```

`TransDecoder` find possible coding regions inside the transcripts. The first step is detecting all the possible ORFs longer than 100 amino acids, the ouput file will be a fasta with the ORFs translated to amino acids sequences. Then the ORFs are annotated on a database, in order to detect a homolgy with proteins inside the db. The last step detects the most likely ORFs.
```
TransDecoder.LongOrfs -t ../cdhit/cdhit_ouput.fasta
diamond blastp --query cdhit_ouput.fasta.transdecoder_dir/longest_orfs.pep --db /var/local/uniprot/uniprot_sprot.fasta-norep_per_diamond.dmnd --evalue 1e-05 --max-targetseqs 1 --threads 5 --outfmt 6 --out blastp.outfmt6
TransDecoder.Predict -t ../cdhit/cdhit_ouput.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.ou
tfmt6

#To shorten the headers
cat cdhit_ouput.fasta.transdecoder.pep | sed 's/\(^>.\+\)\.p.*$/\1/g'  |  sed 's/\(^>.\+\)\.p.*$/\1/g' > transdecoder_final_out.fasta
```
The last file is the input for `panzzer2`, which requires a fasta file of amino acid sequences. Trough Panzzer2, at each ORF will be assigned the GO terms.   
