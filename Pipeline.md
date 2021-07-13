# Transcriptomic project
## SRA selection
The species in the *Strongyloides* genus are soil-transmitted gastrointestinal parasites of human and other animals. The female parasites produce genetically indentical offspring by mitotic parthenogenesis. Eggs leave the host trough feces and develop larvae (iL3). These are able to infect a new host as well as develop into a dioecious free-living adult stage. Therefore, the *Strongyloides* genus has the almost unique feature to produce genetically identical parasite adults and free-living adults (Hunt et al. 2018). In these analysis we retrieve the SRA files from 3 samples of free-living females and 3 samples of parasitic females, in order to test the differences in expression, since the two stages are genetically identical we suggest that the differences between the two ecological life-style are due to a different gene expression, which can be detected in the transcriptomes. 

The SRA code detected for this analysis are:
| SRA code | Ecological feature | Sample | N of raw reads (1) |
| -------- | ------------------ | ------ | --------------- |
| DRR106347 | free living | free_s1 | 16.6 M |
| DRR106349 | free living | free_s2 | 17.7 M |
| DRR106350 | free living | free_s3 | 11.9 M |
| DRR106353 | parasitic | para_s1 | 16.7 M |
| DRR106354 | parasitic | para_s2 | 13.6 M |
| DRR106356 | parasitic | para_s3 | 14.3 M |

(1) Being more precise, it is the number of paired raw reads after the trimming step

## SRA download and quality evaluation
First, with fastqc we download the SRA, which have been previously selected. The reads are pair ends, therefore we will specify it with the flag `<--split-files>`, this command divides the right forward and the reverse reads in two different fastq files. The command:
```
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRACODE
```

Then, to make sure the downloads have been successful, we validate the files in the NCBI folder
```
vdb-validate *.sra
```

To check the reads quality we use `fastqc` for each pair of fastq files: 
```
fastqc *_1.fastq *_2.fastq -o fastqc
```

According to the .html files, you can find them in this [folder](https://github.com/AlessandroFormaggioni/Transcriptomic_Uni_project/tree/main/fastqc), the quality of the raw reads is good, we can safely say that since:
1. The median phred score per base is always under the threshold of 28, even for the last bases (which are known to be the ones called with less reliability)
2.  The frequency of the nucleotides is steady. Despite for the beginning of the read, which is likely due to the adaptors or, more in general, to the building process of the library
3.  The GC content of the reads follows a normal distribution with one peak, meaning a lack of contaminations.
On the other hand, we have to mention that some sequences are duplicated. This could be a computational problem (since during the PCR some sequences were more amplified than others) or a biological one (these sequences are just more expressed).

## Trimming

Through the trimming step some parts of the reads are deleted, in particular: 
* The ends of the reads with low quality
* The adapters, indeed we provide a list of ILLUMINA adapters to allow `trimmomatic` to recognize them in the sequence

In the following command we execute the trimming for 3 samples in a row: 
```
for a in 1 2 3; do
trimmomatic PE -threads 5 -phred33 $a/*_1.fastq $a/*_2.fastq $a/"sp"$a"_pr1" $a/"sp"$a"_unpr1" $a/"sp"$a"_pr2" $a/"sp"$a"_unpr2" ILLUMINACLIP:/usr/local/anaconda3/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75 2> $a/stats_trimmomatic 
done
```
After the trimming, the reads are divided in 2 groups, the ones that are paired and the ones that are not paired, both groups are again divided in forward and reverse reads, for a total of 4 fastqc files. In this analysis we are going to use just the paired reads.
According to the statics related to the trimming, most of the reads were paired, between 88.8% and 90%. The unpaired forward ends were between 6.6%  and 3.9%, whereas the unpaired reverse ends were 2%-1%. Discarded reads were around 1%.  

## Assembly and evaluation of the transcriptome

The de novo assembly allows to obtain an assembled transcriptome from the raw reads. Then, the transcriptome is used as a reference to map the reads. Since we are dealing with a single species and the samples are genetically identical, all the trimmed files will be merged in a single file that contains all the left reads and a file with all the right reads, in order to obtain a single transcriptome from all the samples. The Trinity algorithm performs better when the size of both files is around 10G. The bash script *random_subsampling_PE.sh* samples randomly a number of reads. We used the script to reduce the number of reads and the size of the files, reducing the number of reads to one third. A bash script was written to cycle the process for every sample
```
#!/bin/bash

cd Para
for a in 1 2 3; do
cd $a
cp ../../random_subsampling_PE.sh .
conto=$(wc -l *_pr1 | awk '{print$1}')
#Divided by 12 because I want to reduce the reads to 1/3 and the number of reads is = the lenght of the file / 4
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

Now, with `cat` we merge all the subsamples in two files
```
cat Para/1/*1_subsampl.fastq Para/2/*1_subsampl.fastq Para/3/*1_subsampl.fastq Free/1/*1_subsampl.fastq Free/2/*1_subsampl.fastq Free/3/*1_subsampl.fastq > tot_pr1
cat Para/1/*2_subsampl.fastq Para/2/*2_subsampl.fastq Para/3/*2_subsampl.fastq Free/1/*2_subsampl.fastq Free/2/*2_subsampl.fastq Free/3/*2_subsampl.fastq > tot_pr2
```

The two file are used to assemble the transcriptome with Trinity:
```
Trinity --seqType fq --left tot_pr1 --right tot_pr2 --CPU 6 --max_memory 20G 

```
Different software and statics allow the evaluation of the transcriptome:
* BUSCO is a software that aligns the transcripts to a databases of genes, the aim is to check if in the transcriptome there are a set of core genes that are usually present in all transcriptomes. The percentage of core genes detected in the transcriptome can be considered a percentage of completeness of the transcriptome. We performed different analysis on BUSCO, in order to see whether a different pipeline or database could affect the result:

| Pipeline | Database | Completeness |
| -------- | -------- | ----------- |
| BUSCO v5 | Nematoda | 33.2% |
| BUSCO v5 | Metazoa | 72% |
| BUSCO v3/v2 | Nematoda | 45.2% |

As you can see, the scores are highly dependent on the database selection. In our case, selecting the Metazoa database leads to an higher BUSCO score. This could be due to a compositional bias of the database: the Nematoda database is mainly built on *C. elegans* core genes, which could differ from the *Strongyloides* core genes. The latter could be more similar to a more general Metazoa core genes database.  

* N50 and L50 are measures that define the transcriptome quality in terms of contigs length. Sorting the contigs from the longest to the shortest, we define the group of longest contigs that cover 50% of the total transcriptome length, the N50 is the length of the last contig of that group, whereas L50 is the number of contigs present in that group. <br /><br> **N50=925** <br /> <br>**L50=13765**<br /> The N50 measure can be retrived by the BUSCO output as well as calculated with the command ` TrinityStats.pl Trinity.fasta > Trinity_stats.txt`. However, I realized that the 2 methods calculate two different values for the N50, the one reported above are from BUSCO, while for TrinityStats N50=2331 (considering all contigs)

## Isoforms redundancy

Each gene can be present in the transcriptome with several isoforms. This could be a problem during the mapping, as the reads could map on all the isoforms, mapping on different locations. `cd-hit` clusters the transcripts with 90% of similarity and takes just the most repeated one. Thus we will have less transcripts but the pipeline will perform better.
```
cd-hit-est -i Trinity.fasta -o output.fasta -T 12 -t 1 –c 0.9
```


## Mapping

First the reference, which is the transcriptome, is indexed. Then, the pair ends of each sample are mapped on the indexed reference trough `bowtie2`. The command launched for the first sample is as follows:
```
bowtie2-build cdhit_ouput.fasta indexed_ref
bowtie2 -x ../../cdhit/indexed_ref/references -1 sp1_pr1 -2 sp1_pr2 -S mapped_sp1.sam    > alignment_rate
```
The results of the mapping for each sample. The pair ends reads align discordantly when, during the mapping, they do not maintain the distance that has been set during the sequencing process.

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
for a in 1 2 3; do cd $a; samtools view -h -Sb mapped_sp"$a".sam > mapped_sp"$a".bam; samtools view -h -f 0x2 -F 256 -q 30 -Sb mapped_sp"$a".sam > mapped_sp"$a"_filtered.bam; cd ..; done
```
Then we sort and index the bam file, to finally obtain the raw counts for each transcript
```
for a in 1 2 3; do cd $a; samtools sort mapped_sp"$a"_filtered.bam > mapped_sp"$a"_sortfilt.bam; samtools index mapped_sp"$a"_sortfilt.bam; samtools idxstats mapped_sp"$a"_sortfilt.bam >  sp"$a"_rawcounts.txt; cd ..; done
```
## Annotation

First we annotate the nucleotide sequences of the transcripts on the uniprot database. The database has been already built with `make db`. The output is in the TSV format, each column respectively represents: 1) Query id 2)Name of the target 3) Evalue 4) Bit-score 5) Percentage of identical matches 6) Description of the target
```
diamond blastx --db /var/local/uniprot/uniprot_sprot.fasta-norep_per_diamond.dmnd --query ../cdhit/cdhit_ouput.fasta -p 16 -o output --outfmt 6 qseqid sseqid evalue bitscore pident stitle --max-target-seqs 5 --evalue 0.005
```

`TransDecoder` find possible coding regions inside the transcripts. The first step is detecting all the possible ORFs longer than 100 amino acids, the ouput file will be a fasta with the ORFs translated to amino acids sequences. Then, the ORFs are annotated on a database (trough Diamond and HMMER), in order to detect a homology with proteins inside the db. The last step detects the most likely ORFs.
```
TransDecoder.LongOrfs -t ../cdhit/cdhit_ouput.fasta
diamond blastp --query cdhit_ouput.fasta.transdecoder_dir/longest_orfs.pep --db /var/local/uniprot/uniprot_sprot.fasta-norep_per_diamond.dmnd --evalue 1e-05 --max-targetseqs 1 --threads 5 --outfmt 6 --out blastp.outfmt6
hmmscan --cpu 8 --domtblout pfam.domtblout /var/local/Pfam/Pfam-A.hmm cdhit_ouput.fasta.transdecoder_dir/longest_orfs.pep
TransDecoder.Predict -t ../cdhit/cdhit_ouput.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6

#To shorten the headers
cat cdhit_ouput.fasta.transdecoder.pep | sed 's/\(^>.\+\)\.p.*$/\1/g'  |  sed 's/\(^>.\+\)\.p.*$/\1/g' > transdecoder_final_out.fasta
```
Once we have the ORFs, we can annotate them with HMMER and Diamond 
```
hmmscan --cpu 8 --domtblout hmmer_finalout.domtblout /var/local/Pfam/Pfam-A.hmm transdecoder_final_out.fasta
diamond blastp --query transdecoder_final_out.fasta --db /var/local/uniprot/uniprot_sprot.fasta-norep_per_diamond.dmnd --evalue 1e-05 --max-target-seqs 1 --threads 5 --outfmt 6 --out blastp_finalout.outfmt6
```
The output file of the annotations are in the folder [Annotations](https://github.com/AlessandroFormaggioni/Transcriptomic_Uni_project/tree/main/Annotations)


### KEGG pathways
Submitting the amino acid sequences to **KAAS** (KEGG Automatic Annotation Server), each sequence will be aligned to an ortholog group in the KEGG database, in order to assign the functional classification (KEGG Orthology, KO), each ortholog belongs to one or more KEGG pathways. There are different aligning algorithms, we chose the most performing one (GHOSTZ). Moreover, to define the dataset we selected the representative dataset for Eukaryotes, manually adding all the Nematode species available. In the results we can see for each transcripts at which orthologs it has been assigned. Based on the ortholog/functional assignment, we can see which are the most frequent pathways. For instance, in our case one of the pathways with more orthologs is the [Pathways of neurodegeneration - multiple diseases (ko05022)](https://github.com/AlessandroFormaggioni/Transcriptomic_Uni_project/blob/main/KEGG_path/ko05022_neurodeg.png) with 179 hits, clicking on the code we can graphically see the pathway and also get an idea of where our orthologs are located (the boxes highlighted in green), in this case they are wide spread in the whole pathway. However, other times they are restricted to specific reactions: in the pathway ["2-Oxocarboxylic acid metabolism" (ko01210)](https://github.com/AlessandroFormaggioni/Transcriptomic_Uni_project/blob/main/KEGG_path/ko01210_oxocarb.png) the highlighted reactions are almost restricted to one area, the ones that lead to the transformation of Oxaloacetate into Glutamate. 
Output link: <br />
https://www.genome.jp/kaas-bin/kaas_main?mode=user&id=1625750446&key=KS4L3fEr

## GO terms
The Gene Ontology (GO) terms are properties related to the genes, the properties are ordered hierarchically and divided in 3 main groups: Cellular Component (CE), Molecular Function (MF) and Biological Process (BP). The GO terms are all identified by a code and the more specific terms are contained by the  more general ones. It is important to assign the GO terms to our transcripts to get an idea of which are the more present terms in our transcriptome. We used Pannzer2 for this purpose, submitting a fasta which contains the amino acid sequence of the transcripts, the algorithm will search for GO terms related domains. In our case the amino acid sequences (transdecoder_final_out.fasta) are the input for Panzzer2. Through Panzzer2, at each ORF will be assigned the GO terms.

## Differential expression

To analyze the differential expression between conditions we will use the R package `NOISeq2`. We first have to create a data frame where the rows are the transcripts and the columns the six samples, the intersection between rows and columns is the raw count (how many raw reads map on that transcript):
```
#First import raw counts output, we are interested just in the first and third column (respectively the name of the transcript and the raw count).
sf1=read.table("sf1_rawcounts.txt",col.names=c("transc","","sf1",""))[,c(1,3)] #Here we report the upload of the first output, however all the output were uploaded in six different tables (sf1, sf2, sf3, sp1, sp2, sp3)

#With the function full_join from the dplyr package we merge all the tables according to the "transc" column (the column that contains the transcripts names). Full_join is able to merge two tables at one time, therefore we use the pipe symbol (%>%) to direct the output of the prevoius function as the input of the following function. 
library(dplyr)
stot=full_join(sf1,sf2,by="transc") %>% full_join(.,sf3, by="transc") %>% full_join(.,sp1, by="transc") %>% full_join(.,sp2, by="transc") %>% full_join(.,sp3, by="transc")

#We edit the complete dataframe, in order to have the transcript names as row names. 
stot=data.frame(stot[,2:7],row.names=stot[,1])

#The dataframe was saved and uploaded on GH as "NOISEQ_input_table.txt"

#We create a dataframe that indicates at which group each sample (or column) belongs, as requested from NOISeq2
myfactors=data.frame(LifeStyle=c("Free","Free","Free","Para","Para","Para"))

#Now we pack all the information in a NOISeq object
mydata=readData(data=stot, factors=myfactors)

#We create the saturation plot
mysaturation=dat(mydata,k=0, ndepth=7,type="saturation")
explo.plot(mysaturation, toplot = 1, samples = 1:6, yleftlim = NULL, yrightlim = NULL)

#and the sensitivity plot
mycountsbio = dat(mydata, factor = NULL, type = "countsbio")
explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "barplot") 

#Filtering the loci with low counts, 4767 out of 31727 features have been selected. 
myfilt10 = filtered.data(stot, factor = myfactors$LifeStyle, norm = FALSE, depth = NULL, method = 1, cv.cutoff = 100, cpm = 10, p.adj = "fdr")

#Normalization
mydata_TMM10 = tmm(assayData(mydata)$exprs, long = 1000, lc = 0)

#Differential analysis
mydata = readData(data = mydata_TMM10, factors=myfactors)
mynoiseqbio_para_free_t0=noiseqbio(mydata, k=0.1, norm="n", filter=0, factor="LifeStyle")


#With the function degenes we can select the features that have a high probability to be differentially expressed, firstly we select all the significantly differentially expressed features, then we select the ones more expressed in free samples, and finally the ones more expressed in the parasitic samples.
mynoiseqbio_para_free_t0_degtot=degenes(mynoiseqbio_para_free_t0, q= 0.95, M = NULL) #3444 featues
mynoiseqbio_para_free_t0_degfree=degenes(mynoiseqbio_para_free_t0, q= 0.95, M = "up") #1864 features
mynoiseqbio_para_free_t0_degpara=degenes(mynoiseqbio_para_free_t0, q= 0.95, M = "down") #1580 features

#We plot respectively the expression plot and the MD plot (D: the absolute value of the difference in expression between the two conditions; M: log-fold change). 
DE.plot(mynoiseqbio_para_free_t0,q = 0.95, graphic = "expr", log.scale = TRUE)
DE.plot(mynoiseqbio_para_free_t0,q = 0.95, graphic = "MD")
```
The [saturation plot](https://github.com/AlessandroFormaggioni/Transcriptomic_Uni_project/blob/main/DE_plots/Saturation_plot.pdf) can give us a lot of information. The depth is lower for the samples sf3, sp2, sp3. This agrees with the depth reported on the NCBI page of each SRA and also for the number of paired raw reads in the fastq of each sample. The plot shows that with the same depth we are able to retrieve more features from the free-living samples than in the parasitic samples. However, in parasitic samples the percentage of detected features is lower, suggesting that more features have not been discovered. Moreover, both groups show a low number of features detected, suggesting that with an higher depth the number of features would increase significatively, indeed the curves do not reach the saturation point, especially the ones with a lower depth (sf3, sp2, sp3), which steadily increase, whereas the others show a lower rate of increment after a certain depth. It is also worth mentioning that there is not a clear relationship between depth and features detected, this is evident in the free-living samples where the sf2 is the one with the higher depth as well as the one with lower features detected among free-living samples. In my opinion this could be due to: a highly variability in the transcriptomes of the free living samples (the sf2 has actually transcribed fewer transcripts) or it could be a problem in the transcriptome assembly. In the parasitic samples we can see how the relationship is more linear, where sp2 and sp3 have both the same depth and features detected, whereas sp1 has higher depth and features detected.
<br />
The [sensitivity plot](https://github.com/AlessandroFormaggioni/Transcriptomic_Uni_project/blob/main/DE_plots/sensitivity_plot.pdf) shows the percentage of features that are in a specific range of CPM (count per milion, it means how many raw reads map on that transcript). Before the filtering we see how most of the transcripts are between 0 and 1. The sensitivity plot is useful to get an idea about how to set the threshold for the filtering of the loci with low counts. The [sensitivity plot](https://github.com/AlessandroFormaggioni/Transcriptomic_Uni_project/blob/main/DE_plots/sensitivity_filtered.pdf) after the filtering let me a little bit surprised: since sf3 has a low depth but a high n of features detected (and sf2 is the opposite), I thought that sf3 would have had lots of features with lower CPM and sf2 would have had lower features but with an higher CPM. However, the sensitivity plot after the filtering shows that sf3 is the sample with the highest amount of features with CMP above 10M. In my opinion, this could be a clue that the few features detected from sf2 are a computational problem (and not a biological one) and most of the raw reads of that sample do not map at all (although the statics after bowtie do not indicate a lower percentage of alignment for the sf2 reads. Maybe I am missing something from this personal plot analysis). 
<br />
The results of the DE analysis can be graphically represented with an [expression plot](https://github.com/AlessandroFormaggioni/Transcriptomic_Uni_project/blob/main/DE_plots/expression_plot.pdf), the dots in the upper-left corner are the ones more expressed in the parasitic samples, and the dots in the bottom-right corner are the ones more expressed in the free-living samples. Moreover, we can also [plot](https://github.com/AlessandroFormaggioni/Transcriptomic_Uni_project/blob/main/DE_plots/MD_plot.pdf) the absolute value of the difference in expression between the two conditions (D) with the log2-ratio of the two conditions, positive M value represent transcript more expressed in the free living samples. 

## GO enrichment

In this step we want to see which are the GO terms more significantly present in the differentially expressed transcripts. Since we calculated which are: 1) In general the diff. expr. transcripts 2)The more expressed ones in the parasitic samples 3)The more expressed ones in the free-living samples, we will see how the GO terms change in the two conditions.  <br \>
From the output of Panzzer2 we create a file where for each transcript its GO terms. I realized that most of the transcripts harbor more than one ORF, therefore in these cases the GO terms of one transcript refer to different genes located on the same transcript
```
for b in `awk '{print$1}' GO_filtered.out | sort | uniq`; do stringa=""; for a in `grep -w $b GO_filtered.out | awk '{print$2}'`; do stringa="${stringa}GO:${a}, "; done; echo -e "$b\t${stringa%, }"; done > GO_per_transc
```
Then we use the topGO R package to perform the GO enrichment
```
geneID2GO=readMappings(file="all_genes_annotation_ok")
geneNames = names(geneID2GO)

#First we load the genes that are significatively differentially expressed
tab=read.table("diff_expr_para-free_tot.txt")
genes_int_list=as.vector(rownames(tab))

#We check how many transcripts do not have any GO term, 963 on 3444
length(genes_int_list[!(genes_int_list %in% geneNames)])

#We create the vector geneList, topGO requires a vector of gene with a p-value referred to the probability to be differentially transcribed (1 is significative), so we create a vector of 0 and 1, 1 for the transcript that we have already seen are diff. transc. (the ones in genes_int_list), the 0 are the transcript that are not in genes_int_list, thus not diff. transc. Each element of the vector will be named with the transcript name
geneList <- factor(as.integer(geneNames %in% genes_int_list))
names(geneList) = geneNames

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GOdata, classicFisher = resultFis,ranksOf = "classicFisher", topNodes = 447) #477 is the number of significative terms. 
#I did not quite understand how they calculate the expected number for each GO term
```
Then, the GO enrichment has been performed for the other two categories (CC and MF). Moreover, we analyzed the GO terms of the transcripts up-transcribed in the parasitic condition and in the free-living condition (the pipeline in R is the same as above apart from: `tab=read.table("diff_expr_para-free_para.txt") and tab=read.table("diff_expr_para-free_free.txt")`. The results has been uploaded in the folder [GOenrich](https://github.com/AlessandroFormaggioni/Transcriptomic_Uni_project/tree/main/GOenrich)

## Differences in KEGG pathways
To improve the analysis we performed the **KAAS** analysis for:
1. The transcripts more diff. expr. in the free-living samples
2. The transcripts more diff. expr in the parasitic samples

To obtain this, first we create a fasta file for each, which contains the amino acid sequence of only the diff. expr. transcripts for that condition.
```
for a in `awk '{print$1}' diff_expr_para-free_para.txt`; do grep "$a" -w -A1 SRA/Diamond/transdecoder_final_out_ol.fasta ; done > transecoder_DEpara.fasta
for a in `awk '{print$1}' diff_expr_para-free_free.txt`; do grep "$a" -w -A1 SRA/Diamond/transdecoder_final_out_ol.fasta ; done > transecoder_DEfree.fasta
```
Then, we submitted each fasta file to KAAS, with the same setting as before. The results were not particularly interesting (you can check them from the link above): the only notable pathway is "Splicesome", more present in the parasitic samples (29 hits vs 11). As we will see in the conclusions, the splicesome seems to have a big role in parasitism.

## Conclusions

<br> According to the GO terms, lots of transcripts overexpressed in the PSs (parasitic samples) are related to specific and linked biological processes: "ncRNA metabolic process", "ncRNA processing", "RNA processing", "gene expression". All these transcripts are likely to be involved in small non coding RNA maturation, leading to a different transcription and genetic regulation in the PSs. On the other side, the transcripts overexpressed in the FSs (free-living samples) are assigned to biological processes more related to growth and development: "anatomical structure development", "developmental process", "multicellular organism development", "organelle organization", "cellular component organization", "regulation of locomotion" (these are just some of the GO terms are the top of the FSs list). My personal hypothesis is that in PSs the maturation of sncRNA leads to the inactivation of some genes involved in the development of a free-living animal. This hypothesis is in line with the lower number of features detected in the PSs and with the shared idea that parasitism is a strategy that leads to save energy and cut unnecessary metabolic pathways. <br />
In the first part of the conclusions I tried to give a personal explanation of the results. In this part I try to compare my data to the ones obtained in the original paper (Hunt et al. 2018):
1. Authors claim that Argonaute-like proteins have a putative role in the parasitic life cycle, which agrees with the increased maturation of sncRNA in our data. 
2. There are several molecular functions that have a putative role in *Strongyloides* parasitism. The ones both reported in the paper and found in our analysis are: "acetylcholinesterase activity", "cholinesterase activity", "ubiquitin-like protein ligase binding". The latter is connected to the "Splicesome" pathway detected in the KEGG analysis 
3. Also in the paper they report an higher expression of key genes families in the FSs
4. In the paper 31 genes encoding speckle-type POZ protein-like (SPOP-like) proteins were upregulated in PSs. I searched in the HMMER output file all the target names equal to "Skp1_POZ" (although I am not entirely sure these proteins correspond to the speckle-type POZ protein-like), detecting 11 transcripts that were annotated to that protein. Of those 11 transcripts, only 1 is differentially expressed in the parasitic samples. 
<br />
Overall, I am satisfied by the analysis: the data allowed to guess a biological hypothesis and many results agree with the reference paper. However, it is still unclear how much robust our data are, since:
1. The saturation plot revealed that transcriptomes have a low depth.
2. The completeness of the assembly is highly dependent on which database is chosen.
3. FSs and PSs have different depths. Therefore we are not able to safely say that the lower number of features detected in PSs are due to a lower transcription level. 

## Reference

Hunt, V.L., Hino, A., Yoshida, A. et al. Comparative transcriptomics gives insights into the evolution of parasitism in Strongyloides nematodes at the genus, subclade and species level. Sci Rep 8, 5192 (2018). https://doi.org/10.1038/s41598-018-23514-z
