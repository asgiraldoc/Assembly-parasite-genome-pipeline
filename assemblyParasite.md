---
title: "Assembly parasite genome pipeline"
author: "Axl S. Cepeda"
date: "04/20/2022"
---


---
output:
  html_document: default
  pdf_document: default
---
# Polishing illumina reads

### Quality data
```
fastqc SAMPLE_F1.fq SAMPLE_R1.fq ... SAMPLE_F3.fq SAMPLE_R3.fq

```

### Remove bad quality reads depend on your previous results, eg. (do it for each pair-end library if necessary):

```
F1: 
trimmomatic SE -phred33 SAMPLE_F1.fastq F1_clean.fastq HEADCROP:41 SLIDINGWINDOW:4:15 MINLEN:31

R1:
trimmomatic SE -phred33 SAMPLE_R1.fastq R1_clean.fastq CROP:259 SLIDINGWINDOW:4:15 MINLEN:31


```


# Remove host illumina reads (do it for each pair-end library) ###
see this post for some explanations https://www.metagenomics.wiki/tools/short-read/remove-host-sequences
```
bowtie2-build bird_genome.fna birdDB

bowtie2 -x birdDB -1 F1_clean.fastq -2 R1_clean.fastq --no-unal -p 12 -S unmappedReads.sam

samtools view -bS unmappedReads.sam > unmappedReads.bam

samtools sort -n unmappedReads.bam -o unmappedReads_sorted.bam

samtools fastq -@ 8 unmappedReads_sorted.bam -1 unmappedReads_sorted_F1.fastq.gz -2 unmappedReads_sorted_R1.fastq.gz -0 /dev/null -s /dev/null -n
```

# Remove host illumina PacBio ###
```
makeblastdb -in bird_genome.fna -parse_seqids -dbtype nucl -out birdDB

blastn –query PacBio_reads.fa –db birdDB –out hostMatches.csv –outfmt 6 -evalue 1e-10
```
## remove with a custom script the sequences with match in the blast output e.g.:
```
awk -F "," '{print $2}' hostMatches.csv > listhostMatches.txt

grep -v -f -A1 listhostMatches.txt PacBio_reads.fa > noHost_PacBio_reads.fa
```

# Hybrid Assembly (PacBio + illumina)


### create a config MaSurCa file (parameters should be modify if necessary):

see see MaSurCa manual https://bioinformaticsworkbook.org/dataAnalysis/GenomeAssembly/Assemblers/MaSuRCA.html#gsc.tab=0


Note: I highly recommend to generate different assemblies using different parameter combinations

```
DATA
PE = pa 150 15 /home/asgiraldoc/parasite_genome/libs/unmappedReads_sorted_F1.fastq.gz /home/asgiraldoc/parasite_genome/libs/unmappedReads_sorted_R1.fastq.gz
PE = pb 100 15 /home/asgiraldoc/parasite_genome/libs/unmappedReads_sorted_F2.fastq /home/asgiraldoc/parasite_genome/libs/unmappedReads_sorted_R2.fastq
PE = pc 100 15 /home/asgiraldoc/parasite_genome/libs/unmappedReads_sorted_F3.fastq /home/asgiraldoc/parasite_genome/libs/L008_un_R.fastq
PE = pd 251 15 /home/asgiraldoc/parasite_genome/libs/DRR091073_un_F.fastq /home/asgiraldoc/parasite_genome/libs/unmappedReads_sorted_R2.fastq
PACBIO=/home/asgiraldoc/parasite_genome/libs/noHost_PacBio_reads.fasta
END

PARAMETERS
GRAPH_KMER_SIZE = 99
USE_LINKING_MATES = 0
CLOSE_GAPS=1
NUM_THREADS = 22 
JF_SIZE = 60000000
SOAP_ASSEMBLY=0
CA_PARAMETERS =  cgwErrorRate=0.15
CLOSE_GAPS=1
USE_GRID=0
END
```

### run masurca 
```
masurca masurca_config.txt ### this generate a custom script call assemble.sh

./assemble.sh
```

# QUAST: quality assessment tool for genome assemblies 
see Quast manual to read quality assembly metrics https://github.com/ablab/quast


```
python quast.py assemby/assembly_test1.fasta assemby/assembly_test2.fasta -r assemby/referenceParasite.fasta -o quastComparison_output

```
