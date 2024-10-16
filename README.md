---
authors: "Lukas Dreyling & Henrique Valim"
date: "13-16.11.2023"
---

# 2023 Mycology Master's module (Biodiversity and Ecosystem Health) eDNA barcoding pipeline

The purpose of this pipeline is to process Illumina (MiSeq 2 x 300 bp) metabarcoding data of soil and bark samples taken from the SYMBIODRIVE project of the Biodiversity Exploratories in the frame of the Mycology Masters Course at the Goethe University Frankfurt. 

To that end, we will utilize some command line tools that are built for use in a Unix (i.e. Linux or macOS) environment, as well as make heavy use of R/RStudio for most of the analysis.

Below you can see the full list of programs/tools/packages that comprise our analytical pipeline and some brief descriptions of what those tools are used for.

## Program List

### Command line tools:
Our command line tools are installed using [conda](https://docs.conda.io/en/latest/), which is a package, dependency, and environment management tool for Unix systems (although it can be used on Windows, too). The two main advantages of conda is 1) it allows you to easily install new programs and tools, and 2) it keeps your tools in separate "environments" to avoid dependency-related problems between different versions of certain packages.

If you want to know more about what this means and how to create and use different environments on conda, [you can go here](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-environments).

 * [BLAST+](https://github.com/ncbi/blast_plus_docs): sequence alignment tool
 * [cutadapt](https://cutadapt.readthedocs.io/en/stable/): demultiplexing and primer removal
 * [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/): generating quality reports
 * [multiQC](https://multiqc.info/): aggregating multiple reports into a single report

### R-based tools:

Most of our analyses will be done in [R](https://www.r-project.org/) using the [RStudio](https://posit.co/products/open-source/rstudio/) graphical environment.

Specific tools in R come in "packages" that we need to install and load into our working environment. 

  * [dada2](https://benjjneb.github.io/dada2/): Accurate, high-resolution sample inference from amplicon sequencing data (denoising, filtering, clustering of amplicon sequence variants [ASVs] and assigning taxonomy)
  * [decontam](https://benjjneb.github.io/decontam/): Identifying and removing potential contaminants.
  * [LULU](https://github.com/tobiasgf/lulu): Distribution based post clustering curation of amplicon data.
  * [ShortRead](https://kasperdanielhansen.github.io/genbioconductor/html/ShortRead.html): FASTQ input and manipulation (reading and examining raw sequence reads)
  * [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html): Efficient manipulation of biological strings (reading and examining sequence data)
  * [phyloseq](https://joey711.github.io/phyloseq/): Explore microbiome profiles using R (useful for getting our data into a format that works with many other packages)
  * [tidyverse](https://www.tidyverse.org/): Several key R packages for data science and data manipulation
  * [MicEco](https://github.com/Russel88/MicEco): Various functions for analysis for microbial community data
  * [vegan](https://cran.r-project.org/web/packages/vegan/index.html): Ordination methods, diversity analysis and other functions for community and vegetation ecologists
  * [fantaxtic](https://github.com/gmteunisse/fantaxtic): Nested bar plots for phyloseq data 
  * [here](https://here.r-lib.org/): Package for easy file referencing for reproducible workflows
  * [microbiome](https://microbiome.github.io/tutorials/): Microbiome analytics
  * [ranacapa](https://github.com/gauravsk/ranacapa): We'll use this for some evaluation plots.

# 1. Adapter Trimming 

Normally, you would begin by trimming the Illumina sequencing adapters from the raw sequence data. Luckily for us, the sequencing provider has already done this, and  we can skip this part in our pipeline. The most common tool used for trimming Illumina Sequencing data is [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), which was also used for our data. 

> :memo: **Question 1:** What are sequencing adapters and why do we need to remove them?
  
# 2. Quality Reports

Set the filepath to the location of the raw reads. 
```{bash, eval = F}
cd Data
```

Create a subdirectory to store the FASTQC reports.
```{bash, eval = F}
mkdir ./FASTQC
```

First we have to open the Conda envíronment were the program FASTQC is installed 
```{bash, eval = F}
conda activate fastqc 
```

Create a fastqc report of the read files, in order to examine the read quality.
```{bash, eval = F}
fastqc 211112_SN7280_A_L001_AUXV-6_AdapterTrimmed_R1-0.25.fastq.gz -o ./FASTQC &&
fastqc 211112_SN7280_A_L001_AUXV-6_AdapterTrimmed_R2-0.25.fastq.gz -o ./FASTQC
```

Create a multiqc report to have both quality reports in one html file.
First exit the fastqc conda environment and enter the multiqc environment. 
```{bash, eval = F}
conda deactivate 

conda activate multiqc

multiqc FASTQC -o FASTQC --interactive

conda deactivate
```

> :memo: **Question 2:** What are the important variables in the Quality report? What differences do you notice between the two files?

# 3. Demultiplexing

To demultiplex (distinguish between our samples) we need the octamer tags and primers (reffered to as barcodes from here on) for each sample. We put these on your computers in a folder called BARCODES.

First we need to create a FASTA file from the barcodes so we can use it with a program called CUTADAPT. We take each line of the .txt of the barcodes and add a line containing info on which sample this barcode belongs to.

```{bash, eval = F}
awk '{print">fwd"NR"\n"$0}' ./BARCODES/fungi_fwd_barcodes_big.txt > ./BARCODES/fungi_barcodes_big.fwd.fasta 
```

Take each line of the .txt of the barcodes and add a line containing info on which sample this barcode belongs to.
```{bash, eval = F}
awk '{print">rev"NR"\n"$0}' ./BARCODES/fungi_rev_barcodes_big.txt > ./BARCODES/fungi_barcodes_big.rev.fasta 
```

## Demultiplexing with cutadapt

Create a subdirectory for the file that have the barcodes and primers removed. And enter it.
```{bash, eval = F}
mkdir DEMULTIPLEXED
```

Activate the cutadapt environment through conda. 
```{bash, eval = F}
conda activate cutadapt 
```
Increase the softlimit of the OS because cutadapt will open a lot of files and fill them one after the other. 

One file for each forward and reverse combination. These will be filled with the reads containing the combination.
```{bash, eval = F}
ulimit -S -n 3000
```

## Cutadapt Main Commands
Everything needs to be run twice, because the reads are in mixed orientation, because of the PCR free library preparation. 
Because of the dual indexing approach we need to supply two barcode files. 
```{bash, eval = F}
cutadapt \
-j 0 \
-e 0.15 --no-indels --minimum-length 50 \
-g file:./BARCODES/fungi_barcodes_big.fwd.fasta \
-G file:./BARCODES/fungi_barcodes_big.rev.fasta \
-o ./DEMULTIPLEXED/{name1}-{name2}.round1.1.fastq \
-p ./DEMULTIPLEXED/{name1}-{name2}.round1.2.fastq \
./211112_SN7280_A_L001_AUXV-6_AdapterTrimmed_R1-0.25.fastq.gz \
./211112_SN7280_A_L001_AUXV-6_AdapterTrimmed_R2-0.25.fastq.gz > fungi_round1_cutadapt.txt &&

cutadapt \
-j 0 \
-e 0.15 --no-indels --minimum-length 50 \
-g file:./BARCODES/fungi_barcodes_big.fwd.fasta \
-G file:./BARCODES/fungi_barcodes_big.rev.fasta \
-o ./DEMULTIPLEXED/{name1}-{name2}.round2.1.fastq \
-p ./DEMULTIPLEXED/{name1}-{name2}.round2.2.fastq \
./DEMULTIPLEXED/unknown-unknown.round1.2.fastq \
./DEMULTIPLEXED/unknown-unknown.round1.1.fastq > fungi_round2_cutadapt.txt

# Close the conda environment.
conda deactivate 
```

> :memo: **Question 3:** Which differences do you notice between the two cutadapt commands? And why do we need to use two different commands? 
  
Because of the mixed-orientation of the reads we needed to go trough cutadapt twice. However, the two files of course contain pairs of files belonging to the same biological sample. This is why we need to merge them. 

First we rename the different samples to meaningful names we can attribute to the three regions of the Biodiversity Exploratories and the 50 plots within each of the regions. Additionally it makes sorting easier and we only include reads that are real samples, multiplex controls, blanks and PCR negative controls. The text files necessary to do this indicate which combination of forward and reverse barcodes belong to which sample, and were uploaded to the RENAMING/ directory before. 

We need to enter the directory where we stored the demultiplexed files.
```{bash, eval = F}
cd ./DEMULTIPLEXED
```

Then we rename the files, by taking two text files that contain the old and the new files. Than we loop through the files and change the name of the files by taking the rows content to be the new name. 
```{bash, eval = F}
paste ../RENAMING/renaming_R1_1_old_big.txt \
../RENAMING/renaming_R1_1_new_big.txt \
| while read n k; do mv -v $n $k ; done > ./rename_logR1.1.txt

paste ../RENAMING/renaming_R1_2_old_big.txt \
../RENAMING/renaming_R1_2_new_big.txt \
| while read n k; do mv -v $n $k  ; done > ./rename_logR1.2.txt

paste ../RENAMING/renaming_R2_1_old_big.txt \
 ../RENAMING/renaming_R2_1_new_big.txt \
| while read n k; do mv -v $n $k ; done > ./rename_logR2.1.txt

paste ../RENAMING/renaming_R2_2_old_big.txt \
 ../RENAMING/renaming_R2_2_new_big.txt \
 | while read n k; do mv -v $n $k  ; done > ./rename_logR2.2.txt
```

For a later check of how many reads are passing the pipeline we need to create a FASTA file and count the reads (for a single sample). 

Transform the fastq file to a fasta file. 
```{bash, eval = F}
cat A31_B1.round1.1.sample.fastq | \
awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > A31_B1.round1.1.sample.fa

cat A31_B1.round2.1.sample.fastq | \
awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > A31_B1.round2.1.sample.fa
```

Use grep to get the read numbers by counting the lines beginning with >.
```{bash, eval = F}
grep -c '^>' *.fa | less 
```

> :memo: **Question 4:** What is the difference between FASTQ and FASTA files? 
> :memo: **Question 5:** How many reads are in sample A31_B1? And how many reads are in Sample H12_B3?
  
Now merge the two files from the two passes through cutadapt. 

First we create a new sub-directory for the merged files. 
```{bash, eval = F}
mkdir MERGED
```

Then we create two lists of filenames for the pairs from the cutadapt results. 
```{bash, eval = F}
ls -1 *round1.1.sample.fastq | sed 's/round1.1.sample.fastq//' > listround1.1.

ls -1 *round2.1.sample.fastq | sed 's/round2.1.sample.fastq//' > listround2.1.
```

Now we can merge the pairs by just pasting them together. We start with the files from R1.  
```{bash, eval = F}
paste listround1.1. listround2.1. | while read n k; \
do cat $n"round1.1.sample.fastq" $k"round2.1.sample.fastq" > ./MERGED/$n"sample_demux.1.fastq"; done
```

And then we do the same for the R2 reads. 
```{bash, eval = F}
ls -1 *round1.2.sample.fastq | sed 's/round1.2.sample.fastq//'  > listround1.2.

ls -1 *round2.2.sample.fastq | sed 's/round2.2.sample.fastq//' > listround2.2.

paste listround1.2. listround2.2. | while read n k; \
do cat $n"round1.2.sample.fastq" $k"round2.2.sample.fastq" > ./MERGED/$n"sample_demux.2.fastq"; done
```

In order to check if the merging has worked, we create a FASTA file of the merged sample and check if the read numbers match the paired files we checked before. 

Transform the fastq file to a fasta file. 
```{bash, eval = F}
cat ./MERGED/A31_B1.sample_demux.1.fastq | \
 awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > A31_B1.sample_demux.1.fa
```

Use grep to get the read numbers by counting the lines beginning with >.
```{bash, eval = F}
grep -c '^>' *sample*.fa | less 
```

Exit the DEMULTIPLEXED sub-directory. 
```{bash, eval = F}
cd ..
```

# 4. Checking for leftover primer sequences.

To run DADA2 we need to first activate R. 

We will activate and use R inside the RStudio environment, which is a GUI wrapper for writing and running R code (and other types of code as well).

Once we've opened RStudio, we first load all required packages. 
```{r, eval = F}
library(dada2) ; packageVersion('dada2')
library(ShortRead) ; packageVersion('ShortRead')
library(Biostrings) ; packageVersion('Biostrings')
```

We create objects that contain the primer sequences.  
For the fungi that is:
```{r, eval = F}
FWD <- 'GTGARTCATCGAATCTTTG'
REV <- 'TCCTCCGCTTATTGATATGC'
```

Make a custom function that creates all the possible orientations of the primers e.g. complement, reverse complement.
```{r, eval = F}
allOrients <- function(primer) {
  require(Biostrings)
  # Create all orientations of the input sequence
  dna     <- DNAString(primer)  # turn character to DNAString object
  orients <- c(Forward=dna, Complement=complement(dna), Reverse=reverse(dna),
               RevComp=reverseComplement(dna))
  return(sapply(orients, toString))  # back to character vector
}
```

Make objects that contain the orientations.
```{r, eval = F}
FWD.orients <- allOrients(FWD)
FWD.orients

REV.orients <- allOrients(REV)
REV.orients
```

Set the file paths to the demultiplexed files. 
```{r, eval = F}
fnFs <- sort(list.files(path = './Data/DEMULTIPLEXED/MERGED', pattern = "sample_demux.1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path = './Data/DEMULTIPLEXED/MERGED', pattern = "sample_demux.2.fastq", full.names = TRUE))
```

Dada2 cannot deal with ambiguous bases which are for example contauined in our primers. So we filter out ambiguous Ns with the filterAndTrim function setting maxN to zero.
Place the N-filterd files into a filtN/ subdirectory.
```{r, eval = F}
fnFs.filtN <- file.path(path = './Data/DEMULTIPLEXED/MERGED', "filtN", basename(fnFs)) 
fnRs.filtN <- file.path(path = './Data/DEMULTIPLEXED/MERGED', "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)
```

Check for any leftover primers after the demultiplexing with Cutadapt.

Create a function that counts the number of reads in which the primer is found.
```{r, eval = F}
primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}
```

Search through all the reads and combine in a table.
If the samples come from the same library prep then it is enough to only process one of the files 
(see the [1] at the end of the command). 
```{r, eval = F}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))
```
  
> :memo: **Question 6:** How many reads still contain primer sequences? 

Next we move out of R and back to the terminal. 
We now go through a second pass of cutadapt in order to remove the remaining primers. 

Open the cutadapt environment.
Remove leftover primers with Cutadapt.
```{bash, eval = F}
conda activate cutadapt
```

Create a subdirectory for the first primer removal.
```{bash, eval = F}
mkdir PRIMER_REMOVED1 
```

Enter the MERGED subdirectory.
```{bash, eval = F}
cd DEMULTIPLEXED/MERGED/
```

The command below contains all possible orientations, rc stands for "reverse complement": 
  1) fwd-rcrev + rev-rcfwd; 2) rcfwd-rev + rcrev-fwd; \
  3) fwd + rcfwd; 4) rcfwd + fwd; 5) rev + rcrev; 6) rcrev + rev

Use a for loop to run Cutadapt over all samples. 
```{bash, eval = F}
ls *demu*.fastq | cut -f1 -d'.' > samples

for sample in $(cat samples); do

echo "On sample: $sample"

cutadapt --cores=0 --minimum-length 50 \
 -a ^GTGARTCATCGAATCTTTG...GCATATCAATAAGCGGAGGA -A ^TCCTCCGCTTATTGATATGC...CAAAGATTCGATGAYTCAC \
 -a ^CAAAGATTCGATGAYTCAC...TCCTCCGCTTATTGATATGC -A ^GCATATCAATAAGCGGAGGA...GTGARTCATCGAATCTTTG \
 -a GTGARTCATCGAATCTTTG -A CAAAGATTCGATGAYTCAC \
 -a CAAAGATTCGATGAYTCAC -A GTGARTCATCGAATCTTTG \
 -a TCCTCCGCTTATTGATATGC -A GCATATCAATAAGCGGAGGA \
 -a GCATATCAATAAGCGGAGGA -A TCCTCCGCTTATTGATATGC \
 -o ../../PRIMER_REMOVED1/${sample}.sample_demux_prirm.1.fastq \
 -p ../../PRIMER_REMOVED1/${sample}.sample_demux_prirm.2.fastq \
 ${sample}.sample_demux.1.fastq ${sample}.sample_demux.2.fastq \
 > cutadapt_primer_trimming_stats_{sample}.txt

done 
```

Close the cutadapt environment.
```{bash, eval = F}
conda deactivate 
```

Exit the subdirectory and go to our base directory for the fungal reads.
```{bash, eval = F}
cd ../..
```
> :memo: **Question 7:** How does the cutadapt call for the demultiplexing step differ from the primer removal? What does the flag --cores=0 indicate (hint: Look at the cutadapt website/manual)?
  
Now we check for primers again.  For that we switch to RStudio again.

```{r, eval = F}
# Load in the demultiplexed files. 

fnFs <- sort(list.files(path = './Data/PRIMER_REMOVED1', pattern = "sample_demux_prirm.1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path = './Data/PRIMER_REMOVED1', pattern = "sample_demux_prirm.2.fastq", full.names = TRUE))

# Filter out ambiguous Ns with the filterAndTrim function setting maxN to zero.
# Place the N-filterd files into a filtN/ subdirectory.
fnFs.filtN <- file.path(path = './Data/PRIMER_REMOVED1', "filtN", basename(fnFs)) 
fnRs.filtN <- file.path(path = './Data/PRIMER_REMOVED1', "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

# Check for any leftover primers after the removal with Cutadapt.

# Search through all the reads and combine in a dataframe.
# If the samples come from the same library prep then it is enough to only process one of the files 
# (see the [1] at the end of the command). 
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

```

> :memo: **Question 8:** Was the primer removal successful? Or do we need to go through cutadapt again? 

# 5. Sample inference with DADA2

Load in the files. 
```{r, eval = F}
# For the fungi we need to use the _prirm.fastq files. 

cutFs_1 <- sort(list.files("./Data/PRIMER_REMOVED1", pattern = "1.sample_demux_prirm.1.fastq", full.names = TRUE))
cutRs_1 <- sort(list.files("./Data/PRIMER_REMOVED1", pattern = "1.sample_demux_prirm.2.fastq", full.names = TRUE))

cutFs_2 <- sort(list.files("./Data/PRIMER_REMOVED1", pattern = "2.sample_demux_prirm.1.fastq", full.names = TRUE))
cutRs_2 <- sort(list.files("./Data/PRIMER_REMOVED1", pattern = "2.sample_demux_prirm.2.fastq", full.names = TRUE))

cutFs_3 <- sort(list.files("./Data/PRIMER_REMOVED1", pattern = "3.sample_demux_prirm.1.fastq", full.names = TRUE))
cutRs_3 <- sort(list.files("./Data/PRIMER_REMOVED1", pattern = "3.sample_demux_prirm.2.fastq", full.names = TRUE))

# Create a function to obtain the sample names. "1.s" serves as the point at which the string will be split. 
get.sample.name <- function(fname) strsplit(basename(fname), "1.s")[[1]][1]

# Get the sample names.
sample.names <- unname(sapply(cutFs_1, get.sample.name))
head(sample.names)

```

Filter and trim the reads for quality.

```{r, eval = F}
# Assign filenames for the output of the filtered reads. 
filtFs_1 <- file.path("./Data/PRIMER_REMOVED1", "filtered", basename(cutFs_1))
filtRs_1 <- file.path("./Data/PRIMER_REMOVED1", "filtered", basename(cutRs_1))

filtFs_2 <- file.path("./Data/PRIMER_REMOVED1", "filtered", basename(cutFs_2))
filtRs_2 <- file.path("./Data/PRIMER_REMOVED1", "filtered", basename(cutRs_2))

filtFs_3 <- file.path("./Data/PRIMER_REMOVED1", "filtered", basename(cutFs_3))
filtRs_3 <- file.path("./Data/PRIMER_REMOVED1", "filtered", basename(cutRs_3))

# Apply the filtering parameters , no truncLen because these are ITS reads
# and therefore very variable in length, changed maxEE to 6,6 since the libraries are in mixed orientation.
out_1 <- filterAndTrim(cutFs_1, filtFs_1, cutRs_1, filtRs_1, maxN = 0, maxEE = c(6, 6), 
    truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)

out_2 <- filterAndTrim(cutFs_2, filtFs_2, cutRs_2, filtRs_2, maxN = 0, maxEE = c(6, 6), 
    truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
	
out_3 <- filterAndTrim(cutFs_3, filtFs_3, cutRs_3, filtRs_3, maxN = 0, maxEE = c(6, 6), 
    truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)	

# Check if a good number of reads passed the quality filtering. We are filtering out ~30% which is okay. 
head(out_1)
head(out_2)
head(out_3)

```
> :memo: **Question 9:** What is the function of the truncLen and maxEE parameters in the filterAndTrim function? (Hint: Look at the help file for the function)

Learn the error profiles from the reads. 

```{r, eval = F}
# Learn the error rates for the R1 reads. 
errF_1 <- learnErrors(filtFs_1, multithread = TRUE)

errF_2 <- learnErrors(filtFs_2, multithread = TRUE)

errF_3 <- learnErrors(filtFs_3, multithread = TRUE)

# Learn the error rates for the R2 reads. 
errR_1 <- learnErrors(filtRs_1, multithread = TRUE)

errR_2 <- learnErrors(filtRs_2, multithread = TRUE)

errR_3 <- learnErrors(filtRs_3, multithread = TRUE)

```

> :memo: **Question 10:** Why do we need to do everything three times? 

To have our samples correctly named in the final ASV table we need to rename the filtered files. 
```{r, eval = F}
# Name the filtered objects. 
names(filtFs_1) <- sample.names
names(filtRs_1) <- sample.names

names(filtFs_2) <- sample.names
names(filtRs_2) <- sample.names

names(filtFs_3) <- sample.names
names(filtRs_3) <- sample.names
```

Now we run the DADA2 core for sample inference (i.e. call the ASVs) from sequences within samples. 
```{r, eval = F}
dadaFs_1 <- dada(filtFs_1, err = errF_1, multithread = TRUE)
dadaRs_1 <- dada(filtRs_1, err = errR_1, multithread = TRUE)

dadaFs_2 <- dada(filtFs_2, err = errF_2, multithread = TRUE)
dadaRs_2 <- dada(filtRs_2, err = errR_2, multithread = TRUE)

dadaFs_3 <- dada(filtFs_3, err = errF_3, multithread = TRUE)
dadaRs_3 <- dada(filtRs_3, err = errR_3, multithread = TRUE)
```

Afterwards we need to merge the forward and reverse reads within the technical replicates. The function does this by aligning the reads with each other. Essentially this takes our shorter reads and combines them into a longer sequence that covers the ITS2 region.
```{r, eval = F}
mergers_1 <- mergePairs(dadaFs_1, filtFs_1, dadaRs_1, filtRs_1, verbose=TRUE)

mergers_2 <- mergePairs(dadaFs_2, filtFs_2, dadaRs_2, filtRs_2, verbose=TRUE)

mergers_3 <- mergePairs(dadaFs_3, filtFs_3, dadaRs_3, filtRs_3, verbose=TRUE)
```

We then finally construct the ASV table for each replicate. An ASV table is a sample x ASV table with counts of how often the ASV occurs in each sample (replicate in this case). 
```{r, eval = F}
seqtab_1 <- makeSequenceTable(mergers_1)
dim(seqtab_1)

seqtab_2 <- makeSequenceTable(mergers_2)
dim(seqtab_2)

seqtab_3 <- makeSequenceTable(mergers_3)
dim(seqtab_3)

```
> :memo: **Question 11:** How many ASVs do we find in each of the replicates?

As a last step that is being done to the individual replicates we remove chimeric sequences that sometimes occur during PCR. 

```{r, eval = F}
# Chimera removal for the replicates. 
seqtab_1.nochim <- removeBimeraDenovo(seqtab_1, method="consensus", multithread=TRUE, verbose=TRUE)

seqtab_2.nochim <- removeBimeraDenovo(seqtab_2, method="consensus", multithread=TRUE, verbose=TRUE)

seqtab_3.nochim <- removeBimeraDenovo(seqtab_3, method="consensus", multithread=TRUE, verbose=TRUE)

```

> :memo: **Question 12:** What are chimeric sequences? 

Now we can merge the resulting tables so we have one ASV table for all fungi from the three technical replicates. 

```{r, eval = F}
# Merge the two ASV tables.
# Put the tables in a list before merging them with mergeSequenceTables.
input_tables <- list(seqtab_1.nochim, seqtab_2.nochim, seqtab_3.nochim)

seqtab_merge <- mergeSequenceTables(tables = input_tables, repeats = 'sum', tryRC = TRUE)

# Just for good practice we do an additional run to remove chimeric sequences. 
seqtab_merge.nochim <- removeBimeraDenovo(seqtab_merge, method="consensus", multithread=TRUE, verbose=TRUE)
```

> :memo: **Question 13:** How many ASVs do we have in the merged ASV table?

We want to inspect the sequence lengths to look for anything that seems suspicious. That could be a lot of very short sequences for example, but this is not the case here.
```{r, eval = F}
table(nchar(getSequences(seqtab_merge.nochim)))

# Give the sequence variants more manageable names, e.g. ASVn  
asv_seqs <- colnames(seqtab_merge.nochim)
asv_headers <- vector(dim(seqtab_merge.nochim)[2], mode="character")

for (i in 1:dim(seqtab_merge.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

```

> :memo: **Question 14:** What is the range of the sequence length of our amplicons? 
  
Track our reads through the pipeline so we can see how many reads are filtered out. 
```{r, eval = F}
# Track the reads through the pipeline.
# This is a last sanity check to see if we are not loosing samples at an unexpected step. 

getN <- function(x) sum(getUniques(x))
track_1 <- cbind(out_1, sapply(dadaFs_1, getN), sapply(dadaRs_1, getN), sapply(mergers_1, getN), rowSums(seqtab_1.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track_1) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_1) <- sample.names

track_2 <- cbind(out_2, sapply(dadaFs_2, getN), sapply(dadaRs_2, getN), sapply(mergers_2, getN), rowSums(seqtab_2.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track_2) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_2) <- sample.names

track_3 <- cbind(out_3, sapply(dadaFs_3, getN), sapply(dadaRs_3, getN), sapply(mergers_3, getN), rowSums(seqtab_3.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track_3) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_3) <- sample.names

head(track_1)
head(track_2)
head(track_3)

```

> :memo: **Question 15:** How many sequences did we loose at each step for the first six ASVs occuring in the table? 
  
Now we go on to save our data so we could load it in later. 

```{r, eval = F}
# Make and save a fasta of our final ASV sequences.
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs_fungi.fa")

# Write our chimera screened and merge sequenc table to a text file. This will be our ASV table to work with from here. 
write.table(t(seqtab_merge.nochim), "asv_table_fungi.txt", sep="\t")

# Also save it as an .rds file to not corrupt any structure when we need to reload the table.
saveRDS(seqtab_merge.nochim, 'asv_table_fungi.rds')

```

# 6. Taxonomy assignment 

DADA2 also has built in functionality to assign taxonomy to our ASVs. To do so we need a database that contains sequence information and the corresponding information on taxonomy like the different taxonomic ranks. One of these, and the most commonly used for fungi, is the [UNITE](https://unite.ut.ee/index.php) database. Unfortunately we cannot run these commands on the local computers, because the process is too memory-intensive. We have run the command for you and you can find it in the folder Databases.

Here is how you would do the taxonomy assignment. :no_entry:*PLEASE DO NOT RUN THIS*:no_entry:
```{r, eval = F}
# Read in the UNITE database fasta.  
# unite.ref <- './sh_general_release_dynamic_all_25.07.2023.fasta'  

# Run the taxonomy assignment on the ASV table. 
# taxa <- assignTaxonomy(seqtab_merge.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)

# Save the taxonomy table as an R object.
# saveRDS(taxa, 'tax_table_fungi.rds')

```

# 7. ASV curation and removal of potential contaminants

Before we can start with the curation and decontamination we need to load additional packages. 

```{r, eval = F}
library(here)

library(decontam)

library(phyloseq)

library(lulu)

library(Biostrings)

library(tidyverse)

```

## Decontamination

To look for potential contaminants in our dataset we use a tool called decontam. The algorithm looks for sequences contained in our PCR negative controls and decides if these are contaminants based on the frequency (how often does the sequence occur) and prevalence (in how many samples does it occur).  

Load in the data.

```{r, eval = F}
# Load in the ASV table for the fungi.
fungi_asv <- readRDS(here("asv_table_fungi.rds")) %>%
t() %>%
as.data.frame() %>% 
tibble::rownames_to_column(var = "sequence_fungi")

# Load in the FASTA file so we can get the ASV_IDs.
fungi_seqs_fasta <- Biostrings::readDNAStringSet(here::here('ASVs_fungi.fa'))

# Then we combine the sequences and their corresponding ASV IDs in a dataframe. 
seq_name_fungi <- base::names(fungi_seqs_fasta)
sequence_fungi <- base::paste(fungi_seqs_fasta)
fungi_rep_seqs <- base::data.frame(seq_name_fungi, sequence_fungi)

# Join the ASV table with the representative sequences from the FASTA file.
fungi_asv_IDs <- dplyr::left_join(fungi_rep_seqs, fungi_asv, by = 'sequence_fungi')

# Naming the Samples. 
fungi_asv_IDs <- fungi_asv_IDs %>% 
  dplyr::rename_with(~base::paste0("Sample_", .), -c(1:2))

# Set the ASV ID as the rownames.
base::rownames(fungi_asv_IDs) <- fungi_asv_IDs$seq_name_fungi
fungi_asv_IDs$seq_name_fungi <- NULL
fungi_asv_IDs$sequence_fungi <- NULL

# Load the DNA concentration which we will need to run the decontamination with decontam. 
fungi_conc <- utils::read.csv(here("Data", 'fungi_decontam_conc_big.csv'), header = T, sep =',')
base::rownames(fungi_conc) <- base::paste0("Sample_", fungi_conc$sample_ID)
base::rownames(fungi_conc) <- base::sub("[-]", "_", x = base::rownames(fungi_conc))

```

Many tools for the analysis of microbiome sequencing data are compatible with so-called phyloseq objects. Phlyoseq is a package that manages and combines data such as the ASV table, taxonomy table and metadata. Decontam is one of the packages that accept phyloseq objects. So we combine the data in such an object.

```{r, eval = F}
# Make the parts of the phyloseq object.
ASV_mat_decontam_fun <- base::data.matrix(fungi_asv_IDs[complete.cases(fungi_asv_IDs),])
ASV_fungi_decontam <- phyloseq::otu_table(ASV_mat_decontam_fun, taxa_are_rows = TRUE)
sampledata_fungi_decontam <- phyloseq::sample_data(fungi_conc)

# Combine with phyloseq
ps_fungi_decontam <- phyloseq::phyloseq(ASV_fungi_decontam, sampledata_fungi_decontam)
ps_fungi_decontam
```

> :memo: **Question 16:** How many ASVs does our dataset contain before filtering out contaminants?

In a first step we can check the library sizes of the samples and our controls visually. Library sizes indicate how many reads our samples contain. 

```{r, eval = F}
df_fungi <- base::as.data.frame(phyloseq::sample_data(ps_fungi_decontam)) # Put sample_data into a ggplot-friendly data.frame
df_fungi$LibrarySize <- phyloseq::sample_sums(ps_fungi_decontam)
df_fungi <- df_fungi[base::order(df_fungi$LibrarySize),]
df_fungi$Index <- base::seq(base::nrow(df_fungi))
ggplot2::ggplot(data=df_fungi, ggplot2::aes(x=Index, y=LibrarySize, color=Sample_or_Control)) +
  ggplot2::geom_point()
```

Now we can run the decontamination algorithm. We are using the combined approach of frequency and prevalence.

```{r, eval = F}
phyloseq::sample_data(ps_fungi_decontam)$is.neg <- phyloseq::sample_data(ps_fungi_decontam)$Sample_or_Control == "Control Sample"
contam_fungi_combi <- decontam::isContaminant(ps_fungi_decontam,
                                    method = 'combined',
                                    neg = 'is.neg',
                                    conc = 'quant_reading')
base::table(contam_fungi_combi$contaminant)
base::which(contam_fungi_combi$contaminant)

# Remove the identified contaminants from the phyloseq object 
ps_fungi_noncontam <- phyloseq::prune_taxa(!contam_fungi_combi$contaminant,
                                 ps_fungi_decontam)
ps_fungi_noncontam

# Remove taxa without reads.
ps_fungi_noncontam_pruned <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_fungi_noncontam) > 0,
                                        ps_fungi_noncontam)
```

> :memo: **Question 17:** Which taxa were removed from the original datasets as contaminants? 
> :memo: **Question 18:** How many taxa remain after this step?
> :memo: **Question 19:** What is the difference between the frequency, prevalence and combined approach in the decontam function?
  
## LULU curation

### Match List 

We will run a curation with LULU, a so-called post clustering algorithm. It basically takes similar sequences and decides if they are similar enough to be the same biological entity. It works on their co-occurence pattern in our dataset. 

We don't want it to compare each sequence with eachother. So we need to make a table of which sequences are similar to eacxh other above a certain similarity threshold. This table is called a matchlist.

We move to the terminal.

```{bash, eval = F}
# Open the conda environment containing blast. 
conda activate blast

# Move into the Data folder.
cd Data

# Create a database containing the algal ASVs in fasta format.
makeblastdb -in ASVs_fungi.fa -parse_seqids -dbtype nucl

# Compare all ASVs to each other and create a match list of sequences that are similar to each other.
blastn -db ASVs_fungi.fa \
-outfmt '6 qseqid sseqid pident' \
-out match_list_fungi.txt \
-qcov_hsp_perc 80 \
-perc_identity 84 \
-num_threads 6 \
-query ~/ASVs_fungi.fa

conda deactivate

```

Now that we have obtained our matchlist we move back into RStudio and run the actual LULU curation.

```{r, eval = F}
ASV_table_fungi <- base::as.data.frame(phyloseq::otu_table(ps_fungi_noncontam_pruned))
fungi_matchlist <- utils::read.table(here::here("Data", 'match_list_fungi.txt'),
                              header = F,
                              as.is = T,
                              stringsAsFactors = F)

# Run the LULU algorithm. 

ASV_table_fungi_cur <- lulu::lulu(ASV_table_fungi, fungi_matchlist)

ASV_table_fungi_cur$curated_count
ASV_table_fungi_cur$discarded_count

# Save the table. 
saveRDS(ASV_table_fungi_cur, here("Data", "ASV_table_fungi_cur.rds"))


```

> :memo: **Question 20:** How many reads were merged? 
> :memo: **Question 21:** Which parameters does the LULU algorithm consider in it's merging? 
> :memo: **Question 22:** What are the default parameters? 
  
# Diversity Analysis

Now we have most of the technical processing of the sequencing results done and can finally move on to answering question. In this course we will care mostly about community composition and how it differs between some environmental groups. We will work with the data from one of the Biodiversity Exploratories: The Biosphere Reserve Schorfheide-Chorin. 

## 8. Data initialization 

Before we can actually get into the questions we need to do a few last steps that have to do with cleaning our data so it is nice and presentable. It also makes it easier to work with. 

First we load in the taxonomy table we have created with DADA2. 

```{r, eval = F}
# Load in taxonomy table. 
tax_fungi <- base::readRDS(here::here("Data", "tax_table_fungi.rds")) %>% 
  base::as.data.frame() %>%
  tibble::rownames_to_column('sequence') %>%
  dplyr::rename(sequence_fungi = sequence)

# Load the fungal reads.
fungi_seqs_fasta <- Biostrings::readDNAStringSet(here::here("Data",'ASVs_fungi.fa'))

# Make a dataframe of the sequences and their ASV ID. 
seq_name_fungi <- base::names(fungi_seqs_fasta)
sequence_fungi <- base::paste(fungi_seqs_fasta)
fungi_rep_seqs <- base::data.frame(seq_name_fungi, sequence_fungi)

# Join the taxonomy table and the representative sequences
tax_clean_fungi <- dplyr::left_join(tax_fungi, fungi_rep_seqs, by = 'sequence_fungi') %>% drop.na(seq_name_fungi)

```

To make the taxonomy table nicer to look at and easier to work with we need to do some cleaning. The syntax of the UNITE database we used includes the information on taxonomy in a different format than we would like. We want to remove the underscores and have one column per taxonomic level. 

```{r, eval = F}
# Have a look at the original table.
View(tax_clean_fungi)

# Split the taxonomy into different columns of taxonomic levels.
fungi_tax_fin <- tidyr::separate(tax_clean_fungi, Kingdom, c(NA, 'Kingdom') , sep = '__') %>% 
  tidyr::separate(Phylum, c(NA, 'Phylum') , sep = '__') %>% 
  tidyr::separate(Class, c(NA, 'Class') , sep = '__') %>% 
  tidyr::separate(Order, c(NA, 'Order') , sep = '__') %>% 
  tidyr::separate(Family, c(NA, 'Family') , sep = '__') %>% 
  tidyr::separate(Genus, c(NA, 'Genus') , sep = '__') %>% 
  tidyr::separate(Species, c(NA, 'Species') , sep = '__')

# Rename the ASV_ID column. 
fungi_tax_fin <- dplyr::rename(fungi_tax_fin, ASV_ID = seq_name_fungi)

# Set rownames.
base::row.names(fungi_tax_fin) <- fungi_tax_fin$ASV_ID

fungi_tax_fin$sequence_fungi <- NULL
fungi_tax_fin$ASV_ID <- NULL
```

The ITS is a genetic marker not only found in fungi, so we need to remove amplicons we are not interested in. 
Before you do the removal have a look at which other things we assigned.
```{r, eval = F}
fungi_tax_fin <- fungi_tax_fin %>% 
  dplyr::filter(Kingdom == "Fungi")
```
> :memo: **Question 23:** List three other things that you found that are not fungi. 

We also need some metadata with information on the exploratories, tree species etc.. 

```{r, eval = F}
metadata <- utils::read.csv(here::here("Data",'sample_data_3_AHS.csv'), sep = ';')

# Set the sample name as the rowname for the phyloseq creation.
base::row.names(metadata) <- metadata$sample
metadata$sample <- NULL

```
In a next step we load in the ASV table we curated with LULU. 

```{r, eval = F}
ASV_table_fungi_cur <- base::readRDS(here::here("Data", "ASV_table_fungi_cur.rds"))

# Keep only samples that do represent real tree swabs. Cut Controls. 
asv_fungi <- ASV_table_fungi_cur$curated_table %>% 
  dplyr::select(all_of(base::rownames(metadata)))
```

Combine the data into a phyloseq object.  

```{r, eval = F}
# First create a phyloseq object with all samples included.
# Create the matrices needed for the phyloseq functions.  
asv_mat <- base::as.matrix(asv_fungi)
taxmat <- base::as.matrix(fungi_tax_fin) 
  
# Create the phyloseq object. 
ASV <- phyloseq::otu_table(asv_mat, taxa_are_rows = TRUE)
TAX <- phyloseq::tax_table(taxmat)
sampledata <- phyloseq::sample_data(metadata)

full_physeq <- phyloseq::phyloseq(ASV, TAX, sampledata)

```

> :memo: **Question 24:** How many taxa does our phyloseq object contain now? 
  
For our analysis we want to focus on the differences between soil and bark samples as well as differences between coniferous and deciduous trees. So we need to do a lot of splitting and filtering of the phyloseq objects. We will only work within one region of the biodiversity exploratories, Schorfheide-Chorin in North-East Germany.  

```{r, eval = F}
# Filter out samples that were sampled on trees that are not from our two target species and from the other . 

physeq_sch <- phyloseq::subset_samples(full_physeq, 
                                            dominant_tree %in% c("Fagus_sylvatica",
                                                                 "Pinus_sylvestris")) %>% 
  phyloseq::subset_samples(exploratory == "Schorfheide") %>%
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

# Split the data into soil and bark samples.
physeq_sch_bark <- phyloseq::subset_samples(physeq_sch, substrate == "bark") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_sch_soil <- phyloseq::subset_samples(physeq_sch, substrate == "soil") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

# Split the phyloseqobject by tree species. 
# Schorfheide Chorin
physeq_sch_fagus <- phyloseq::subset_samples(physeq_sch, dominant_tree == "Fagus_sylvatica") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_sch_pinus <- phyloseq::subset_samples(physeq_sch, dominant_tree == "Pinus_sylvestris") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

# Split the phyloseq object into the the substrate & the tree species.
physeq_sch_bark_fagus <- phyloseq::subset_samples(physeq_sch_bark, dominant_tree == "Fagus_sylvatica") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_sch_soil_fagus <- phyloseq::subset_samples(physeq_sch_soil, dominant_tree == "Fagus_sylvatica") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_sch_bark_pinus <- phyloseq::subset_samples(physeq_sch_bark, dominant_tree == "Pinus_sylvestris") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

physeq_sch_soil_pinus <- phyloseq::subset_samples(physeq_sch_soil, dominant_tree == "Pinus_sylvestris") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

```

> :memo: **Question 25:** How many taxa do the different splits contain? Make a little table. 
  
# 9. Analysis of library sizes

We already have had a short look at how many reads we have per sample when we did the decontamination. To have a deeper look at the library sizes, and to judge if we have sequenced deep enough to find rare species, we can look at so-called rarefaction curves.

```{r; eval = F}
# Create a sample_data column containing a column that corresponds to tree species and substrate.
physeq_sch_curve <- physeq_sch
phyloseq::sample_data(physeq_sch_curve) <- phyloseq::sample_data(physeq_sch) %>% 
  base::data.frame() %>%  
  dplyr::mutate(tree_substrate = base::paste(dominant_tree, substrate, sep = "-"))

# Create rarefaction curve and color the lines by substrate (bark/soil) and host tree species. 
rare_sch <- ranacapa::ggrare(physeq_sch_curve, step = 50,
                             color = "tree_substrate", se = FALSE) +
  theme(legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "cm"))

rare_sch

```

> :memo: **Question 26:** What differences can we see from these rarefaction curves?

# 10. Alpha Diversity

A first interesting diversity question might be if substrates or tree species differ in their number of species (richness) or diversity (e.g. Shannon diversity index). We can easily calculate both of these measures with the build-in phyloseq functions and then visually compare them with simple boxplots. 

First, we calculate the diversity measures and add them to the phyloseq object so that we have access to all other sample information as well.
```{r, eval = F}
# Calculate the diversity measures.
div_fun <- estimate_richness(physeq_sch, split = T, measures = c('Observed', 'Shannon'))

# Make a new sample_data object and merge with existing phyloseq object.
new_var_fun <- sample_data(div_fun)
physeq_sch <- merge_phyloseq(physeq_sch, new_var_fun)

fungi_sampledata <- sample_data(physeq_sch)
```

Now we can use the data to create several boxplots. First, we will look at the differences in tree species.
```{r, eval = F}
fun_box_tree_shannon <- ggplot(fungi_sampledata, aes(x = dominant_tree, y = Shannon, fill = dominant_tree, colour = dominant_tree)) + 
  geom_boxplot(
    width = .15, 
    ## remove outliers
    outlier.color = NA, ## `outlier.shape = NA` works as well 
    alpha = 0.5,
    size = 0.1
  )
fun_box_tree_shannon  

fun_box_tree_richness <- ggplot(fungi_sampledata, aes(x = dominant_tree, y = Observed, fill = dominant_tree, colour = dominant_tree)) + 
  geom_boxplot(
    width = .15, 
    ## remove outliers
    outlier.color = NA, ## `outlier.shape = NA` works as well 
    alpha = 0.5,
    size = 0.1
  )
fun_box_tree_richness  
```  
We can also look at the differences between substrates.
```{r, eval = F}
fun_box_substrate_shannon <- ggplot(fungi_sampledata, aes(x = substrate, y = Shannon, fill = substrate, colour = substrate)) + 
  geom_boxplot(
    width = .15, 
    ## remove outliers
    outlier.color = NA, ## `outlier.shape = NA` works as well 
    alpha = 0.5,
    size = 0.1
  )
fun_box_substrate_shannon  

fun_box_substrate_richness <- ggplot(fungi_sampledata, aes(x = substrate, y = Observed, fill = substrate, colour = substrate)) + 
  geom_boxplot(
    width = .15, 
    ## remove outliers
    outlier.color = NA, ## `outlier.shape = NA` works as well 
    alpha = 0.5,
    size = 0.1
  )
fun_box_substrate_richness  
```
> :memo: **Question 27:** What differences can you see? Does host species or substrate have a greater effect on species richness or diversity?

  
# 11. Community Composition

Some fundamental questions in community ecology are: Which taxa occur in a particular dataset? Which taxa are more specific (e.g. unique to) each sample? And are there differences? 

To get a first clue as to the answer to these questions, we can create a community composition bar plot. This will give us a visual representation of what is in our sample, and whether or not one group occurs more often than another. We plot this visualization at the "Order" taxonomic level, because we would not have enough information to clearly resolve differences otherwise.

For this purpose, we need to transform our read counts into relative abundances (the proportion of how often a taxon occurs in the data) and find the 25 most abundant orders. 

First we create a phyloseq object that is aggregated to the top 25 orders. 

```{r, eval = F}
physeq_sch_barplot <- physeq_sch
phyloseq::sample_data(physeq_sch_barplot) <- phyloseq::sample_data(physeq_sch) %>% 
  base::data.frame() %>%  
  dplyr::mutate(tree_substrate = base::paste(dominant_tree, substrate, sep = "-"))


# Subset the phyloseq object to the top 24 orders and put the rest in 
# a category "Others", based on relative abundance. 
phy_sch_ord_top25 <- fantaxtic::top_taxa(physeq_sch_barplot,
                                         tax_level = 'Order',
                                         n_taxa =  24,
                                         by_proportion = TRUE,
                                         merged_label = "Other",
                                         include_na_taxa = T) 
phy_sch_ord_top25_named <- fantaxtic::name_na_taxa(phy_sch_ord_top25$ps_obj, include_rank = T)
```

Afterwards, we transform the read counts from our ASV table to relative abundances. 
```{r, eval = F}
# Transform the subset dataset to compositional (relative) abundances.
phy_sch_ord_top25_named_plot <-  phy_sch_ord_top25_named %>%
  microbiome::aggregate_taxa(level = "Order") %>%  
  microbiome::transform(transform = "compositional") 
```

Then we need to extract the names of the 25 orders.

```{r, eval = F}
# Extract the names of the Orders.
phyloseq::taxa_names(phy_sch_ord_top25_named_plot) <- phyloseq::tax_table(phy_sch_ord_top25_named_plot)[, 4]

# Sort the taxa names alphabetically. 
taxa_names_sch_ord <- sort(phyloseq::taxa_names(phy_sch_ord_top25_named_plot))
```

To get our desired plotting order and group names, we need to change the exploratory names and order them as factors.

```{r, eval = F}
sampledata_sch <- data.frame(phyloseq::sample_data(phy_sch_ord_top25_named_plot))
sampledata_sch <- sampledata_sch %>% 
  mutate(across("tree_substrate", stringr::str_replace, "Fagus_sylvatica-bark", "F. sylvatica bark")) %>% 
  mutate(across("tree_substrate", stringr::str_replace, "Fagus_sylvatica-soil", "F. sylvatica soil")) %>% 
  mutate(across("tree_substrate", stringr::str_replace, "Pinus_sylvestris-bark", "P. sylvestris bark")) %>% 
  mutate(across("tree_substrate", stringr::str_replace, "Pinus_sylvestris-soil", "P. sylvestris soil"))
          
sampledata_sch$tree_substrate <- factor(sampledata_sch$tree_substrate, 
                                        levels = c("F. sylvatica bark", "P. sylvestris bark",
                                                   "F. sylvatica soil", "P. sylvestris soil"))  

phyloseq::sample_data(phy_sch_ord_top25_named_plot) <- phyloseq::sample_data(sampledata_sch)
```

To create some nice plots, we will need to define some colors and a function that fills the bars that we create.

```{r, eval = F}
my_cols <- paletteer::paletteer_d('ggsci::default_igv')

unique_orders <- sort(unique(taxa_names_sch_ord)) 

full_cols <- data.frame(order = unique_orders, color = my_cols[1:25])

sch_cols <- full_cols %>% 
  dplyr::filter(order %in% taxa_names_sch_ord)

scale_fill_sch <- function(...){
  ggplot2:::manual_scale(
    'fill', 
    values = setNames(sch_cols$color, sch_cols$order), 
    ...
  )
}
```

Now we can use functionality from a package called ggplot2. This will let us manipulate a lot of the elements in the plot to suit our visualization needs. We create two plots: One for the soil and one for the bark microbiome. We will then have two subplots within them, one for the beech trees and one for the pine trees.

```{r, eval = F}
sch_ord_soil_plots <- phyloseq::subset_samples(phy_sch_ord_top25_named_plot, substrate == "soil") %>% 
  microbiome::plot_composition(group_by =  'tree_substrate', otu.sort = "abundance") +
  scale_fill_sch() +
  guides(fill = guide_legend(title.position = 'top', ncol = 10)) +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = 'black', size = 0.5),
        axis.text.x =  element_blank(),
        axis.text.y =  element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 10),
        legend.position = 'bottom', 
        plot.title = element_text(vjust = -4, hjust = 0.03), 
        legend.text = element_text(colour = 'black', size = 7),
        legend.title =  element_text(size = 10),
        legend.key.size = unit(2.5, 'mm'),
        axis.ticks.length.x = unit(-0.2, "cm"), 
        legend.box.spacing = unit(-4, 'mm'),
        legend.background = element_rect(fill = 'transparent'),
        text = element_text(colour = 'black', size = 20),
        strip.text = element_text(face = "italic")) + 
  xlab('Sample') +
  ylab('Relative Abundance') 
sch_ord_soil_plots

sch_ord_bark_plots <- phyloseq::subset_samples(phy_sch_ord_top25_named_plot, substrate == "bark") %>% 
  microbiome::plot_composition(group_by =  'tree_substrate', otu.sort = "abundance") +
  scale_fill_sch() +
  guides(fill = guide_legend(title = 'Order',title.position = 'top', ncol = 10)) +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = 'black', size = 0.5),
        axis.text.x =  element_blank(),
        axis.text.y =  element_text(colour = "black", size = 10),
        axis.title.x = element_text(colour = "black", size = 10),
        axis.title.y = element_text(colour = "black", size = 10),
        legend.position = 'bottom',
        text = element_text(colour = 'black', size = 20), 
        plot.title = element_text(vjust = -4, hjust = 0.03), 
        legend.text = element_text(colour = 'black', size = 7),
        legend.title =  element_text("Order",size = 8),
        legend.key.size = unit(1, 'mm'),
        axis.ticks.length.x = unit(-0.2, "cm"), 
        legend.box.spacing = unit(-4, 'mm'),
        strip.text = element_text(face = "italic")) + 
  xlab('Sample') +
  ylab('Relative Abundance') 
sch_ord_bark_plots
```

> :memo: **Question 28:** What are the most abundant orders in each tree species in a) soil and b) bark?

# 12. Differences in community composition/beta diversity between bark and soil

We can see from the previous plots that there are differences between tree species and substrates. We don't know yet if they are really meaningful (i.e. statistically significant), or if just a trend/pattern that we can visually distinguish.

Before we get an answer to that question, however, first let's look at another method to visually represent differences between groups: Ordination. 

Ordinations are based on a distance metric between samples and groups. The distance measure indicates how different the trees are in the abundance and number of ASVs. A common measure that we are going to use here is Bray-Curtis dissimilarity, which looks at the abundances of taxa shared by two samples and the number of taxa found in each sample. If both samples share the same number of ASVs, with the same abundances, the dissimilarity (distance) will be zero. We can than ordinate the samples with an approach called Non-Metric Multidimensional Scaling. This approach ranks how dissimilar samples are to each other and plots them in a 2D space.

```{r, eval = F}
# Ordinate the Schorfheide-Chorin phyloseq using an NMDS with Bray-Curtis distance.
nmds_sch <- phyloseq::ordinate(physeq_sch, method = "NMDS", distance = "bray")

# Plot the ordination and enclose the 95% confidence intervals around the group of tree species.
ordination_sch_tree <- phyloseq::plot_ordination(physeq_sch, nmds_sch, type="samples", color="dominant_tree", shape="substrate") + 
  ggplot2::geom_point(size = 4) +
  ggplot2::stat_ellipse(ggplot2::aes(group = dominant_tree), linetype = 2) +
  ggplot2::scale_colour_manual(values = c("green","darkolivegreen4"), name = "dominant tree species",
                               labels = c("Fagus sylvatica", "Pinus sylvestris")) +
  ggplot2::labs(subtitle = "(B) Schorfheide-Chorin") +
  ggplot2::theme(legend.position = "right",
        axis.text = ggplot2::element_text(size = 15), 
        axis.title = ggplot2::element_text(size = 15))

ordination_sch_tree

# Plot the ordination and enclose the 95% confidence intervals around the group of substrate.
ordination_sch_substrate <- phyloseq::plot_ordination(physeq_sch, nmds_sch, type="samples", color="dominant_tree", shape="substrate") + 
  ggplot2::geom_point(size = 4) +
  ggplot2::stat_ellipse(ggplot2::aes(group = substrate), linetype = 2) +
  ggplot2::scale_colour_manual(values = c("green","darkolivegreen4"), name = "dominant tree species",
                               labels = c("Fagus sylvatica", "Pinus sylvestris")) +
  ggplot2::labs(subtitle = "(B) Schorfheide-Chorin") +
  ggplot2::theme(legend.position = "right",
        axis.text = ggplot2::element_text(size = 15), 
        axis.title = ggplot2::element_text(size = 15))

ordination_sch_substrate

```

With the dissimilarity/distance measure we can now also test if our groups (tree species & substrate) are really different statistically. To do so we will run a PERMANOVA: a permutational multivariate analysis of variance. This will decompose the dissimilarity matrix into variation within and between the groups we are interested in. 

```{r, eval = F}
# We need to calculate a distance matrix first.
bray_mat_sch <- phyloseq::distance(physeq_sch, method = "bray")

# Run the PERMANOVA analysis through vegans adonis2() including both effects of habitat and tree species.
vegan::adonis2(bray_mat_sch ~ factor(phyloseq::sample_data(physeq_sch)$substrate) +
                 factor(phyloseq::sample_data(physeq_sch)$dominant_tree), by = 'margin')

# We can see that both have significant differences. 

```

PERMANOVA results might be confounded if the groups were to have very different distributions around their centers. In order to see if our results are reliable, we therefore have to look at the dispersion.  

```{r, eval = F}
# Test the within group dispersion for the substrate.
dispr_substrate_sch <- vegan::betadisper(bray_mat_sch, 
                                         factor(phyloseq::sample_data(physeq_sch)$substrate))
dispr_substrate_sch

vegan::permutest(dispr_substrate_sch)

# Test the within group dispersion for the tree Species.
dispr_tree_sch <- vegan::betadisper(bray_mat_sch, 
                                    factor(phyloseq::sample_data(physeq_sch)$dominant_tree))
dispr_tree_sch

vegan::permutest(dispr_tree_sch)

```

> :memo: **Question 29:** What other options can you find for distance measures between samples? 
> :memo: **Question 30:** Are there other ordination methods apart from NMDS?  
> :memo: **Question 31:** Are there statistically meaningful differences between substrates? Between tree species? 
