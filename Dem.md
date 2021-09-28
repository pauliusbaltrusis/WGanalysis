## POOL1-5
```shell
lima ../raw_data/ps_365_001.fastq.gz ../tags_primers/pool1_tags/9tags.fasta p1.fastq --split-named --min-score 80 --min-length 150 --same # demultiplexing pool 1 using pool1 tags
...
```
## Optional when creating a fasta DB of sequences
``` shell
sed '2n;s/pattern/replacement/;n' ## In order to add ";" to the ending of every second line starting from the first
```
## Manipulations in R using DADA2
```shell
# DADA2 pipeline for pools starts here

##Install DADA2
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2", version = "3.13")

##Load it and others
library(dada2)
library(ShortRead)
library(Biostrings)

## data import/prep for primer removal/trimming
path1<-'D:/1_pool'
path2<-'D:/2_pool'
path3<-'D:/3_pool'
path4<-'D:/4_pool'
path5<-'D:/5_pool'

## Assign names to dataset
fNFs1<-sort(list.files(path1,pattern='.fastq',full.names = TRUE))
fNFs2<-sort(list.files(path2,pattern='.fastq',full.names = TRUE))
fNFs3<-sort(list.files(path3,pattern='.fastq',full.names = TRUE))
fNFs4<-sort(list.files(path4,pattern='.fastq',full.names = TRUE))
fNFs5<-sort(list.files(path5,pattern='.fastq',full.names = TRUE))

## define primers
FWD<-"ACGTCTGGTTCAGGGTTGTT"
REV<-"TTAGTTTCTTTTCCTCCGCT"

## Function to find primer orientation
allOrients<- function(primer) {
  require(Biostrings)
  dna<- DNAString(primer)
  orients<-c(Forward=dna, Complement=complement(dna), Reverse=reverse(dna),
             RevComp=reverseComplement(dna))
  return(sapply(orients,toString))
}

## Make orientations of both forw and rev

FWD.orients<- allOrients(FWD)
REV.orients<-allOrients(REV)

## Remove reads containing "Ns"

fNFs1.filtN<-file.path(path1,"filtN", basename(fNFs1)) ##Put N-filtered files in filtN/ subdirectory
filterAndTrim(fNFs1, fNFs1.filtN, maxN=0, multithread = F)
fNFs2.filtN<-file.path(path2,"filtN", basename(fNFs2)) ##Put N-filtered files in filtN/ subdirectory
filterAndTrim(fNFs2, fNFs2.filtN, maxN=0, multithread = F)
fNFs3.filtN<-file.path(path3,"filtN", basename(fNFs3)) ##Put N-filtered files in filtN/ subdirectory
filterAndTrim(fNFs3, fNFs3.filtN, maxN=0, multithread = F)
fNFs4.filtN<-file.path(path4,"filtN", basename(fNFs4)) ##Put N-filtered files in filtN/ subdirectory
filterAndTrim(fNFs4, fNFs4.filtN, maxN=0, multithread = F)
fNFs5.filtN<-file.path(path5,"filtN", basename(fNFs5)) ##Put N-filtered files in filtN/ subdirectory
filterAndTrim(fNFs5, fNFs5.filtN, maxN=0, multithread = F)

## Identifying how many primer seqs are in pools1-5
primerHits<-function(primer, fn) {
  # counts no of reads in which the primer is found
  nhits<-vcountPattern(primer, sread(readFastq(fn)), fixed=FALSE)
  return(sum(nhits>0))
}
rbind(FWD.ForwardReads=sapply(FWD.orients, primerHits, fn=fNFs1.filtN[1]),
      FWD.ReverseReads=sapply(FWD.orients, primerHits, fn=fNFs1.filtN[1]),
      REV.ForwardReads=sapply(REV.orients, primerHits, fn=fNFs1.filtN[1]),
      REV.ReverseReads=sapply(REV.orients, primerHits, fn=fNFs1.filtN[1])
      )
rbind(FWD.ForwardReads=sapply(FWD.orients, primerHits, fn=fNFs2.filtN[1]),
      FWD.ReverseReads=sapply(FWD.orients, primerHits, fn=fNFs2.filtN[1]),
      REV.ForwardReads=sapply(REV.orients, primerHits, fn=fNFs2.filtN[1]),
      REV.ReverseReads=sapply(REV.orients, primerHits, fn=fNFs2.filtN[1])
      )
rbind(FWD.ForwardReads=sapply(FWD.orients, primerHits, fn=fNFs3.filtN[1]),
      FWD.ReverseReads=sapply(FWD.orients, primerHits, fn=fNFs3.filtN[1]),
      REV.ForwardReads=sapply(REV.orients, primerHits, fn=fNFs3.filtN[1]),
      REV.ReverseReads=sapply(REV.orients, primerHits, fn=fNFs3.filtN[1])
      )
rbind(FWD.ForwardReads=sapply(FWD.orients, primerHits, fn=fNFs4.filtN[1]),
      FWD.ReverseReads=sapply(FWD.orients, primerHits, fn=fNFs4.filtN[1]),
      REV.ForwardReads=sapply(REV.orients, primerHits, fn=fNFs4.filtN[1]),
      REV.ReverseReads=sapply(REV.orients, primerHits, fn=fNFs4.filtN[1])
      )
rbind(FWD.ForwardReads=sapply(FWD.orients, primerHits, fn=fNFs5.filtN[1]),
      FWD.ReverseReads=sapply(FWD.orients, primerHits, fn=fNFs5.filtN[1]),
      REV.ForwardReads=sapply(REV.orients, primerHits, fn=fNFs5.filtN[1]),
      REV.ReverseReads=sapply(REV.orients, primerHits, fn=fNFs5.filtN[1])
      )
## Note to self: since it's pacbio - a single read per contig - orientation of primers in reads can be bidirectional
### Using cutadapt on Windows

cutadapt<-'D:/cutadapt-3.4.exe'
system2(cutadapt, args='--version')

## Create paths for cutadapt
path1.cut<-file.path(path1,'cutadapt')
if(!dir.exists(path1.cut)) dir.create(path1.cut)

path2.cut<-file.path(path2,'cutadapt')
if(!dir.exists(path2.cut)) dir.create(path2.cut)

path3.cut<-file.path(path3,'cutadapt')
if(!dir.exists(path3.cut)) dir.create(path3.cut)

path4.cut<-file.path(path4,'cutadapt')
if(!dir.exists(path4.cut)) dir.create(path4.cut)

path5.cut<-file.path(path5,'cutadapt')
if(!dir.exists(path5.cut)) dir.create(path5.cut)

## Asign samples to the path for cutadapt
fnFs1.cut<-file.path(path1.cut, basename(fNFs1))
fnFs2.cut<-file.path(path2.cut, basename(fNFs2))
fnFs3.cut<-file.path(path3.cut, basename(fNFs3))
fnFs4.cut<-file.path(path4.cut, basename(fNFs4))
fnFs5.cut<-file.path(path5.cut, basename(fNFs5))

## Reverse complements for primers
FWD.RC<- dada2::rc(FWD)
REV.RC<- dada2::rc(REV)
# Variable for: Trimming FWD and the Reverse-complement of REV off of forward reads
R1.flags<-paste('-g', FWD, '-a', REV.RC)
# Variable for: Trimming REV and the Reverse-complement of FWD off of reverse reads
R2.flags<-paste('-g', REV, "-a", FWD.RC)

# Run cutadapt #1 to remove FWD and REVcomp primers from pools
for(i in seq_along(fNFs1)) {
  system2(cutadapt, args=c(R1.flags, R2.flags, '-n',2, # -n 2 required to remove FWD and REV from reads
                           '-o', fnFs1.cut[i], #output file
                           fNFs1.filtN[i] # input file
                           )) 
}

## Check up on pool1

rbind(FWD.ForwardReads=sapply(FWD.orients, primerHits, fn=fnFs1.cut[[1]]),
      FWD.ReverseReads=sapply(FWD.orients, primerHits, fn=fnFs1.cut[[1]]),
      REV.ForwardReads=sapply(REV.orients, primerHits, fn=fnFs1.cut[[1]]),
      REV.ReverseReads=sapply(REV.orients, primerHits, fn=fnFs1.cut[[1]]))
      

##cutadapt #2

for(i in seq_along(fNFs2)) {
  system2(cutadapt, args=c(R1.flags, R2.flags, '-n',2, # -n 2 required to remove FWD and REV from reads
                           '-o', fnFs2.cut[i], #output file
                           fNFs2.filtN[i] # input file
                           )) 
}

##cutadapt #3

for(i in seq_along(fNFs3)) {
  system2(cutadapt, args=c(R1.flags, R2.flags, '-n',2, # -n 2 required to remove FWD and REV from reads
                           '-o', fnFs3.cut[i], #output file
                           fNFs3.filtN[i] # input file
                           )) 
}

##cutadapt #4

for(i in seq_along(fNFs4)) {
  system2(cutadapt, args=c(R1.flags, R2.flags, '-n',2, # -n 2 required to remove FWD and REV from reads
                           '-o', fnFs4.cut[i], #output file
                           fNFs4.filtN[i] # input file
                           )) 
}

##cutadapt #5

for(i in seq_along(fNFs5)) {
  system2(cutadapt, args=c(R1.flags, R2.flags, '-n',2, # -n 2 required to remove FWD and REV from reads
                           '-o', fnFs5.cut[i], #output file
                           fNFs5.filtN[i] # input file
                           )) 
}

## Check ups for pool2-5
rbind(FWD.ForwardReads=sapply(FWD.orients, primerHits, fn=fnFs2.cut[[1]]),
      FWD.ReverseReads=sapply(FWD.orients, primerHits, fn=fnFs2.cut[[1]]),
      REV.ForwardReads=sapply(REV.orients, primerHits, fn=fnFs2.cut[[1]]),
      REV.ReverseReads=sapply(REV.orients, primerHits, fn=fnFs2.cut[[1]]))

rbind(FWD.ForwardReads=sapply(FWD.orients, primerHits, fn=fnFs3.cut[[1]]),
      FWD.ReverseReads=sapply(FWD.orients, primerHits, fn=fnFs3.cut[[1]]),
      REV.ForwardReads=sapply(REV.orients, primerHits, fn=fnFs3.cut[[1]]),
      REV.ReverseReads=sapply(REV.orients, primerHits, fn=fnFs3.cut[[1]]))

rbind(FWD.ForwardReads=sapply(FWD.orients, primerHits, fn=fnFs4.cut[[1]]),
      FWD.ReverseReads=sapply(FWD.orients, primerHits, fn=fnFs4.cut[[1]]),
      REV.ForwardReads=sapply(REV.orients, primerHits, fn=fnFs4.cut[[1]]),
      REV.ReverseReads=sapply(REV.orients, primerHits, fn=fnFs4.cut[[1]]))

rbind(FWD.ForwardReads=sapply(FWD.orients, primerHits, fn=fnFs5.cut[[1]]),
      FWD.ReverseReads=sapply(FWD.orients, primerHits, fn=fnFs5.cut[[1]]),
      REV.ForwardReads=sapply(REV.orients, primerHits, fn=fnFs5.cut[[1]]),
      REV.ReverseReads=sapply(REV.orients, primerHits, fn=fnFs5.cut[[1]]))
      
## giving names to cut pools

cutFs1<- sort(list.files(path1.cut, pattern = ".fastq", full.names = TRUE))
cutFs2<- sort(list.files(path2.cut, pattern = ".fastq", full.names = TRUE))
cutFs3<- sort(list.files(path3.cut, pattern = ".fastq", full.names = TRUE))
cutFs4<- sort(list.files(path4.cut, pattern = ".fastq", full.names = TRUE))
cutFs5<- sort(list.files(path5.cut, pattern = ".fastq", full.names = TRUE))

## Extract sample names

get.sample.name<-function(fname) strsplit(basename(fname), "--")[[1]][1]

sample.names1<-unname(sapply(cutFs1, get.sample.name))
sample.names2<-unname(sapply(cutFs2, get.sample.name))
sample.names3<-unname(sapply(cutFs3, get.sample.name))
sample.names4<-unname(sapply(cutFs4, get.sample.name))
sample.names5<-unname(sapply(cutFs5, get.sample.name))

## plot nt scores in reads
plotQualityProfile(cutFs5[1:2]) ## Could be any number of samples (e.g. 1-2, 3-16 etc.)

## Filter and trim the reads
### Assign names to output samples
filtrFs1<-file.path(path1.cut, 'filtered',basename(cutFs1))
filtrFs2<-file.path(path2.cut, 'filtered',basename(cutFs2))
filtrFs3<-file.path(path3.cut, 'filtered',basename(cutFs3))
filtrFs4<-file.path(path4.cut, 'filtered',basename(cutFs4))
filtrFs5<-file.path(path5.cut, 'filtered',basename(cutFs5))
## standard filtering parameters
out1<-filterAndTrim(cutFs1, filtrFs1, rev=NULL, filt.rev = NULL,minLen = 50, maxN = 0,maxEE = 2,truncQ = 2,rm.phix = TRUE,compress = TRUE,multithread = FALSE)

out2<-filterAndTrim(cutFs2, filtrFs2, rev=NULL, filt.rev = NULL,minLen = 50, maxN = 0,maxEE = 2,truncQ = 2,rm.phix = TRUE,compress = TRUE,multithread = FALSE)

out3<-filterAndTrim(cutFs3, filtrFs3, rev=NULL, filt.rev = NULL,minLen = 50, maxN = 0,maxEE = 2,truncQ = 2,rm.phix = TRUE,compress = TRUE,multithread = FALSE)

out4<-filterAndTrim(cutFs4, filtrFs4, rev=NULL, filt.rev = NULL,minLen = 50, maxN = 0,maxEE = 2,truncQ = 2,rm.phix = TRUE,compress = TRUE,multithread = FALSE)

out5<-filterAndTrim(cutFs5, filtrFs5, rev=NULL, filt.rev = NULL,minLen = 50, maxN = 0,maxEE = 2,truncQ = 2,rm.phix = TRUE,compress = TRUE,multithread = FALSE)

head(out5) ## Check the output

table(file.exists(filtrFs1))

## D:/1_pool/cutadapt/filtered/p1.NC1_tag_026--NC1_tag_026.fastq has been removed !!! So it needs to be removed here:

filtrFs1_no26<-filtrFs1[filtrFs1!='D:/1_pool/cutadapt/filtered/p1.NC1_tag_026--NC1_tag_026.fastq']
table(file.exists(filtrFs_no26))

### worked?

## Learn error rates
errF1<-learnErrors(filtrFs1_no26, multithread = TRUE)
errF2<-learnErrors(filtrFs2, multithread = TRUE)
errF3<-learnErrors(filtrFs3, multithread = TRUE)
errF4<-learnErrors(filtrFs4, multithread = TRUE)
errF5<-learnErrors(filtrFs5, multithread = TRUE)

## Check out some plots
plotErrors(errF5, nominalQ=T)

## Dereplicating fastqs
derepFs1<- derepFastq(filtrFs1_no26, verbose=T)
derepFs2<- derepFastq(filtrFs2, verbose=T)
derepFs3<- derepFastq(filtrFs3, verbose=T)
derepFs4<- derepFastq(filtrFs4, verbose=T)
derepFs5<- derepFastq(filtrFs5, verbose=T)

## removing 26th sample that was empty from pool1
sample.names1<-sample.names1[sample.names1!='p1.NC1_tag_026']

## Name the derep-class objects by sample names
names(derepFs1)<-sample.names1
names(derepFs2)<-sample.names2
names(derepFs3)<-sample.names3
names(derepFs4)<-sample.names4
names(derepFs5)<-sample.names5

## Sample inference and applying it to derep data
dadaF1<-dada(derepFs1, err=errF1)
dadaF2<-dada(derepFs2, err=errF2)
dadaF3<-dada(derepFs3, err=errF3)
dadaF4<-dada(derepFs4, err=errF4)
dadaF5<-dada(derepFs5, err=errF5)

## Skip merging paired end-reads since it's pacbio...

## Constructing Sequence Table
seqtab1<-makeSequenceTable(dadaF1)
seqtab2<-makeSequenceTable(dadaF2)
seqtab3<-makeSequenceTable(dadaF3)
seqtab4<-makeSequenceTable(dadaF4)
seqtab5<-makeSequenceTable(dadaF5)

## Remove Chimeras
seqtab.nochim1<- removeBimeraDenovo(seqtab1,method='consensus', verbose=T)
seqtab.nochim2<- removeBimeraDenovo(seqtab2,method='consensus', verbose=T)
seqtab.nochim3<- removeBimeraDenovo(seqtab3,method='consensus', verbose=T)
seqtab.nochim4<- removeBimeraDenovo(seqtab4,method='consensus', verbose=T)
seqtab.nochim5<- removeBimeraDenovo(seqtab5,method='consensus', verbose=T)

# export (pools1-5) ASV table:
write.table(seqtab.nochim5, "pool5_frequencies", sep = '\t')
write.table(seqtab.nochim4, "pool4_frequencies", sep = '\t')
write.table(seqtab.nochim3, "pool3_frequencies", sep = '\t')
write.table(seqtab.nochim2, "pool2_frequencies", sep = '\t')
write.table(seqtab.nochim1, "pool1_frequencies", sep = '\t')

## Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab.nochim1)))
table(nchar(getSequences(seqtab.nochim2)))
table(nchar(getSequences(seqtab.nochim3)))
table(nchar(getSequences(seqtab.nochim4)))
table(nchar(getSequences(seqtab.nochim5)))

## Track reads through the pipeline
getN<-function(x) sum(getUniques(x))
track1<-cbind(out1, sapply(dadaF1, getN),rowSums(seqtab.nochim1))
track2<-cbind(out2, sapply(dadaF2, getN),rowSums(seqtab.nochim2))
track3<-cbind(out3, sapply(dadaF3, getN),rowSums(seqtab.nochim3))
track4<-cbind(out4, sapply(dadaF4, getN),rowSums(seqtab.nochim4))
track5<-cbind(out5, sapply(dadaF5, getN),rowSums(seqtab.nochim5))

colnames(track1)<-c('input','filtered','denoisedF','nochim')
rownames(track1) <- sample.names1

colnames(track2)<-c('input','filtered','denoisedF','nochim')
rownames(track2) <- sample.names2

colnames(track3)<-c('input','filtered','denoisedF','nochim')
rownames(track3) <- sample.names3

colnames(track4)<-c('input','filtered','denoisedF','nochim')
rownames(track4) <- sample.names4

colnames(track5)<-c('input','filtered','denoisedF','nochim')
rownames(track5) <- sample.names5

## Assign taxonomy

ref_seq<-'D:/dada2.fasta'

taxa1<-assignTaxonomy(seqtab.nochim1, ref_seq, tryRC = T)
taxa2<-assignTaxonomy(seqtab.nochim2, ref_seq, tryRC = T)
taxa3<-assignTaxonomy(seqtab.nochim3, ref_seq, tryRC = T)
taxa4<-assignTaxonomy(seqtab.nochim4, ref_seq, tryRC = T)
taxa5<-assignTaxonomy(seqtab.nochim5, ref_seq, tryRC = T)
```
```{r}
## Inspect taxonomy
taxa1.print<- taxa1
taxa2.print<- taxa2
taxa3.print<- taxa3
taxa4.print<- taxa4
taxa5.print<- taxa5

## Print the tables out
write.table(x=data.frame(taxa1.print), file="pool1_species",row.names = TRUE)
write.table(x=data.frame(taxa2.print), file="pool2_species", row.names = TRUE)
write.table(x=data.frame(taxa3.print), file="pool3_species",row.names = T)
write.table(x=data.frame(taxa4.print), file="pool4_species",row.names = T)
write.table(x=data.frame(taxa5.print), file="pool5_species",row.names = T)
```
## In Excel
```shell
Remove singleton reads (i.e. ASVs that only appear once! in pool5 that's only 1 ASV)
Remove ASVs (per sample) that comprise less than 0.5% of total reads (per sample)
Remove samples wherein total number of reads below 100
```
