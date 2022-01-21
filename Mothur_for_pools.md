## Trying to assign taxonomy to sequenced parasite pool from sheep with Mothur

### According to nemabiome.ca:
```shell
#Not needed. make.contigs(file=nemabiome.files, processors=4)
1. screen.seqs(fasta=p5.fasta, minlength=200, maxlength=450, maxambig=0, group=p5.groups, processors=2)

2.align.seqs(candidate=p5.good.fasta, template=nematode1_3.fasta) ## 120585 of your sequences generated alignments that eliminated too many bases, a list is provided in D:\Mothur_pipeline\p5.good.flip.accnos.

3.screen.seqs(fasta=p5.good.align, alignreport=p5.good.align.report, minsim=90, minscore=10, group=p5.good.groups) # removed 12537 seqs

4.classify.seqs(fasta=p5.good.good.align, template=nematode1_3.fasta, taxonomy=nematode1_3.tax, method=knn, processors=2, numwanted=3)

5. summary.tax(taxonomy=p5.good.good.nematode1_3.knn.taxonomy, group=p5.good.good.groups) ## Warning! Should not use '--' these in the group names

6.split.groups(fasta=p5.good.fasta, group=p5.good.groups)

7. system(mv p5.good.good.nematode1_3.knn.tax.summary p5_results.summary) ## did it manually (i.e. the name changing)

NOTES:
.groups file should have the first column ending in a unique identifier (i.e. 525634, 525635 etc.) and not contain symbols '-', ':'...
```
### Optional commands to facilitate format conversions
``` shell
sed -n '1~4s/^@/>/p;2~4p' in.fastq > out.fasta ## convert .fastq to .fasta
THEN CONCATENATE ALL .FASTA FILES into 1

sed -n '1~2p' file > file.out ## print every 2nd line starting from 1st
Then generate headers and add a column with tag pairs to generate the .groups file
```
### Handling many .fasta files in R
``` R
path<-''
files<- list.files(path=path,pattern='.fasta')
## Import many files from the same dir
### Create data frames with a second column containing sample names(cut)
for (f in 1:length(files))
{
  file_name<-str_sub(string=files[f], start=58, end=-15)
  file_df<-read.delim(files[f], header = F)
  file_df$V2<- file_name
  assign(x=file_name, value=file_df,envir=.GlobalEnv)
  
}


## make a huge list of all data.frames in the global environment space
l.df<-lapply(ls(), function(x) if (class(get(x))=='data.frame') get(x))
## alternatively l.df<-lapply(ls(pattern="df[0-9+]"), function(x) get(x))

## bind rows in a list
p5_groups<-bind_rows(l.df)

## write it out
write.table(p5_groups, file="p5_groups_file.GROUPS", col.names = F, row.names = F, sep='\t', quote = F)
```
### Assessing Mothur's output in Excel
#### Criteria
```shell
Remove singleton reads for each species/genus 
Remove ASVs (per sample) that comprise less than 0.5% of total reads (per sample)
Remove samples wherein total number of reads below 100
```
