## Trying to assign taxonomy to sequenced parasite pool from sheep with Mothur

### According to nemabiome.ca:
```shell
make.contigs(file=nemabiome.files, processors=4) ### skip this, since it's pacbio reads. Use cutadapt files from the previous DADA2 pipeline


screen.seqs(fasta=nemabiome.trim.contigs.fasta, minlength=200, maxlength=450, maxambig=0, group=nemabiome.contigs.groups, processors=2) ### the same as filterandtrim in DADA2\


align.seqs(candidate=nemabiome.trim.contigs.good.fasta, template=mothur.fasta) ### No dereplication, chimera removal prior to this??\


screen.seqs(fasta=nemabiome.trim.contigs.good.align, alignreport=nemabiome.trim.contigs.good.align.report, minsim=90, minscore=10, group=nemabiome.contigs.good.groups) ## min similarity of 90% and score of 10?\


classify.seqs(fasta=nemabiome.trim.contigs.good.good.align, template=mothur.fasta, taxonomy=mothur.tax, method=knn, processors=2, numwanted=3)\


summary.tax(taxonomy=nemabiome.trim.contigs.good.good.mothur.tax, group=nemabiome.contigs.good.good.groups)\


split.groups(fasta=nemabiome.trim.contigs.good.fasta, group=nemabiome.contigs.good.groups)\


system(mv nemabiome.trim.contigs.good.good.nematode_taxonomy_1_3.knn.tax.summary nemabiome_results.summary)\
```
``` shell
sed -n '1~4s/^@/>/p;2~4p' in.fastq > out.fasta ## convert .fastq to .fasta


sed -n '1~2p' file > file.out ## print every 2nd line starting from 1st
```
``` R
path<-''
files<- list.files(path=path,pattern='.fasta')
## Import many files from the same dir
## create data frames with a second column containing sample names(cut)
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
