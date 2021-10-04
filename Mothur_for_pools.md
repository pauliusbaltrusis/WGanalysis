## Trying to assign taxonomy to sequenced parasite pool from sheep with Mothur

### According to nemabiome.ca:
```shell
make.contigs(file=nemabiome.files, processors=4) ### skip this, since it's pacbio reads. Use cutadapt files from the previous DADA2 pipeline\
screen.seqs(fasta=nemabiome.trim.contigs.fasta, minlength=200, maxlength=450, maxambig=0, group=nemabiome.contigs.groups, processors=2) ### the same as filterandtrim in DADA2\
align.seqs(candidate=nemabiome.trim.contigs.good.fasta, template=mothur.fasta) ### No dereplication, chimera removal prior to this??\
screen.seqs(fasta=nemabiome.trim.contigs.good.align, alignreport=nemabiome.trim.contigs.good.align.report, minsim=90, minscore=10, group=nemabiome.contigs.good.groups) ## min similarity of \
90% and score of 10?\
classify.seqs(fasta=nemabiome.trim.contigs.good.good.align, template=mothur.fasta, taxonomy=mothur.tax, method=knn, processors=2, numwanted=3)\
summary.tax(taxonomy=nemabiome.trim.contigs.good.good.mothur.tax, group=nemabiome.contigs.good.good.groups)\
split.groups(fasta=nemabiome.trim.contigs.good.fasta, group=nemabiome.contigs.good.groups)\
system(mv nemabiome.trim.contigs.good.good.nematode_taxonomy_1_3.knn.tax.summary nemabiome_results.summary)\
```
