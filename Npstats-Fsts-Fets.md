
# Nucleotide diversity, Taj D / Fst & FET analysis
## Merging separate aligned .bam files into two "I" and "P"
``` shell
module load bioinfo-tools
module load samtools

samtools merge merged.I.bam -b list_bam_I.txt
samtools merge merged.P.bam -b list_bam_P.txt
```
## Generating mpileups
``` shell
samtools mpileup -d 500 --min-MQ 30 --min-BQ 30 --adjust-MQ 50 -f /domus/h1/pauliusb/Haemonchus_2018_genome/haemonchusnewest.fa \
merged.I.bam \
merged.P.bam \
-o IP2.unchecked.mpileup
sed 's/\t\t/\t!\t!/g' IP2.unchecked.mpileup > IP2.final.mpileup
```
## Spliting mpileups by sample (I or P) and chr
``` shell
less IP2.final.mpileup | cut -f 1,2,3,4,5,6 > I.only.mpileup
less IP2.final.mpileup | cut -f 1,2,3,7,8,9 > P.only.mpileup
```
### By chr for I
``` shell
grep "hcontortus_chr4_Celeg_TT_arrow_pilon" I.only.mpileup > I_CHR4.mpileup
grep "hcontortus_chr5_Celeg_TT_arrow_pilon" I.only.mpileup > I_CHR5.mpileup
grep "hcontortus_chr1_Celeg_TT_arrow_pilon" I.only.mpileup > I_CHR1.mpileup
grep "hcontortus_chr2_Celeg_TT_arrow_pilon" I.only.mpileup > I_CHR2.mpileup
grep "hcontortus_chr3_Celeg_TT_arrow_pilon" I.only.mpileup > I_CHR3.mpileup
grep "hcontortus_chrX_Celeg_TT_arrow_pilon" I.only.mpileup > I_CHRX.mpileup
grep "hcontortus_chr_mtDNA_arrow_pilon" I.only.mpileup > I_CHRMT.mpileup
```
### By chr for P
``` shell
grep "hcontortus_chr4_Celeg_TT_arrow_pilon" P.only.mpileup > P_CHR4.mpileup
grep "hcontortus_chr5_Celeg_TT_arrow_pilon" P.only.mpileup > P_CHR5.mpileup
grep "hcontortus_chr1_Celeg_TT_arrow_pilon" P.only.mpileup > P_CHR1.mpileup
grep "hcontortus_chr2_Celeg_TT_arrow_pilon" P.only.mpileup > P_CHR2.mpileup
grep "hcontortus_chr3_Celeg_TT_arrow_pilon" P.only.mpileup > P_CHR3.mpileup
grep "hcontortus_chrX_Celeg_TT_arrow_pilon" P.only.mpileup > P_CHRX.mpileup
grep "hcontortus_chr_mtDNA_arrow_pilon" P.only.mpileup > P_CHRMT.mpileup
```
## Npstat/1
``` shell
module load NPStat/1
for i in *_CHR*
do
npstat -n 200 -l 100000 -maxcov 500 -minqual 20 $i
done
```
### Tajima's D vs 100kbp windows and nucleotide diversity (Pi) for "I"

### Tajima's D vs 100kbp windows and nucleotide diversity (Pi) for "P"


## Syncing
``` shell
module load bcftools
module load popoolation2
java -ea -Xmx7g -jar /sw/bioinfo/popoolation2/1201/rackham/mpileup2sync.jar --input IP2.final.mpileup --min-qual 20 --output IP2_filtered.sync
```
## Finding +5 bp indels and removing them from the synced files
``` shell
perl /sw/bioinfo/popoolation2/1201/rackham/indel_filtering/identify-indel-regions.pl --min-count 2 --indel-window 5 --input IP2.final.mpileup --output indels.gtf
perl /sw/bioinfo/popoolation2/1201/rackham/indel_filtering/filter-sync-by-gtf.pl --input IP2_filtered.sync --output IP_noindels.sync --gtf indels.gtf
```
## Fst and Fet calculations
``` shell
perl /sw/bioinfo/popoolation2/1201/rackham/fst-sliding.pl --input IP_noindels.sync --output IP2.fst --pool-size 1000 --window-size 10000 --step-size 5000 --min-count 4 --min-coverage 20 --max-coverage 2%
perl /sw/bioinfo/popoolation2/1201/rackham/fisher-test.pl --input IP_noindels.sync --output IP2.fet  --min-count 4 --min-coverage 20 --max-coverage 2% --suppress-noninformative
```
## Genewise Fst calculations
### annotations file filtering/ retaining only whole genes
``` shell
more PRJEB506.WBPS16.gtf | egrep "gene" haemonchus_contortus.PRJEB506.WBPS16.canonical_geneset.gtf.gz > PRJEB506.WBPS16.modified.gtf

```
### annotations file filtering/ retaining only exons of genes
``` shell
# or retaining only exons

module load bioinfo-tools
module load GenomeTools
gt gff3 -tidy -sort -retainids hcontortus_original.gff3 > hcontortus_original.fixed.gff3
gt gff3_to_gtf hcontortus_original.fixed.gff3 > haemonchus.gtf
less haemonchus.gtf | egrep "exon" > haemonchus_edited.gtf
```
### Genewise Fst
``` shell
perl /sw/bioinfo/popoolation2/1201/rackham/create-genewise-sync.pl --input IP2_filtered.sync --gtf PRJEB506.WBPS16.modified.gtf #or haemonchus_edited.gtf# --output IP2.allgeneFst2021.sync

perl /sw/bioinfo/popoolation2/1201/rackham/fst-sliding.pl --min-count 4 --min-coverage 20 \
--max-coverage 2% --pool-size 1000 --min-covered-fraction 0.0 \
--window-size 1000000 --step-size 1000000 --input IP2.allgeneFst2021.sync --output IP2_allgeneFst2021.fst
```
## Exporting into easier-to-handle formats
``` shell
perl /sw/bioinfo/popoolation2/1201/rackham/export/pwc2igv.pl --input IP2.fst --output IP2.fst.igv
perl /sw/bioinfo/popoolation2/1201/rackham/export/pwc2igv.pl --input IP2.fet --output IP2.fet.igv
```
## FST results / Rstudio
### importing, plotting
``` shell
library(tidyverse)
library(ggpubr)
library(patchwork)
# IP2.fst -> IP2_10ksize.fst.edited.txt (renamed the file) 
fst<-read_tsv("L:/NPstat1/IP2_10ksize.fst.edited.txt")

fst_final <-fst %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chr4_Celeg_TT_arrow_pilon", "4")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chr5_Celeg_TT_arrow_pilon", "5")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chrX_Celeg_TT_arrow_pilon", "6")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chr_mtDNA_arrow_pilon", "7")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chr1_Celeg_TT_arrow_pilon", "1")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chr2_Celeg_TT_arrow_pilon", "2")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chr3_Celeg_TT_arrow_pilon", "3"))
fst_final<- fst_final[order(fst_final$Chromosome),]
colnames(fst_final)[5] <- "Fst"

fst_final<-dplyr::filter(fst_final,Chromosome<7)
f5t<-fst_final %>% 
  ggplot(aes(x=window, y=Fst, color=Chromosome))+
  ggtitle('')+
  xlab(label="Genomic windows (10kbp) throughout genome")+
  ylab(expression(Pairwise~genetic~difference~(F[ST])))+
  scale_y_continuous(breaks = seq(0,0.15, by=0.05))+
  scale_x_continuous(breaks=seq(0,5e+07, by=1e+07),labels = c('0','10 Mbp', '20 Mbp', '30 Mbp', '40 Mbp', '50 Mbp'))+
   scale_color_manual(name = "Chromosome",
values = c( "1" = "#F8766D", "2" = "#C49A00", "3" = "#53B400", "4"="#00C094", "5"="#00B6EB", "6"="#A58AFF", "7"="#FB61D7"),
labels = c("1", "2", "3", "4", "5", "X", "mt"),
breaks = c())+
  geom_point(size=1)+
  geom_hline(yintercept = mean(fst_final$Fst)+3*sd(fst_final$Fst), linetype='dashed', color="black")+
  geom_hline(yintercept = mean(fst_final$Fst)+5*sd(fst_final$Fst), linetype='dotted', color='black')+
  geom_point(data=threeabove5sdpoints, aes(x=window,y=Fst), color="red")+
  theme_bw()+
  theme(strip.text.y = element_blank())+
  facet_grid(Chromosome~.)
f5t
```
### Importing SNP/window overlap over 5SD/3SD
``` shell
snps5sdfst<-read_tsv('L:/NPstat1/joined5sd.bed', col_names = FALSE)
colnames(snps5sdfst)<-c("Chromosome", "start", "end", "CHR", "startpos", "window","Fst")
snps5sdfst_final<- snps5sdfst %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chr4_Celeg_TT_arrow_pilon", "4")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chr5_Celeg_TT_arrow_pilon", "5")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chrX_Celeg_TT_arrow_pilon", "6")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chr_mtDNA_arrow_pilon", "7")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chr1_Celeg_TT_arrow_pilon", "1")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chr2_Celeg_TT_arrow_pilon", "2")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chr3_Celeg_TT_arrow_pilon", "3"))
snps5sdfst_final<- snps5sdfst_final[order(snps5sdfst_final$Chromosome),]

snps3sdfst<-read_tsv('L:/NPstat1/joined3sd.bed', col_names = FALSE)
colnames(snps3sdfst)<-c("Chromosome", "start", "end", "CHR", "startpos", "window","Fst")
snps3sdfst_final<- snps3sdfst %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chr4_Celeg_TT_arrow_pilon", "4")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chr5_Celeg_TT_arrow_pilon", "5")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chrX_Celeg_TT_arrow_pilon", "6")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chr_mtDNA_arrow_pilon", "7")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chr1_Celeg_TT_arrow_pilon", "1")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chr2_Celeg_TT_arrow_pilon", "2")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chr3_Celeg_TT_arrow_pilon", "3"))
snps3sdfst_final<- snps3sdfst_final[order(snps3sdfst_final$Chromosome),]
```
### Identifying three consecutive 10kbp windows with fst > 3 or 5SD
``` shell
threesdfst <- mean(fst_final$Fst)+3*sd(fst_final$Fst)
fst_final_3sd<-dplyr::filter(fst_final,Fst>threesdfst)
write_delim(fst_final_3sd, file="L:/NPstat1/fst_final_3sd.txt", delim="\t")
fst_final_3sd$numerical<-c(1:770)

numbers = fst_final_3sd$number
x1<-c(diff(numbers),0)
x2<-c(0,diff(numbers[-1]),0)
x3<-c(0,diff(numbers[c(-1,-2)]),0,0)
colSums(rbind(x1,x2,x3))==3

sum(colSums(rbind(x1,x2,x3))==3)

Fstwindows41<-which(colSums(rbind(x1,x2,x3))==3)
Fstwindows41<-as.data.frame(Fstwindows41)
colnames(Fstwindows41) <- 'numerical'

Fst3sdlocations<-dplyr::semi_join(fst_final_3sd,Fstwindows41, by = "numerical")
Fst3sdlocations

## 5SD

fivesdfst<-mean(fst_final$Fst)+5*sd(fst_final$Fst)
fst_final_5sd<-dplyr::filter(fst_final,Fst>fivesdfst)
write.table(fst_final_5sd)

numbers5sd = fst_final_5sd$number
x1b<-c(diff(numbers5sd),0)
x2b<-c(0,diff(numbers5sd[-1]),0)
x3b<-c(0,diff(numbers5sd[c(-1,-2)]),0,0)
colSums(rbind(x1b,x2b,x3b))==3

sum(colSums(rbind(x1b,x2b,x3b))==3)

which(colSums(rbind(x1b,x2b,x3b))==3)

threeabove5sdpoints<-dplyr::filter(fst_final,number==17315 | number== 17316| number== 17317)
```
## Genewise Fst analysis
```shell
library(ggrepel)
library(tidyverse)
library(ggpubr)

genewise<-read_tsv("L:/NPstat1/IP2_allIVMgenes.fst")

genewiseplot<-genewise %>% 
  ggplot(aes(x=reorder(gene,-Fst), y=log10(Fst)))+
  xlab(label="Genes")+
  ylab(label="Genewise differentiation log10(Fst)")+
  geom_point(aes(), color="gray")+
  geom_point(data=candidate_list,aes(x=reorder(gene, -Fst), y=log10(Fst)), color="black")+
  geom_text_repel(data=candidate_list,aes(label=name), size=3, max.overlaps = Inf)+
  geom_point(data=topgene, color="#00C094")+
  geom_text_repel(data=topgene,aes(label=name),size=3)+
  geom_hline(yintercept = log10(mean(genewise$Fst)+3*sd(genewise$Fst)), linetype='dashed', color="black")+
  geom_hline(yintercept = log10(mean(genewise$Fst)+5*sd(genewise$Fst)), linetype='dotted', color='black')+
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(),panel.background =element_blank(), axis.line = element_line(colour = "black"), panel.border = element_rect(fill=NA))+
  scale_x_discrete(expand = expansion(mult = c(.05,.05)))
genewiseplot
ggsave("Genewise_comparison2.png")

tablegenewisehigherthanpgp1<-dplyr::filter(genewise, Fst>0.04251413)
```
```{r}
candidate_list<-dplyr::filter(genewise, gene==18930 | gene==2927 |gene== 2055 |gene== 7837| gene== 16528 | gene==947
                              | gene==1504 | gene==6560 | gene==16528 | gene==101 | gene==15246 | gene==14130 |
                                gene==16296 | gene==3653 | gene==12402 | gene==2699 | gene==10029 | gene==4237 |
                                gene==1820 | gene==19215 | gene==3663 | gene==15296 | gene==4388 | gene==12903 |
                                gene==13314 | gene==13315 | gene==19292)
candidate_list$gene
candidate_list$name<-c('pgp-1',	'osm-3',	'haf-6',	'pgp-9',	'pgp-9.1',	'lgc-36',	'glc-3',	'osm-6',	'lgc-55',	'avr-15',	'osm-1',	'che-11',	'pgp-12',	'pgp-3',	'ggr-3',	'dyf-7',	'unc-9',	'mrp-1',	'glc-2',	'che-3',	'dyf-11',	'unc-38',	'avr-14',	'osm-5',	'glc-5',	'lgc-37')
candidate_list<-dplyr::arrange(candidate_list, desc(Fst))
write.table(candidate_list, sep="\t", dec=".", row.names = FALSE)
topgene<-dplyr::filter(genewise, Fst>0.15)
topgene$name<-"HCON_00129180"
```
## FET results
### importing, plotting etc.
``` shell
fet<-read_tsv("L:/NPstat1/IP2.fet.igv")

fet_final <-fet %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chr4_Celeg_TT_arrow_pilon", "4")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chr5_Celeg_TT_arrow_pilon", "5")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chrX_Celeg_TT_arrow_pilon", "6")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chr_mtDNA_arrow_pilon", "7")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chr1_Celeg_TT_arrow_pilon", "1")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chr2_Celeg_TT_arrow_pilon", "2")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome=="hcontortus_chr3_Celeg_TT_arrow_pilon", "3"))
fet_final<- fet_final[order(fet_final$Chromosome),]
colnames(fet_final)[5] <- "P"

fet_final<-dplyr::filter(fet_final,Chromosome<7)
F3t<-fet_final %>% 
  ggplot(aes(x=Start, y=P, color=Chromosome))+
  ggtitle('')+
  xlab(label="Position in chromosome (bp)")+
  ylab(label="SNP frequency difference (-log10(P))")+
   scale_color_manual(name = "Chromosome",
values = c( "1" = "#F8766D", "2" = "#C49A00", "3" = "#53B400", "4"="#00C094", "5"="#00B6EB", "6"="#A58AFF", "7"="#FB61D7"),
labels = c(),
breaks = c())+
  geom_point(aes(), size=1,alpha=0.5)+
  #geom_hline(yintercept = mean(fet_final$P)+3*sd(fet_final$P), linetype='dashed', color="black")
  #geom_hline(yintercept = mean(fet_final$P)+5*sd(fet_final$P), linetype='dotted', color='black')
  geom_hline(yintercept=-log10(0.05/2709683), color="firebrick")+
  theme_bw()+
  theme(strip.text.y = element_blank())+
  facet_grid(Chromosome~.)
F3t
```
### Annotating FET SNPs (>3 or >5 SDs) with called known impact having SNPs (from bcftools call)
``` shell
fet_annot<-read_tsv("L:/CMH-test/subsetted_datasets/Jul22.subsetted.final.FET.vcf", col_names = FALSE)
colnames(fet_annot) <- c("Chromosome", "BP", "X3", "X4", "X5")
fet_annot_split<-split(fet_annot, fet_annot$Chromosome)
fet_annot_split1<-fet_annot_split$hcontortus_chr1_Celeg_TT_arrow_pilon
fet_annot_split2<-fet_annot_split$hcontortus_chr2_Celeg_TT_arrow_pilon
fet_annot_split3<-fet_annot_split$hcontortus_chr3_Celeg_TT_arrow_pilon
fet_annot_split4<-fet_annot_split$hcontortus_chr4_Celeg_TT_arrow_pilon
fet_annot_split5<-fet_annot_split$hcontortus_chr5_Celeg_TT_arrow_pilon
fet_annot_split6<-fet_annot_split$hcontortus_chrX_Celeg_TT_arrow_pilon

fet_ann_3sd1<-dplyr::inner_join(threeSDFETsplit1, fet_annot_split1, by='BP')
fet_ann_3sd2<-dplyr::inner_join(threeSDFETsplit2, fet_annot_split2, by='BP')
fet_ann_3sd3<-dplyr::inner_join(threeSDFETsplit3, fet_annot_split3, by='BP')
fet_ann_3sd4<-dplyr::inner_join(threeSDFETsplit4, fet_annot_split4, by='BP')
fet_ann_3sd5<-dplyr::inner_join(threeSDFETsplit5, fet_annot_split5, by='BP')
fet_ann_3sd6<-dplyr::inner_join(threeSDFETsplit6, fet_annot_split6, by='BP')

fet_ann_5sd1<-dplyr::inner_join(fiveSDFETsplit1, fet_annot_split1, by='BP')
fet_ann_5sd2<-dplyr::inner_join(fiveSDFETsplit2, fet_annot_split2, by='BP')
fet_ann_5sd3<-dplyr::inner_join(fiveSDFETsplit3, fet_annot_split3, by='BP')
fet_ann_5sd4<-dplyr::inner_join(fiveSDFETsplit4, fet_annot_split4, by='BP')
fet_ann_5sd5<-dplyr::inner_join(fiveSDFETsplit5, fet_annot_split5, by='BP')
fet_ann_5sd6<-dplyr::inner_join(fiveSDFETsplit6, fet_annot_split6, by='BP')

fet_ann_3SD<-dplyr::bind_rows(fet_ann_3sd1,fet_ann_3sd2,fet_ann_3sd3,fet_ann_3sd4,fet_ann_3sd5,fet_ann_3sd6)

fet_ann_5SD<-dplyr::bind_rows(fet_ann_5sd1,fet_ann_5sd2,fet_ann_5sd3,fet_ann_5sd4,fet_ann_5sd5,fet_ann_5sd6)
```
### Plotting
``` shell
fet_35SD<- fet_ann_3SD %>%
  ggplot(aes())+
  geom_point(aes(x=Start, y=P), color='gray')+
  geom_point(data=fet_ann_5SD, aes(x=Start, y=P),color='black')+
  labs(y="SNP frequency difference, -log10(P)", x="Position in chromosome (bp)")+
  facet_grid(Chromosome.x~.)+
  theme_bw()
fet_35SD
```
## Combining plots with ggarange()
### For supplementary
``` shell
x1<-ggarrange(F3t,CMH.plot,ncol = 2,nrow=1, common.legend = TRUE, labels=c("a", "b"))
x1
ggsave("Supplementary1.FET.Cmh-plot.png")

y0<-ggarrange(fet_35SD, cmh_35SD, ncol=2, nrow=1,labels=c("b", "c"))
y1<-ggarrange(f5t_35sd, y0, labels=c("a"), ncol=1, nrow=2)
y1
ggsave('Supplementary2.F5t,FET.CMH.35SDs.png')
```
### for figure 2.
```shell
x2<-ggarrange(genewiseplot,fet_cmh_plot, ncol=2, nrow=1, labels=c("b","c"))
x2             
ggarrange(f5t,x2, ncol = 1, nrow = 2, common.legend = TRUE,legend = 'right', labels=c("a"))
```
