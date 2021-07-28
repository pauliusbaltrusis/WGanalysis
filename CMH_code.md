# Alignment, sorting and duplicate removal
The first step is aligning, sorting, indexing and removing duplicates from our alignment data
## Code
### Alignment, sorting and indexing
``` shell
#!/bin/bash -l
module load bioinfo-tools
module load bwa
module load samtools
module load picard
rec=/home/pauliusb/Cleandata/All_reads
bwa_db=/home/pauliusb/Haemonchus_2018_genome/BWA_all_genomes/haemonchus_cc
new_dir=/home/pauliusb/snic2020-16-116/alignment/WORKING_FOLDER
for sample in $rec/*R1.fq.gz
do
base=$(basename $sample R1.fq.gz)
echo "$base"
bwa mem -t 16 $bwa_db $rec/${base}R1.fq.gz $rec/${base}R2.fq.gz |
samtools view -b |
samtools sort --threads 8 -T temp > $new_dir/${base}R.bam
samtools index $new_dir/${base}R.bam

done
```
### Picard_tools and duplicate removal
```shell
for sample in $new_dir/*.bam
do
base=$(basename $sample .bam)
java -jar $PICARD_ROOT/picard.jar MarkDuplicates -I $sample -O $base.cleanreads.bam -REMOVE_DUPLICATES true -M .txt
samtools view -b -f 2 $base.cleanreads.bam > $base.final.bam
samtools sort $base.final.bam -o $base.final.sorted.bam
done
```
### Indexing reference genome
``` shell
samtools faidx /domus/h1/pauliusb/Haemonchus_2018_genome/haemonchusnewest.fa
```
### Making a joint mpileup
``` shell
samtools mpileup -d 500 --min-MQ 30 --min-BQ 30 --adjust-MQ 50 -f /domus/h1/pauliusb/Haemonchus_2018_genome/haemonchusnewest.fa \
I1_R.final.sorted.bam \
I2_R.final.sorted.bam \
I3_R.final.sorted.bam \
I4_R.final.sorted.bam \
P1_R.final.sorted.bam \
P2_R.final.sorted.bam \
P3_R.final.sorted.bam \
P4_R.final.sorted.bam \
-o $new_dir/IP.unchecked.mpileup
```
### Substituting hidden tabulations
``` shell
sed 's/\t\t/\t!\t!/g' $new_dir/IP.unchecked.mpileup > IP.final.mpileup
```
### Creating synchronized files
``` shell
module load popoolation2
java -ea -Xmx7g -jar /sw/bioinfo/popoolation2/1201/rackham/mpileup2sync.jar --input IP.final.mpileup --min-qual 20 --output IP_filtered.sync
```
### Removing indels
``` shell 
perl /sw/bioinfo/popoolation2/1201/rackham/indel_filtering/identify-indel-regions.pl --min-count 2 --indel-window 5 --input IP.final.mpileup --output indelsCMH.gtf
perl /sw/bioinfo/popoolation2/1201/rackham/indel_filtering/filter-sync-by-gtf.pl --input IP_filtered.sync --output IP_noindelscmh.sync --gtf indelsCMH.gtf
```
### CMH-test 
#### populations 1,2,3,4 are replicates (of untreated group), populations 5,6,7,8 are replicates (of the treated group)
``` shell
perl /sw/bioinfo/popoolation2/1201/rackham/cmh-test.pl --input IP_noindelscmh.sync --output IP_cmh2.cmh --min-count 4 --min-coverage 20 --max-coverage 2% --population 1-5,2-6,3-7,4-8
perl /sw/bioinfo/popoolation2/1201/rackham/cmh2gwas.pl --input IP_cmh2.cmh --output IP_cmh.gwas --min-pvalue 1.0e-20
```
### Displaying the data in Rstudio
``` shell
CMH<- read_tsv("L:/CMH-test/IP_cmh.gwas")

CMH_cleaner<-CMH %>%
  mutate(CHR=replace(CHR, CHR=="hcontortus_chr4_Celeg_TT_arrow_pilon", "4")) %>%
  mutate(CHR=replace(CHR, CHR=="hcontortus_chr5_Celeg_TT_arrow_pilon", "5")) %>%
  mutate(CHR=replace(CHR, CHR=="hcontortus_chrX_Celeg_TT_arrow_pilon", "6")) %>%
  mutate(CHR=replace(CHR, CHR=="hcontortus_chr_mtDNA_arrow_pilon", "7")) %>%
  mutate(CHR=replace(CHR, CHR=="hcontortus_chr1_Celeg_TT_arrow_pilon", "1")) %>%
  mutate(CHR=replace(CHR, CHR=="hcontortus_chr2_Celeg_TT_arrow_pilon", "2")) %>%
  mutate(CHR=replace(CHR, CHR=="hcontortus_chr3_Celeg_TT_arrow_pilon", "3"))

CMH_cleaner<- CMH_cleaner[order(CMH_cleaner$CHR),]
CMH_cleaner<-dplyr::filter(CMH_cleaner, CHR<7)

CMH.plot<-CMH_cleaner %>%
  ggplot(aes(x=BP, y=-log10(P), colour=CHR))+
  xlab(label="Position in chromosome (bp)")+
  ggtitle(label="")+
  ylab(label="")+
  geom_point(aes(), size=1,alpha=0.5)+
  geom_line(aes(x=BP, y=Bonf),color="firebrick", show.legend=NA)+
  scale_color_manual(name = "Chromosome",
values = c( "1" = "#F8766D", "2" = "#C49A00", "3" = "#53B400", "4"="#00C094", "5"="#00B6EB", "6"="#A58AFF", "7"="#FB61D7"),
labels = c(),
breaks = c())+
  theme_bw()+
  #scale_x_continuous(expand=c(0,0), limits=c(0,NA))+
  theme(strip.text.y = element_blank())+
  facet_grid(CHR~.)

CMH.plot
```
### Extracting mean P-value +3SD/5SD
``` shell
threeSDCMH<-dplyr::filter(CMH_cleaner, -log10(P)>=mean(-log10(CMH_cleaner$P))+3*sd(-log10(CMH_cleaner$P)))
threeSDCMHsplit<-split(threeSDCMH,threeSDCMH$CHR)
threeSDCMHsplit1<-threeSDCMHsplit$`1`
threeSDCMHsplit2<-threeSDCMHsplit$`2`
threeSDCMHsplit3<-threeSDCMHsplit$`3`
threeSDCMHsplit4<-threeSDCMHsplit$`4`
threeSDCMHsplit5<-threeSDCMHsplit$`5`
threeSDCMHsplit6<-threeSDCMHsplit$`6`

fiveSDCMH<-dplyr::filter(CMH_cleaner, -log10(P)>=mean(-log10(CMH_cleaner$P))+5*sd(-log10(CMH_cleaner$P)))
fiveSDCMHsplit<-split(fiveSDCMH,fiveSDCMH$CHR)
fiveSDCMHsplit1<-fiveSDCMHsplit$`1`
fiveSDCMHsplit2<-fiveSDCMHsplit$`2`
fiveSDCMHsplit3<-fiveSDCMHsplit$`3`
fiveSDCMHsplit4<-fiveSDCMHsplit$`4`
fiveSDCMHsplit5<-fiveSDCMHsplit$`5`
fiveSDCMHsplit6<-fiveSDCMHsplit$`6`
```
### Overlapping annotated SNPs with values >3SD/5SD obtained in CMH test
``` shell
cmh_annot<-read_tsv("L:/CMH-test/subsetted_datasets/Jul22.subsetted.final.CMH.vcf", col_names = FALSE)
colnames(cmh_annot) <- c("Chromosome", "BP", "X3", "X4", "X5")

cmh_annot_split<-split(cmh_annot, cmh_annot$Chromosome)
cmh_annot_split1<-cmh_annot_split$hcontortus_chr1_Celeg_TT_arrow_pilon
cmh_annot_split2<-cmh_annot_split$hcontortus_chr2_Celeg_TT_arrow_pilon
cmh_annot_split3<-cmh_annot_split$hcontortus_chr3_Celeg_TT_arrow_pilon
cmh_annot_split4<-cmh_annot_split$hcontortus_chr4_Celeg_TT_arrow_pilon
cmh_annot_split5<-cmh_annot_split$hcontortus_chr5_Celeg_TT_arrow_pilon
cmh_annot_split6<-cmh_annot_split$hcontortus_chrX_Celeg_TT_arrow_pilon

cmh_ann_3sd1<-dplyr::inner_join(threeSDCMHsplit1, cmh_annot_split1, by='BP')
cmh_ann_3sd2<-dplyr::inner_join(threeSDCMHsplit2, cmh_annot_split2, by='BP')
cmh_ann_3sd3<-dplyr::inner_join(threeSDCMHsplit3, cmh_annot_split3, by='BP')
cmh_ann_3sd4<-dplyr::inner_join(threeSDCMHsplit4, cmh_annot_split4, by='BP')
cmh_ann_3sd5<-dplyr::inner_join(threeSDCMHsplit5, cmh_annot_split5, by='BP')
cmh_ann_3sd6<-dplyr::inner_join(threeSDCMHsplit6, cmh_annot_split6, by='BP')

cmh_ann_5sd1<-dplyr::inner_join(fiveSDCMHsplit1, cmh_annot_split1, by='BP')
cmh_ann_5sd2<-dplyr::inner_join(fiveSDCMHsplit2, cmh_annot_split2, by='BP')
cmh_ann_5sd3<-dplyr::inner_join(fiveSDCMHsplit3, cmh_annot_split3, by='BP')
cmh_ann_5sd4<-dplyr::inner_join(fiveSDCMHsplit4, cmh_annot_split4, by='BP')
cmh_ann_5sd5<-dplyr::inner_join(fiveSDCMHsplit5, cmh_annot_split5, by='BP')
cmh_ann_5sd6<-dplyr::inner_join(fiveSDCMHsplit6, cmh_annot_split6, by='BP')

cmh_ann_3SD<-dplyr::bind_rows(cmh_ann_3sd1,cmh_ann_3sd2,cmh_ann_3sd3,cmh_ann_3sd4,cmh_ann_3sd5,cmh_ann_3sd6)

cmh_ann_5SD<-dplyr::bind_rows(cmh_ann_5sd1,cmh_ann_5sd2,cmh_ann_5sd3,cmh_ann_5sd4,cmh_ann_5sd5,cmh_ann_5sd6)
```
### Plotting annotated snps above 3/5 SDs in CMH-test
``` shell
cmh_35SD<- cmh_ann_3SD %>%
  ggplot(aes())+
  geom_point(aes(x=BP, y=-log10(P)), color='gray')+
  geom_point(data=cmh_ann_5SD, aes(x=BP, y=-log10(P)),color='black')+
  labs(y="SNP frequency difference, -log10(P)", x="Position in chromosome (bp)")+
  facet_grid(CHR~.)+
  theme_bw()
cmh_35SD
```
## FET vs CMH plot
### Data prep
``` shell
colnames(fet_final)[3] <- "BP"
# It is critical to split the data into chromosomes and then join by 'BP', otherwise incorrect data points are joined!
fet_split<-split(fet_final, fet_final$Chromosome)
CMH_split<-split(CMH_cleaner, CMH_cleaner$CHR)

fet_splitchr1<-fet_split$`1`
fet_splitchr2<-fet_split$`2`
fet_splitchr3<-fet_split$`3`
fet_splitchr4<-fet_split$`4`
fet_splitchr5<-fet_split$`5`
fet_splitchr6<-fet_split$`6`

CMH_splitchr1<-CMH_split$`1`
CMH_splitchr2<-CMH_split$`2`
CMH_splitchr3<-CMH_split$`3`
CMH_splitchr4<-CMH_split$`4`
CMH_splitchr5<-CMH_split$`5`
CMH_splitchr6<-CMH_split$`6`

fet_CMH_split1<-dplyr::inner_join(fet_splitchr1,CMH_splitchr1, by="BP")
fet_CMH_split2<-dplyr::inner_join(fet_splitchr2,CMH_splitchr2, by="BP")
fet_CMH_split3<-dplyr::inner_join(fet_splitchr3,CMH_splitchr3, by="BP")
fet_CMH_split4<-dplyr::inner_join(fet_splitchr4,CMH_splitchr4, by="BP")
fet_CMH_split5<-dplyr::inner_join(fet_splitchr5,CMH_splitchr5, by="BP")
fet_CMH_split6<-dplyr::inner_join(fet_splitchr6,CMH_splitchr6, by="BP")

fet_cmh<-dplyr::bind_rows(fet_CMH_split1,fet_CMH_split2,fet_CMH_split3,fet_CMH_split4,fet_CMH_split5,fet_CMH_split6)
```
### Actual plotting
``` shell
# Note! 1388262 - is the number of SNPs shared between the two tests
## Otherwise it could be 0.05/number of total snps
fet_cmh_plot<- fet_cmh_ann %>% 
  ggplot(aes())+
  xlab(label='FET -log10(P)-values')+
  ylab(label='CMH -log10(P)-values')+
  geom_hline(yintercept =-log10(0.05/1388262) , color='firebrick', show.legend = NA)+
  geom_vline(xintercept =-log10(0.05/1388262), color='firebrick', show.legend = NA)+
  geom_point(aes(y=-log10(P.y), x=P.x), color='gray')+
  geom_point(data=ann_points, aes(y=-log10(P.y), x=P.x), color='black')+
  theme_bw()+
  scale_fill_discrete(name="Putative SNP impact", labels=c("Moderate/High", "Unknown"))
fet_cmh_plot
```
### Testing correlation between the FET and CMH shared points
``` shell
cor.test(fet_cmh$P.x, -log10(fet_cmh$P.y), method='pearson', alternative='two.sided')
```
### Supplementing the CMH vs FET plot with annotated SNPs
``` shell
modhigh<-read_tsv("L:/CMH-test/subsetted_datasets/Jul22.subsetted.final.FET.vcf",col_names = FALSE)
colnames(modhigh) <- c("Chromosome", "BP", "X3", "X4", "X5")
modhigh$marker<-c('M')

modhigh_split<-split(modhigh, modhigh$Chromosome)
modhigh_splitchr1<-modhigh_split$hcontortus_chr1_Celeg_TT_arrow_pilon
modhigh_splitchr2<-modhigh_split$hcontortus_chr2_Celeg_TT_arrow_pilon
modhigh_splitchr3<-modhigh_split$hcontortus_chr3_Celeg_TT_arrow_pilon
modhigh_splitchr4<-modhigh_split$hcontortus_chr4_Celeg_TT_arrow_pilon
modhigh_splitchr5<-modhigh_split$hcontortus_chr5_Celeg_TT_arrow_pilon
modhigh_splitchrX<-modhigh_split$hcontortus_chrX_Celeg_TT_arrow_pilon
# Again splitting and joining annotated SNPs with fet_cmh by base positions per every chrom
fet_cmh_split<-split(fet_cmh, fet_cmh$Chromosome)
fet_cmhchr1<-fet_cmh_split$`1`
fet_cmhchr2<-fet_cmh_split$`2`
fet_cmhchr3<-fet_cmh_split$`3`
fet_cmhchr4<-fet_cmh_split$`4`
fet_cmhchr5<-fet_cmh_split$`5`
fet_cmhchr6<-fet_cmh_split$`6`

fet_CMH_ann1<-dplyr::left_join(fet_cmhchr1,modhigh_splitchr1, by="BP")
fet_CMH_ann2<-dplyr::left_join(fet_cmhchr2,modhigh_splitchr2, by="BP")
fet_CMH_ann3<-dplyr::left_join(fet_cmhchr3,modhigh_splitchr3, by="BP")
fet_CMH_ann4<-dplyr::left_join(fet_cmhchr4,modhigh_splitchr4, by="BP")
fet_CMH_ann5<-dplyr::left_join(fet_cmhchr5,modhigh_splitchr5, by="BP")
fet_CMH_ann6<-dplyr::left_join(fet_cmhchr6,modhigh_splitchrX, by="BP")

fet_cmh_ann<-dplyr::bind_rows(fet_CMH_ann1,fet_CMH_ann2,fet_CMH_ann3,fet_CMH_ann4,fet_CMH_ann5,fet_CMH_ann6)
# making a data frame with SNPs that were identified by annotation
ann_points<-dplyr::filter(fet_cmh_ann, marker=="Moderate/High")
```
