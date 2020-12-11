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
samtools sort --threads 8 > $new_dir/${base}R.bam
samtools index $new_dir/${base}R.bam

done
```
### Picard_tools and duplicate removal
```shell
for sample in $new_dir/*R.bam
do
base=$(basename $sample .bam)
java -jar $PICARD_HOME/picard.jar MarkDuplicates I=$sample O=$base.cleanreads.bam REMOVE_DUPLICATES=true M=$base_metrics.txt
samtools view -b -f 2 $base.cleanreads.bam |
samtools sort --threads 8 -T temp > $base.final.sorted.bam
samtools index $base.final.sorted.bam
done
```
### Indexing reference genome
``` shell
samtools faidx /domus/h1/pauliusb/Haemonchus_2018_genome/haemonchusnewest.fa
```
### Making mpileup
``` shell
for sample in $new_dir/*.bam
do
base=$(basename $sample .bam)
samtools mpileup -d 500 -f /domus/h1/pauliusb/Haemonchus_2018_genome/haemonchusnewest.fa \
$sample > ${base}.raw.mpileup
done
```
### Substituting hidden tabulations
``` shell
for sample in $new_dir/*.mpileup
do
sed 's/\t\t/\t!\t!/g' $sample > $sample.new
done
```
