## POOL1
```shell
lima ../raw_data/ps_365_001.fastq.gz ../tags_primers/8tags.fasta p1.fastq --split-named --min-score 80 --same # demultiplexing pool 1
```
```shell
for i in ../p1.*.fastq
  do
  base=$(basename $i .fastq)
  lima $i ../../tags_primers/primers.fasta ${base}_dep.fastq --different --min-score 80 # Removing primers from pool 1
  done
  ```
