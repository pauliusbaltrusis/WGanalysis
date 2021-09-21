```shell
lima INPUT BARCODES(.FASTA) OUTPUT --split-named --min-score 80
```
```shell
for i in p1.*.fastq
  do
  base=$(basename $i .fastq)
  lima $i ../tags_primers/primers.fasta ${base}_dep.fastq --different --min-score 80
  done
  ```
