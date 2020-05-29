#!/bin/bash/

> data/bamfiles/mirgenedb.log
for filepath in data/fasta/*.fas.gz
do
  filename=${filepath##*/}
  echo $filename
  echo $filename >> data/bamfiles/mirgenedb.log
  echo 'gunzip'
  gunzip -c $filepath | fastx_uncollapser > $filepath.uncollapsed.fasta

  echo 'bowtie align'
  bowtie -f --norc -l 18 -n 0 -p 4 data/mirna_index/hsa-hg38-pri-30-30_no_multiple_versions \
  $filepath.uncollapsed.fasta -S > data/bamfiles/$filename.sam 2>> data/bamfiles/mirgenedb.log
  echo 'samtools .sam to .bam'
  samtools view -bS data/bamfiles/$filename.sam > data/bamfiles/$filename.bam
  echo 'remove uncollapsed fasta and .sam files'
  rm $filepath.uncollapsed.fasta
  rm data/bamfiles/$filename.sam

done
