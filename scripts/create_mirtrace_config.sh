#!/bin/bash

# usage: bash scripts/create_mirtrace_config.sh <directory_of_fastqs> <study_name> <adapter_sequence>

echo "" > data/"$2"_config

for line in "$1"/*.fastq.gz
do
  filename=$line
  sample_name="${line##*/}"
  sample_name=${sample_name%%_*}
  echo $filename,$sample_name,$3 >> data/"$2"_config
done
