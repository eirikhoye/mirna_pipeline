#!bin/bash/

cd /home/jcdenton/
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/hsa/ --species-abbrev hsa --species-verbosename homo_sapiens --mirna-seqs ../../mirna_reference/hsa-pre.fas
