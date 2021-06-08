rm config
python scripts/create_config.py $1
echo "run miRTrace"
bash scripts/run_mirtrace.sh config
echo "map to MirGeneDB with bowtie 1.2"
bash scripts/bowtie_mirgenedb_human_fastafiles.sh
echo "count miRNAs"
Rscript scripts/featureCounts.R
