#!/bin/bash

rm /home/jcdenton/projects/mirna_pipeline/custom_databases/*
cd /home/jcdenton/projects/mirna_pipeline/scripts/generate_custom_databases/

python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev aae --species-verbosename yellow_fever_mosquito --mirna-seqs ../../mirna_reference/aae-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev aca --species-verbosename green_anole_lizard --mirna-seqs ../../mirna_reference/aca-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev ami --species-verbosename american_alligator --mirna-seqs ../../mirna_reference/ami-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev asu --species-verbosename large_roundworm --mirna-seqs ../../mirna_reference/asu-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev bfl --species-verbosename florida_lancelet --mirna-seqs ../../mirna_reference/bfl-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev bge --species-verbosename cockroach --mirna-seqs ../../mirna_reference/bge-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev bta --species-verbosename cow --mirna-seqs ../../mirna_reference/bta-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev cbr --species-verbosename roundworm2 --mirna-seqs ../../mirna_reference/cbr-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev cel --species-verbosename roundworm1 --mirna-seqs ../../mirna_reference/cel-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev cfa --species-verbosename dog --mirna-seqs ../../mirna_reference/cfa-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev cgi --species-verbosename pacific_oyster --mirna-seqs ../../mirna_reference/cgi-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev cin --species-verbosename sea_squirt --mirna-seqs ../../mirna_reference/cin-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev cli --species-verbosename rock_pigeon --mirna-seqs ../../mirna_reference/cli-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev cpi --species-verbosename western_painted_turtle --mirna-seqs ../../mirna_reference/cpi-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev cpo --species-verbosename guinea_pig --mirna-seqs ../../mirna_reference/cpo-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev cte --species-verbosename polychaete_worm --mirna-seqs ../../mirna_reference/cte-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev dan --species-verbosename fruit_fly2 --mirna-seqs ../../mirna_reference/dan-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev dme --species-verbosename fruit_fly1 --mirna-seqs ../../mirna_reference/dme-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev dmo --species-verbosename fruit_fly3 --mirna-seqs ../../mirna_reference/dmo-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev dno --species-verbosename nine-banded_armadillo --mirna-seqs ../../mirna_reference/dno-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev dpu --species-verbosename common_water_flea --mirna-seqs ../../mirna_reference/dpu-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev dre --species-verbosename zebrafish --mirna-seqs ../../mirna_reference/dre-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev efe --species-verbosename common_brandling_worm --mirna-seqs ../../mirna_reference/efe-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev ete --species-verbosename lesser_hedgehog_tenrec --mirna-seqs ../../mirna_reference/ete-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev gga --species-verbosename chicken --mirna-seqs ../../mirna_reference/gga-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev hme --species-verbosename longwing_butterfly --mirna-seqs ../../mirna_reference/hme-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev hsa --species-verbosename homo_sapiens --mirna-seqs ../../mirna_reference/hsa-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev isc --species-verbosename deer_tick --mirna-seqs ../../mirna_reference/isc-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev lan --species-verbosename lingula --mirna-seqs ../../mirna_reference/lan-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev lgi --species-verbosename owl_limpet --mirna-seqs ../../mirna_reference/lgi-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev mdo --species-verbosename gray_short-tailed_opossum --mirna-seqs ../../mirna_reference/mdo-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev mml --species-verbosename rhesus_monkey --mirna-seqs ../../mirna_reference/mml-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev mmu --species-verbosename mus_musculus --mirna-seqs ../../mirna_reference/mmu-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev oan --species-verbosename platypus --mirna-seqs ../../mirna_reference/oan-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev ocu --species-verbosename rabbit --mirna-seqs ../../mirna_reference/ocu-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev pfl --species-verbosename ptychodera --mirna-seqs ../../mirna_reference/pfl-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev pmi --species-verbosename bat_starfish --mirna-seqs ../../mirna_reference/pmi-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev rno --species-verbosename norway_rat --mirna-seqs ../../mirna_reference/rno-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev sha --species-verbosename tasmanian_devil --mirna-seqs ../../mirna_reference/sha-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev sko --species-verbosename saccoglossus --mirna-seqs ../../mirna_reference/sko-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev spu --species-verbosename purple_sea_urchin --mirna-seqs ../../mirna_reference/spu-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev sto --species-verbosename cloudy_catshark --mirna-seqs ../../mirna_reference/sto-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev tca --species-verbosename red_flour_beetle --mirna-seqs ../../mirna_reference/tca-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev tgu --species-verbosename zebra_finch --mirna-seqs ../../mirna_reference/tgu-pre.fas
python3 generate-mirtrace-rnatype-database.py --out-dir ../../custom_databases/ --species-abbrev xtr --species-verbosename tropical_clawed_frog --mirna-seqs ../../mirna_reference/xtr-pre.fas

rm ../../custom_databases/*.rrna.db.gz
rm ../../custom_databases/*.trna.db.gz
rm ../../custom_databases/*.artifacts.db.gz
