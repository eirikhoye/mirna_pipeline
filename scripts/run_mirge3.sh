#/bin/bash

miRge3.0 -s $1 -lib miRge3_Lib -on human -db mirgenedb -o output_dir -tcf -cpu 4 -a illumina
