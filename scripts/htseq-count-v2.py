import HTSeq
import itertools
import pandas as pd
import os

gtf_file = HTSeq.GFF_Reader( "data/mirna_reference/hsa_mature.gff3" )
mature_miRNA = HTSeq.GenomicArrayOfSets( "auto", stranded=False )

for feature in gtf_file:
    if feature.type == "miRNA":
        mature_miRNA[ feature.iv ] += feature.attr["ID"]#

import collections
counts = collections.Counter( )

for bam_file in os.listdir('data/bamfiles/'):
    if not bam_file.endswith('.bam'):
        continue
    print(bam_file)
    almnt_file = HTSeq.BAM_Reader( "/home/jcdenton/projects/mirna_pipeline/data/bamfiles/" + bam_file )
    for almnt in almnt_file:
        if not almnt.aligned:
            counts[ "_unmapped" ] += 1
            continue
        gene_ids = set()
        for iv, val in mature_miRNA[ almnt.iv ].steps():
            gene_ids |= val
        if len(gene_ids) == 1:
            gene_id = list(gene_ids)[0]
            counts[ gene_id ] += 1
        elif len(gene_ids) == 0:
            counts[ "_no_feature" ] += 1
        else:
            counts[ "_ambiguous" ] += 1
    print('writing_' + bam_file + 'to data/mirna_counts/')
    with open('/home/jcdenton/projects/mirna_pipeline/data/mirna_counts/' + bam_file + '.csv', 'w') as outfile:
        for gene_id in counts:
            outfile.write(str(gene_id))
            outfile.write(',')
            outfile.write(str(counts[ gene_id ]))
            outfile.write('\n')
