#!/bin/bash

set -e

source /sfs/lustre/bahamut/scratch/jho5ze/bionets/BGCs/scripts/venv/bin/activate

# python3 anti_gids_to_tsv.py $@

outprefix=$1
shift

python /sfs/lustre/bahamut/scratch/jho5ze/bionets/patric_genera/fork/pangenome_graphs/fam_to_graph.py --layout --output ${outprefix}_with_pfams.gexf --alpha pgfam_id ${outprefix}_with_pfams.tsv

python3 anti_annotate_graph.py ${outprefix}_with_pfams.gexf
