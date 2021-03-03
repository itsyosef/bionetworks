#!/bin/bash

#SBATCH -A nssac_students
#SBATCH -t 04:00:00
#SBATCH --cpus-per-task=10
#SBATCH -p bii

set -e

module load singularity anaconda parallel
module load gcc bwa/0.7.17 samtools
source activate operon

cd /scratch/jho5ze/bionets/operons/SRA

info_file=/scratch/jho5ze/bionets/operons/SRA/genera/${genus}/${genus}_query_info.csv

script_path="/scratch/jho5ze/bionets/operons/SRA"
patric_mash="/project/biocomplexity/fungcat/anwarren_locker/mash/patric_all.msh"

readarray -t projects < <(cat $info_file | cut -f 2 -d "," | sort -u | grep -v "study_accession")
for (( j=0; j<${#projects[@]}; j++)) do
    set -e

    project=${projects[$j]}
    
    readarray -t orgs < <(grep $project $info_file | cut -f 5 -d "," | sort -u | grep -v "sample_taxon_id")
    
    path=/scratch/jho5ze/bionets/operons/SRA/salmon_pipeline_devel/$genus/$project
    mkdir -p $path

    for (( i=0; i<${#orgs[@]}; i++)) do 
        cd /scratch/jho5ze/bionets/operons/SRA
        org=${orgs[$i]}
        #To make sure that the taxon id is actually present in the list of refseq genomes (sometimes the pysradb results include mixed populations with taxids not represented in the refseq dir (or in any single reference), see query "Escherichia RNA-Seq" first 300 results)
#         echo $org
#         echo $project
#         echo "RELEVANT REFSEQ LINE FROM patric_genome_metadata.txt"
    #     grep -P "\t$org" /scratch/jho5ze/bionets/operons/assembly_summary_refseq.txt
#         grep -P "^$org.[0-9]*\t" /sfs/lustre/bahamut/scratch/jho5ze/bionets/operons/patric_genome_metadata.txt
#         if ! grep --quiet -P "^$org.[0-9]*\t" /sfs/lustre/bahamut/scratch/jho5ze/bionets/operons/patric_genome_metadata.txt ; then
#             echo "Failed search"
#             continue
#         fi
#         echo "Getting line:"
#         grep -P "^$org.[0-9]*\t" /sfs/lustre/bahamut/scratch/jho5ze/bionets/operons/patric_genome_metadata.txt | head -n 1
#     #     read genome_accession ftp_path <<<$(get_reftp "$org")
#         genome_accession=$(get_patric_id "$org")
#         echo "$path/$genome_accession"
#         echo $genome_accession
#         if [ ! -d "$path/$genome_accession" ]; then
#             $annogesic create --project_path $path/$genome_accession
#     #         $annogesic get_input_files --project_path $path/$genome_accession --ftp_path $ftp_path/ --ref_gff --ref_gbk --ref_fasta
#             wget -qN -P $path/$genome_accession/input/references/fasta_files ftp://ftp.patricbrc.org/genomes/$genome_accession/$genome_accession.fna
#             wget -qN -P $path/$genome_accession/input/references/annotations ftp://ftp.patricbrc.org/genomes/$genome_accession/$genome_accession.gff || (wget -qN -O $path/$genome_accession/input/references/annotations/$genome_accession.gff ftp://ftp.patricbrc.org/genomes/$genome_accession/$genome_accession.PATRIC.gff && sed -i 's/^accn|//g' $path/$genome_accession/input/references/annotations/$genome_accession.gff)
#             awk '{Q=$0;gsub(/CDS/, "gene", Q);gsub(/ID=/, "ID=gene_for-", Q); if ($0 ~ "CDS") print $0 ORS Q; else print $0}' $path/$genome_accession/input/references/annotations/$genome_accession.gff > $path/$genome_accession/input/references/annotations/$genome_accession.gff.temp; mv $path/$genome_accession/input/references/annotations/$genome_accession.gff.temp $path/$genome_accession/input/references/annotations/$genome_accession.gff
            
#             cd $path/$genome_accession/input/references/fasta_files
#     #         ls *.bwt
#             if [ "$(ls -1 *.bwt 2>/dev/null | wc -l)" -eq "0" ]; then
#                 ref_seq=$(ls *.f*a | head -n 1)
#                 bwa index $ref_seq
#             fi

#             echo "Moving!"
#             pwd
#             cd $path/$genome_accession/input/reads
#             pwd
#             echo "grep ${genome_accession%.*} $info_file | cut -f 15 -d "," | parallel prefetch "
            
#             grep ${genome_accession%.*} $info_file | cut -f 15 -d "," | parallel prefetch 
#     #         for srr_file in $(ls -d $(pwd)/SRR*/*); do
#     #             parallel-fastq-dump -s $srr_file --threads 10 --outdir . --gzip 
#             parallel-fastq-dump -s $(ls -d $(pwd)/*RR*/*) --threads 10 --outdir . --gzip  --split-e
#             #ToDo: Handle paired end reads
# #             for xRR in $(ls *RR*.gz | cut -f 1 -d "." | cut -f 1 -d "_" | sort -u); do
# #                 echo $xRR
# #                 if [ "$(ls ${xRR}*.gz | wc -l)" -eq "3" ]; then
# #                     rm $xRR.fastq.gz
                    
# #                 fi
# #             done
#             for file in $(ls *.gz | grep -v "trimmed"); do
#                 fastp -i $file -o ${file%.fastq*}_trimmed.fastq.gz -q 25

#             done
#         fi

        break
    done
    break
    
done