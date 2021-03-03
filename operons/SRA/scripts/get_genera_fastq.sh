#!/bin/bash

#SBATCH -A nssac_students
#SBATCH -t 04:00:00
#SBATCH --cpus-per-task=10
#SBATCH -p bii

module load singularity anaconda parallel
module load gcc bwa/0.7.17 samtools
source activate operon

set -e

get_reftp () {
    grep -P "\t$1" /scratch/jho5ze/bionets/operons/assembly_summary_refseq.txt | cut -f 6,20 | head -n 1
    exit $?
}

#Grabs the patric id associated with the taxon id from the SRA result. Takes the patric genome id with the smallest number after the period for consistency...
get_patric_id () {
    grep -P "^$1.[0-9]*\t" /sfs/lustre/bahamut/scratch/jho5ze/bionets/operons/patric_genome_metadata.txt | cut -f 1 | sort -u | head -n 1
}

cd /scratch/jho5ze/bionets/operons/SRA


genera=($(ls operon_prediction/)) #Just work with what we have
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    SLURM_ARRAY_TASK_ID=0
fi
genus=${genera[$SLURM_ARRAY_TASK_ID]}

info_file=/scratch/jho5ze/bionets/operons/SRA/genera/${genus}/${genus}_query_info.csv
annogesic="singularity exec /scratch/jho5ze/bionets/operons/SRA/annogesic_latest.sif annogesic"
# reademption="singularity exec /scratch/jho5ze/bionets/operons/SRA/reademption_latest.sif reademption"
script_path="/scratch/jho5ze/bionets/operons/SRA"

readarray -t projects < <(cat $info_file | cut -f 2 -d "," | sort -u | grep -v "sample_taxon_id")

for (( j=0; j<${#projects[@]}; j++)) do
    
    readarray -t orgs < <(cat $info_file | cut -f 5 -d "," | sort -u | grep -v "sample_taxon_id")
    
    project=${projects[$j]}
    path=/scratch/jho5ze/bionets/operons/SRA/pipeline_devel/$genus/$project
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
        if ! grep --quiet -P "^$org.[0-9]*\t" /sfs/lustre/bahamut/scratch/jho5ze/bionets/operons/patric_genome_metadata.txt ; then
#             echo "Failed search"
            continue
        fi
    #     read genome_accession ftp_path <<<$(get_reftp "$org")
        genome_accession=$(get_patric_id "$org")
        if [ ! -d "$path/$genome_accession" ]; then
            $annogesic create --project_path $path/$genome_accession
    #         $annogesic get_input_files --project_path $path/$genome_accession --ftp_path $ftp_path/ --ref_gff --ref_gbk --ref_fasta
            wget -qN -P $path/$genome_accession/input/references/fasta_files ftp://ftp.patricbrc.org/genomes/$genome_accession/$genome_accession.fna
            wget -qN -P $path/$genome_accession/input/references/annotations ftp://ftp.patricbrc.org/genomes/$genome_accession/$genome_accession.gff || wget -qN -P $path/$genome_accession/input/references/annotations ftp://ftp.patricbrc.org/genomes/$genome_accession/$genome_accession.PATRIC.gff

            echo "Moving!"
            pwd
            cd $path/$genome_accession/input/reads
            pwd
            grep ${genome_accession%.*} $info_file | cut -f 15 -d "," | parallel prefetch 
    #         for srr_file in $(ls -d $(pwd)/SRR*/*); do
    #             parallel-fastq-dump -s $srr_file --threads 10 --outdir . --gzip 
            parallel-fastq-dump -s $(ls -d $(pwd)/SRR*/*) --threads 10 --outdir . --gzip --split-e
            for SRR in $(ls SRR* | cut -f 1 -d "." | cut -f 1 -d "_" | sort -u); do
                echo $SRR
                if [ "$(ls ${SRR}* | wc -l)" -eq "3" ]; then
                    rm $SRR.fastq.gz
                    mv ${SRR}_1.fastq.gz ${SRR}_forward.fastq.gz
                    mv ${SRR}_2.fastq.gz ${SRR}_reverse.fastq.gz
                fi
            done
            echo "Here"
#             cd $path/$genome_accession/input/reads && rm -rf $(ls -d */)
            cd $path/$genome_accession/output
            for missing_dir in {alignments,processed_reads,unaligned_reads,index,reports_and_stats/stats_data_json}; do
                mkdir -p align/$missing_dir
            done

            for missing_dir in {coverage-raw,coverage-tnoar_mil_normalized,coverage-tnoar_min_normalized}; do
                mkdir -p coverage/$missing_dir
            done

            ln -s $path/$genome_accession/input/references/fasta_files/ $path/$genome_accession/input/reference_sequences
            ln -s $path/$genome_accession/input/references/annotations/ $path/$genome_accession/input/annotation_files
    #         $reademption align -f $path/$genome_accession -g -p 10
    #         done
        fi
        #To run after downloading the SRA info
        cd $path/$genome_accession/input/references/fasta_files
        ls *.bwt
        if [ ! ls *.bwt ]; then
            ref_seq=$(ls *.f*a | head -n 1)
            bwa index $ref_seq
        fi

        cd $path/$genome_accession/input/reads
        for file in $(ls | grep -v "trimmed"); do
            fastp -i $file -o ${file%.fastq*}_trimmed.fastq.gz -q 25
            rm $file
    #         bbduk in=$file out=${file%.fastq*}_trimmed.fastq.gz ref=adapters minbasequality=25 ziplevel=2
            bwa mem -t 10 ../../reference_sequences/$ref_seq ${file%.fastq*}_trimmed.fastq.gz | samtools view -@ 10 -bh -q 25 -F 4 - | samtools sort > ${file%.fastq*}.sorted.bam
            python $script_path/bam_to_wiggle.py ${file%.fastq*}.sorted.bam -o ${file%.fastq*}.wig
#             rm $file
#             rm ${file%.fastq*}.sorted.bam
                mv ${file%.fastq*}.wig ../wigs/fragment/
        done
        cd $path
        
#         $reademption align -f $genome_accession -g -p 10
        break
    done
    break
    
done