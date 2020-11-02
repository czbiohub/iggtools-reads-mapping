#! /usr/bin/bash
# Chunyu Zhao

set -e
set -x

if [ $# -ne 2 ]; then
    echo "Usage: $0 GENOME SAMPLE_NAME"
    exit 1
fi

genome="$1"
sample="$2"

sample_name="${sample%_1.fastq}"
genome_name="${genome%.fna}"


project_dir="/mnt/chunyu_6TB/snps/20200819_eval_simreads/${genome_name}"
genome_fp="${project_dir}/data/${genome}"
models_fp="/mnt/chunyu_6TB/snps/prototype_1_reads_simulator/models_ecoli_trained"


mkdir -p $project_dir
outdir="${project_dir}/reads"
mkdir -p ${outdir}
mkdir -p ${project_dir}/log



###################### Compare multiple reads simulator
# InSilicoSeq
sim_genlen=`cat ${genome_fp} | awk '$0 ~ ">" {if (NR > 1) {print c; } c=0; printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' | awk -F'\t' '{s+=$2} END {print s}'`
sim_read_counts=`awk "BEGIN {read_length=125; coverage=50; genome_length=${sim_genlen}; read_counts = genome_length * coverage / read_length; print read_counts}"`
# simulate reads without sequencing error
iss generate --draft ${genome_fp} --mode perfect --output ${outdir}/iss_perfect --cpus 16 --abundance uniform --n_reads ${sim_read_counts} # --compress
# simulate HiSeq reads
iss generate --draft ${genome_fp} --model HiSeq --output ${outdir}/iss_hiseq --cpus 16 --abundance uniform --n_reads ${sim_read_counts} # --compress
# simulate based on raw ecoli genome
iss generate --draft ${genome_fp} --model ${models_fp}/iss/iss_ecoli_trained.npz --output ${outdir}/iss_trained --cpus 16 --abundance uniform --n_reads ${sim_read_counts}


# neat
# use the prebuilt models from ecoli
python genReads.py -r ${genome_fp} -p 1 -R 126 -c 50 -e ${models_fp}/neat/model_seq_error.p --pe-model ${models_fp}/neat/model_fraglen.p -o ${outdir}/neat_trained


# art
art_bin_MountRainier/art_illumina -ss HS25 -i ${genome_fp} -l 125 -f 50 -p -m 400 -s 40 -sp -sam -o ${outdir}/art_hiseq

# reseq


###################### Moved to scripts_alnreads.sh
# sequencing reads
python -m iggtools midas_run_snps --sample ${sample_name} \
  -1 ${outdir}/${sample_name}_1.fastq \
  -2 ${outdir}/${sample_name}_2.fastq \
  --midas_iggdb ${project_dir}/midas_db \
  --prebuilt_bowtie2_indexes ${project_dir}/data/repgenomes \
  --prebuilt_bowtie2_species ${project_dir}/data/repgenomes.species \
  --marker_depth=-1 --num_cores 16 \
  ${project_dir}/midas_out &> ${project_dir}/log/${sample_name}.log



# prebuild the bowtie2 database
python -m iggtools build_bowtie2_indexes \
  --midas_iggdb ${project_dir}/midas_db \
  --bt2_indexes_name repgenomes --species_csv ${species_id} \
  --bt2_indexes_dir ${project_dir}/data
