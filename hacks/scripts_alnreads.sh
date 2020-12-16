#! /usr/bin/bash
# Chunyu Zhao 2020-08-25
set -e
set -x

if [ $# -ne 2 ]; then
    echo "Usage: $0 SPECIES_ID BT2_COMBO"
    exit 1
fi

species_id="$1"
bt2_combo="$2"
sim_cov=20


base_dir="/mnt/chunyu/prep_s3/${species_id}"
reads_dir="${base_dir}/reads"
genome_under_investigation=`cat ${base_dir}/genome_under_investigation`

project_dir="${base_dir}/midas/${bt2_combo}"
logs_dir="${project_dir}/logs"
midas_db_dir="${project_dir}/db"

mkdir -p ${reads_dir}
mkdir -p ${logs_dir}


################ build the bowtie2 database
genomes_toc="${midas_db_dir}/genomes.tsv"
species_ids="${midas_db_dir}/species_ids"
tail -n +2 ${genomes_toc} | awk '{print $2}' > ${species_ids}

python -m iggtools build_bowtie2_indexes \
  --midas_iggdb ${midas_db_dir} --num_cores 16 \
  --bt2_indexes_name repgenomes --species_list ${species_ids} \
  --bt2_indexes_dir ${project_dir}/db


################## genome sequence of genome_under_investigation
genome_name="${genome_under_investigation%.fna}"
genome_fp="${base_dir}/${genome_under_investigation}"
if [ ! -f ${genome_fp} ]; then
  aws s3 cp s3://microbiome-igg/2.0/cleaned_imports/${species_id}/${genome_name}/${genome_name}.fna.lz4 - | lz4 -dc > ${genome_fp}
fi

################## simulate reads
sim_read_length=125
sim_frag_length=400
sim_frag_std=40

sim_genome_len=`cat ${genome_fp} | awk '$0 ~ ">" {if (NR > 1) {print c; } c=0; printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' | awk -F'\t' '{s+=$2} END {print s}'`
sim_read_counts=`awk "BEGIN {read_length=${sim_read_length}; coverage=${sim_cov}; genome_length=${sim_genome_len}; read_counts = genome_length * coverage / read_length; print read_counts}"`


reads_name="art_${sim_cov}"
r1="${reads_dir}/${reads_name}_1.fastq"
r2="${reads_dir}/${reads_name}_2.fastq"
if [ ! -f ${r1} ]; then
  art_illumina -ss HS25 -i ${genome_fp} -l ${sim_read_length} -f ${sim_cov} -p -m ${sim_frag_length} -s ${sim_frag_std} -sp -sam -o ${reads_dir}/${reads_name}
  mv ${reads_dir}/${reads_name}1.fq ${r1}
  mv ${reads_dir}/${reads_name}2.fq ${r2}
fi


################## midas_run_snps for
sample_name="art_${bt2_combo}_${sim_cov}X"
python -m iggtools midas_run_snps --sample ${sample_name} \
  -1 ${r1} -2 ${r2} --midas_iggdb ${midas_db_dir} \
  --prebuilt_bowtie2_indexes ${midas_db_dir}/repgenomes \
  --prebuilt_bowtie2_species ${midas_db_dir}/repgenomes.species \
  --marker_depth=-1 --num_cores 16 \
  ${project_dir}/out &> ${logs_dir}/${sample_name}_snps.log
