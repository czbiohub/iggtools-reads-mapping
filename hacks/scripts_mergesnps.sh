#! /usr/bin/bash
# Chunyu Zhao 2020-08-25
set -e
set -x

if [ $# -ne 2 ]; then
    echo "Usage: $0 SPECIES_ID SIM_COV"
    exit 1
fi

species_id="$1"
sim_cov="$2"


base_dir="/mnt/chunyu/prep_s3/${species_id}"
merge_dir="${base_dir}/merge/${sim_cov}X"
mkdir -p ${merge_dir}

midas_basedir="${base_dir}/midas"
midas_db_dir="${base_dir}/midas/ani.100/db"


logs_dir="${project_dir}/logs"


sample_list="${merge_dir}/samples_list.tsv"
echo -e  "sample_name\tmidas_outdir" > ${sample_list}
paste <(ls ${midas_basedir} | awk -v sc=$sim_cov 'BEGIN {OFS=""} {print "art_", $1, "_", sc, "X"}') \
  <(ls ${midas_basedir} | awk -v bd=$midas_basedir 'BEGIN {OFS=""} {print bd, "/", $1, "/out"}') >> ${sample_list}



################## midas_run_snps for given sim_cov
python -m iggtools midas_merge_snps \
  --samples_list ${sample_list} --species_list ${species_id} \
  --midas_iggdb ${midas_db_dir} \
  --site_depth 1 --site_ratio 500 \
  --site_prev 0 --site_type common \
  --snp_pooled_method abundance --debug \
  --snp_maf 0.001 --snp_type=any ${merge_dir} # &> ${merge_dir}/merge_snps.log
