set -e

if [ $# -ne 2 ]; then
    echo "Usage: $0 SPECIES_ID THREADS"
    exit 1
fi

species_id="$1"
threads="$2"


datadir="/mnt/chunyu/20201030_species_w_isolates/${species_id}"

snakemake --configfile /mnt/chunyu/snps_eval/config.yml --config species_id=${species_id} -p --cores ${threads} _prepare

snakemake --configfile /mnt/chunyu/snps_eval/config.yml --config species_id=${species_id} -p --cores ${threads} _reads

snakemake --configfile /mnt/chunyu/snps_eval/config.yml --config species_id=${species_id} -p --cores ${threads} _snps

snakemake --configfile /mnt/chunyu/snps_eval/config.yml --config species_id=${species_id} -p --cores ${threads} _tsv --rerun-incomplete


#rm -r "${datadir}/4_midas/*/out/*/temp"
#rm -r "${datadir}/4_midas/*/db/*bt2*"


#aws s3 cp --recursive ${datadir} s3://czhao-bucket/2020-bt2-reads/20201030_species_w_isolates/${species_id}
