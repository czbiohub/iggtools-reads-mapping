#!/bin/bash
set -e
set -x

if [ $# -ne 2 ]; then
    echo "Usage: $0 SPECIES_ID THREADS"
    exit 1
fi

species_id="$1"
threads="$2"


datadir=/mnt/cz/20201110_species_w_isolates/${species_id}

snakemake --configfile /mnt/cz/iggtools-reads-mapping/config.yml --config species_id=${species_id} -p --cores ${threads} _prepare

snakemake --configfile /mnt/cz/iggtools-reads-mapping/config.yml --config species_id=${species_id} -p --cores ${threads} _reads

snakemake --configfile /mnt/cz/iggtools-reads-mapping/config.yml --config species_id=${species_id} -p --cores ${threads} _snps

snakemake --configfile /mnt/cz/iggtools-reads-mapping/config.yml --config species_id=${species_id} -p --cores ${threads} _major --rerun-incomplete


rm -r ${datadir}/4_midas/*/out/*/temp
rm -r ${datadir}/4_midas/*/db/*bt2*


aws s3 cp --recursive ${datadir} s3://czhao-bucket/2020-bt2-reads/20201110_species_w_isolates/${species_id}

rm -r ${datadir}/[123]*
