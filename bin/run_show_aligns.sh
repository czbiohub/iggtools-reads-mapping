#! /usr/bin/bash
# Chunyu Zhao 2020-10-21
set -e

if [ $# -ne 3 ]; then
    echo "Usage: $0 DELTA_FILE REFERENCE_GENOME QUERY_GENOME "
    exit 1
fi

delta="$1"
ref="$2"
qry="$3"


workdir=`dirname "${delta}"`
mkdir -p ${workdir}/${ref}

show-aligns -r -m 0 ${delta} ${ref} ${qry} | sed '/^$/d' > ${workdir}/${ref}/${qry}.aln
