####
# Evaluate genome selections effect for alignment based SNP calling
# author: Chunyu Zhao chunyu.zhao@czbiohub.org
# time: 2020-09-16
####

import os
import sys
import configparser
from collections import defaultdict
from Bio import SeqIO
from math import ceil


PROJ_DIR = config["base_dir"] + "/" + str(config["species_id"])
LIB_DIR = config["lib_dir"]

with open(PROJ_DIR + "/" + config["bt2_combo_fp"]) as stream:
    BT2_COMBO = [line.rstrip() for line in stream]

SIM_COV_LIST = range(1, int(config["sim_cov"]) + 1)


RC_DICT = defaultdict(int)
if os.path.exists(PROJ_DIR + "/sim_rc.tsv"):
    with open(PROJ_DIR + "/sim_rc.tsv") as stream:
        for line in stream:
            line = line.rstrip().split('\t')
            RC_DICT[line[0]] = int(line[1])


include: "reads.rules"
include: "snps.rules"
include: "nucmer.rules"


rule _tsv_v2:
    input:
        PROJ_DIR + "/6_rdata/" + str(config["species_id"]) + "_sites_summary.tsv",
        PROJ_DIR + "/6_rdata/" + str(config["species_id"]) + "_reads_summary.tsv"


rule reference_sets:
    input:
        PROJ_DIR + "/5_nucmer/" + str(config["species_id"]) + "_aligned_sites.tsv"
    output:
        PROJ_DIR + "/6_rdata/" + str(config["species_id"]) + "_true_sites.tsv"
    params:
        R = LIB_DIR + "/aln_to_truesites.R",
    shell:
        "Rscript {params.R} {input} {output}"


rule sites_summary:
    input:
        expand(PROJ_DIR + "/4_midas/{bt2_combo}/out/art_{bt2_combo}_{{sim_cov}}X/snps/" + str(config["species_id"]) + ".snps.tsv.lz4", bt2_combo=BT2_COMBO),
        true_sites = PROJ_DIR + "/6_rdata/" + str(config["species_id"]) + "_true_sites.tsv"
    output:
        PROJ_DIR + "/6_rdata/" + str(config["species_id"]) + "/sites_summary_{sim_cov}X.tsv",
    params:
        R = LIB_DIR + "/pileup_to_summary.R"
    threads:
        1
    shell:
        """
        Rscript {params.R} {PROJ_DIR} {config[species_id]} {input.true_sites} {wildcards.sim_cov}
        """


rule per_species_summary:
    input:
        expand(PROJ_DIR + "/6_rdata/" + str(config["species_id"]) + "/sites_summary_{sim_cov}X.tsv", sim_cov=SIM_COV_LIST)
    output:
        PROJ_DIR + "/6_rdata/" + str(config["species_id"]) + "_sites_summary.tsv",
        PROJ_DIR + "/6_rdata/" + str(config["species_id"]) + "_reads_summary.tsv"
    params:
        R = LIB_DIR + "/xs_cov_summary.R",
    threads:
        1
    shell:
        """
        Rscript {params.R} {PROJ_DIR} {config[species_id]}
        """
