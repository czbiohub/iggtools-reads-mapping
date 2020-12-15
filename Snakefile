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

SIM_COV_LIST = range(1, int(config["sim_cov"]))


RC_DICT = defaultdict(int)
if os.path.exists(PROJ_DIR + "/sim_rc.tsv"):
    with open(PROJ_DIR + "/sim_rc.tsv") as stream:
        for line in stream:
            line = line.rstrip().split('\t')
            RC_DICT[line[0]] = int(line[1])


include: "rules/reads.rules"
include: "rules/snps.rules"
include: "rules/nucmer.rules"
include: "rules/summary.rules"


rule _stats:
    input:
        PROJ_DIR + "/7_rdata/" + str(config["species_id"]) + "_sites_summary.tsv",
        PROJ_DIR + "/7_rdata/" + str(config["species_id"]) + "_reads_summary.tsv",
        #PROJ_DIR + "/7_rdata/" + str(config["species_id"]) + "_sites_summary_nofilter",
        #PROJ_DIR + "/7_rdata/" + str(config["species_id"]) + "_reads_summary_nofilter",


rule _dbs:
    input:
        expand(PROJ_DIR + "/4_midas/{bt2_combo}/db/repgenomes.species", bt2_combo = BT2_COMBO)


rule _snps:
    input:
        [expand(PROJ_DIR + "/4_midas/{bt2_combo}/out/art_{bt2_combo}_{sim_cov}X/snps/" + str(config["species_id"]) + ".snps.tsv.lz4",
            bt2_combo = BT2_COMBO, sim_cov = SIM_COV_LIST),
        expand(PROJ_DIR + "/4_midas_nofilter/{bt2_combo}/out/art_{bt2_combo}_{sim_cov}X/snps/" + str(config["species_id"]) + ".snps.tsv.lz4",
            bt2_combo = BT2_COMBO, sim_cov = SIM_COV_LIST)]


rule _bamaln:
    input:
        expand(PROJ_DIR + "/7_rdata/" + str(config["species_id"]) + "/{bt2_combo}/bamalns_{sim_cov}X.tsv", bt2_combo = BT2_COMBO, sim_cov = range(15, 21))


rule _true_sites:
    input:
        PROJ_DIR + "/7_rdata/" + str(config["species_id"]) + "_true_sites.tsv"


rule _sites:
    input:
        PROJ_DIR + "/5_nucmer/" + str(config["species_id"]) + "_aligned_sites.tsv"


rule _aligns:
    input:
        PROJ_DIR + "/5_nucmer/DONE-show-aligns"


rule _prepare:
    input:
        PROJ_DIR + "/sim_rc.tsv"


rule _genomes:
    input:
        [PROJ_DIR + "/genome_under_investigation.fna", PROJ_DIR + "/uhgg_rep.fna"]


rule _reads:
    input:
        expand(PROJ_DIR + "/3_reads/cov_{sim_cov}_1.fastq", sim_cov=SIM_COV_LIST)
