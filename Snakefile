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

#include: "snps.rules"
#include: "nucmer.rules"


rule _tsv:
    input:
        PROJ_DIR + "/7_rdata/" + str(config["species_id"]) + "_detection_power.tsv",
        PROJ_DIR + "/7_rdata/" + str(config["species_id"]) + "_sites_presence.tsv",
        PROJ_DIR + "/7_rdata/" + str(config["species_id"]) + "_sites_abundance.tsv",
        PROJ_DIR + "/7_rdata/" + str(config["species_id"]) + "_sites_abuncorr.tsv",


rule _ana:
    input:
        expand(PROJ_DIR + "/7_rdata/" + str(config["species_id"]) + "/sites_allelefreq_{sim_cov}X.tsv", sim_cov=SIM_COV_LIST)


rule _major:
    input:
        expand(PROJ_DIR + "/6_rtemp/" + "maj.wide_{sim_cov}X.tsv", sim_cov=SIM_COV_LIST)


rule major_pass_one:
    input:
        expand(PROJ_DIR + "/4_midas/{bt2_combo}/out/art_{bt2_combo}_{{sim_cov}}X/snps/" + str(config["species_id"]) + ".snps.tsv.lz4", bt2_combo=BT2_COMBO),
        aligned_sites = PROJ_DIR + "/5_nucmer/" + str(config["species_id"]) + "_aligned_sites.tsv"
    output:
        PROJ_DIR + "/6_rtemp/" + "maj.wide_{sim_cov}X.tsv",
        PROJ_DIR + "/7_rdata/" + str(config["species_id"]) + "/maj.stats_{sim_cov}X.tsv",
    params:
        maj_pass_one_R = LIB_DIR + "/maj_pass_one.R",
    threads:
        2
    shell:
        """
        Rscript {params.maj_pass_one_R} {PROJ_DIR} {config[species_id]} {input.aligned_sites} {wildcards.sim_cov}
        """


rule analysis_pass_one:
    input:
        PROJ_DIR + "/6_rtemp/" + "maj.wide_{sim_cov}X.tsv"
    output:
        PROJ_DIR + "/7_rdata/" + str(config["species_id"]) + "/sites_abundance_{sim_cov}X.tsv",
        PROJ_DIR + "/7_rdata/" + str(config["species_id"]) + "/sites_abuncorr_{sim_cov}X.tsv",
        PROJ_DIR + "/7_rdata/" + str(config["species_id"]) + "/sites_presence_{sim_cov}X.tsv",
        PROJ_DIR + "/7_rdata/" + str(config["species_id"]) + "/sites_allelefreq_{sim_cov}X.tsv"
    params:
        maj_pass_one_R = LIB_DIR + "/maj_ana_one.R",
    threads:
        2
    shell:
        """
        Rscript {params.maj_pass_one_R} {PROJ_DIR} {config[species_id]} {wildcards.sim_cov}
        """



rule major_pass_two:
    input:
        expand(PROJ_DIR + "/6_rtemp/" + "maj.wide_{sim_cov}X.tsv", sim_cov=SIM_COV_LIST)
    output:
        PROJ_DIR + "/7_rdata/" + str(config["species_id"]) + "_detection_power.tsv",
        PROJ_DIR + "/7_rdata/" + str(config["species_id"]) + "_sites_presence.tsv",
        PROJ_DIR + "/7_rdata/" + str(config["species_id"]) + "_sites_abundance.tsv",
        PROJ_DIR + "/7_rdata/" + str(config["species_id"]) + "_sites_abuncorr.tsv",
    params:
        gen_major_R = LIB_DIR + "/maj_pass_two.R",
    threads:
        4
    shell:
        """
        Rscript {params.gen_major_R} {PROJ_DIR} {config[species_id]}
        """
