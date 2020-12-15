import json
from collections import defaultdict
from operator import itemgetter
import multiprocessing
from math import ceil
import argparse

from iggtools.common.utils import tsprint, command, InputStream, OutputStream, select_from_tsv
from iggtools.params.schemas import snps_pileup_schema


DEFAULT_ALLELE_DEPTH = 2
DEFAULT_SITE_TYPES = ["ref", "snp"]


true_sites_schema = {
    "refid": str,
    "qryid": str,
    "pos1": str,
    "pos2": str,
    "refallele": str,
    "altallele": str,
    "sitetype": str,
}



sites_stats_schema = {
    "refid": str,
    "sitetype": str,
    "sites_total": int,
    "sites_uncalled": int,
    "sites_correct": int,
    "sites_incorrect": int,
}


def get_site_id(ref_id, ref_pos, ref_allele):
    return f"{ref_id}|{ref_pos}|{ref_allele}"


def gen_myid(ref_id, site_type):
    return f"{ref_id}|{site_type}"


def read_true_sites(truesites_fp):
    # Read-only global variables
    global true_sites_dict
    true_sites_dict = defaultdict(dict)
    sites_summary_dict = defaultdict(dict)
    with InputStream(truesites_fp) as stream:
        for row in select_from_tsv(stream, selected_columns=true_sites_schema, result_structure=dict):
            ref_id = row["refid"]
            ref_pos = str(int(float(row["pos1"])))
            site_id = get_site_id(ref_id, ref_pos, row["refallele"])
            true_sites_dict[site_id] = row
            site_type = row["sitetype"]
            my_id = gen_myid(ref_id, site_type)
            if my_id in sites_summary_dict:
                sites_summary_dict[my_id]["sites_total"] += 1
            else:
                sites_stats = {
                    "refid": ref_id,
                    "sitetype":site_type,
                    "sites_total": 1,
                    "sites_uncalled": 0,
                    "sites_correct": 0,
                    "sites_incorrect": 0,
                }
                sites_summary_dict[my_id] = sites_stats
    return (true_sites_dict, sites_summary_dict)


def call_major_allele(tuple_of_alleles, min_allele_depth):
    # Keep alleles passing the min allele depth
    alleles_above_cutoff = tuple(al for al in tuple_of_alleles if al[1] >= min_allele_depth)
    # un-callable sites
    number_alleles = len(alleles_above_cutoff)
    if number_alleles == 0:
        return (None, 0)
    alleles_above_cutoff = sorted(alleles_above_cutoff, key=itemgetter(1), reverse=True)[:2]
    major_allele = alleles_above_cutoff[0][0]
    major_allele_readcounts = alleles_above_cutoff[0][1]
    if number_alleles > 1 and major_allele_readcounts == alleles_above_cutoff[-1][1]:
        # in the event of a tie => undetermined
        return (major_allele, -major_allele_readcounts)
    # return major allele and read counts
    return (major_allele, major_allele_readcounts)


def process_pileup(pileup_fp, truesites_fp, min_allele_depth):
    true_sites_dict, sites_summary_dict = read_true_sites(truesites_fp)
    with InputStream(pileup_fp) as stream:
        for row in select_from_tsv(stream, selected_columns=snps_pileup_schema, result_structure=dict):
            ref_id = row["ref_id"]
            ref_pos = row["ref_pos"]
            site_id = get_site_id(ref_id, ref_pos, row["ref_allele"])
            if site_id in true_sites_dict:
                site_type = true_sites_dict[site_id]["sitetype"]
                my_id =gen_myid(ref_id, site_type)
                # current site is covered by at least one reads
                acgt_tuple = (('A', row["count_a"]), ('C', row["count_c"]), ('G', row["count_g"]), ('T', row["count_t"]))
                major_allele, major_allele_readcounts = call_major_allele(acgt_tuple, min_allele_depth)
                if major_allele_readcounts <= 0:
                    sites_summary_dict[my_id]["sites_uncalled"] += 1
                else:
                    alt_allele = true_sites_dict[site_id]["altallele"]
                    if major_allele == alt_allele:
                        sites_summary_dict[my_id]["sites_correct"] += 1
                    else:
                        sites_summary_dict[my_id]["sites_incorrect"] += 1
    return sites_summary_dict


def main():
    p = argparse.ArgumentParser(prog="python summarize_pileup.py", description='compute sites level summary per pileup file')
    p.add_argument(
        "--pileup_fp", required=True,
        type=str,
        help=f"Path of pileup file.")
    p.add_argument(
        "--truesites_fp", required=True,
        type=str,
        help=f"Path to reference set files.")
    p.add_argument(
        "--summary_fp", required=True,
        type=str,
        help=f"Path to output summary file.")
    p.add_argument(
        "--min_allele_depth",
        type = int, default = DEFAULT_ALLELE_DEPTH,
        help=f"Number of physical cores to use ({DEFAULT_ALLELE_DEPTH})")

    args = p.parse_args()
    sites_summary_dict = process_pileup(args.pileup_fp, args.truesites_fp, args.min_allele_depth)
    with OutputStream(args.summary_fp) as stream:
        stream.write("\t".join(sites_stats_schema.keys()) + "\n")
        for _, record in sites_summary_dict.items() :
            stream.write("\t".join(map(str, record.values())) + "\n")

main()
