###imports
from get_aln_blocks import main_func as get_aln_blocks
import subprocess
import multiprocessing as mp
import pickle
from math import *
import unicodedata
from hmm import HMM
import redis
import join_clust
############
###MAIN
def main (params):
###parameters
    num_threads = params["num_threads"]
    in_f = params["input_file"]
    matrix_f = "draft_matrix.p"
    db_name = params["db_number"]
    out_f = "subs_matrix.p"
    stops_list = params["spec_symbs"]
    aln_min_prob = params["aln_min_prob"]
    mismatch = params["mismatch"]
    match_s_prob = params["match_s_prob"]
    match_to_match = params["match_to_match"]
    av_match = params["av_match"]
    gap_open = params["gap_open"]
    gap_ext = params["gap_ext"]
    to_exit_prob = params["to_exit_prob"]
    start_prob = params["start_prob"]
    ends_bonus = params["ends_bonus"]
    gap_fine = params["gap_fine"]
    no_aln_fine = params["no_aln_fine"]
    exit_prob = params["exit_prob"]
    delta = params["delta"]
    n_iter = params["n_iter"]
    len_bonus_co = params["aln_glocality"]#The main parameter
    len_bonus = ((2.0**len_bonus_co)/((match_s_prob**3.0)*(match_to_match**2.0)*(av_match**3.0)*gap_open*no_aln_fine))**0.25
    par_list = [num_threads, in_f, matrix_f, db_name, out_f, stops_list, aln_min_prob, mismatch, match_s_prob, match_to_match, av_match, gap_open, gap_ext, to_exit_prob, start_prob, ends_bonus, gap_fine, no_aln_fine, exit_prob, len_bonus]
############
    get_aln_blocks(*par_list)
    join_clust.main(params)