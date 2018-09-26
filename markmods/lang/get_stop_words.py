###imports
import pickle
import unicodedata
from hmm import HMM
from math import *
import re
############
###functions
def reader (f, cl, s_symb, s_segr, matr, par_l, r_f):
    log_emis = log_matrix(matr)
    log_params = get_params_log(*par_l)
    hmm = HMM(log_emis, *log_params)
    clust_num_d = {}
    splitter = get_splitter(s_segr)
    read = open(f)
    for lines in read:
        line = lines.strip().decode('cp1251').lower()
        if len(line) == 0:
            line = None
            continue
        sent_list = re.split(splitter, line)
        for sents in sent_list:
            sent = sents.strip()
            if len(sent) == 0:
                sent = None
                continue 
            clust_num_d = sentance_an(sent, clust_num_d, hmm, s_symb)
            sent = None
    return clust_num_d
###
def get_splitter (spl_list):
    spl = []
    for s in spl_list:      
        spl.append(s)
    return "|".join(spl)
###
def sentance_an (s, cl_num_d, hmm, s_s):
    s_list = s.split()
    if "|" in s_list[0]:
        a = s_list.pop(0)
        a = None
    words_l = []
    splitter = get_splitter(s_s)
    for w in s_list:
        w_list = re.split(splitter, w)
        for words in w_list:
            word = words.strip()
            if len(word) <= 3:#this words shoud be analyzed manually
                word = None
                continue  
            words_l.append(word)
            word = None
    for word in words_l:
        if word in cl_num_d.keys():
            cl_num_d[word][word] = cl_num_d[word][word] + 1
        else:
            cl_num_d[word] = {word:1}
        for w in words_l:
            if w == word:
                continue
            if w in cl_num_d[word].keys():
                cl_num_d[word][w] = cl_num_d[word][w] + 1
            else:
                cl_num_d[word][w] = 1
    return cl_num_d                
###
def log_matrix (matr):
    log_matr = {}
    for k in matr.keys():
        for ki in matr[k].keys():
            if k in log_matr.keys():
                log_matr[k][ki] = log(matr[k][ki],2)
            else:
                log_matr[k] = {ki:log(matr[k][ki],2)}
    return log_matr
###
def get_params_log (*par_l):
    log_par = []
    for p in par_l:
        log_par.append(log(p,2))
    return log_par
###
def define_s_words (s_dict, s_ws_n, sig_cl, r_f):
    surr_size = {}
    cl_k = s_dict.keys()
    for k in cl_k:
        n = 0
        cl_freq = s_dict[k][k]
        if cl_freq < sig_cl:
            continue
        surr_k = s_dict[k].keys()
        for ki in surr_k:
            if ki == k:
                continue
            if s_dict[k][ki] >= sig_cl:
                n = n + 1
        if n == 0:
            continue
        surr_size[k] = float(n)/(1.0+log10(cl_freq))
    hash_dict = hashing(surr_size)
    keys_l = hash_dict.keys()
    keys_l.sort(reverse=True)
    i = 0
    while i < s_ws_n:
        if i >= len(keys_l):
            break
        j = 0
        for cl in hash_dict[keys_l[i]]:
            if i+j >= s_ws_n:
                break
            writer(r_f, cl.encode("cp1251"))
            j = j + 1
        i = i + j
###
def hashing (in_dict):
    hash_dict = {}    
    ss_k = in_dict.keys()
    for k in ss_k:
        if in_dict[k] in hash_dict.keys():
            hash_dict[in_dict[k]].append(k) 
        else:
            hash_dict[in_dict[k]] = [k]
    return hash_dict
###
def writer (f, s):
    w = open (f, "a")
    w.write(s + "\n")
    w.close()
############
###MAIN
def main (params):
###parameters
    clust_f = "clusters.p"
    out_f = "stop_words.txt"
    stop_ws_number = params["stop_words_number"]
    signif_clust = params["signif_clust"]
    in_text = params["input_file"]
    matr_f = "draft_matrix.p"
    stop_symb = params["spec_symbs"]
    sent_segr = params["phr_sep"]
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
    len_bonus_co = params["aln_glocality"]
    len_bonus = ((2.0**len_bonus_co)/((match_s_prob**3.0)*(match_to_match**2.0)*(av_match**3.0)*gap_open*no_aln_fine))**0.25
    par_list = [match_s_prob, match_to_match, gap_open, gap_ext, to_exit_prob, start_prob, ends_bonus, gap_fine, no_aln_fine, exit_prob, len_bonus]
############
    wr = open(out_f, 'w')
    wr.close()
    matr_read = open(matr_f, 'rb')
    matrix = pickle.load(matr_read)
    cl_read = open(clust_f, 'rb')
    clusters = pickle.load(cl_read)
    surr_dict = reader(in_text, clusters, stop_symb, sent_segr, matrix, par_list, out_f)
    define_s_words(surr_dict, stop_ws_number, signif_clust, out_f)