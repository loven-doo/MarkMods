###imports
import subprocess
import multiprocessing as mp
import pickle
from math import *
import unicodedata
from hmm import HMM
import redis
import re
############
###functions
def reader (f, stops_l):
    words_l = []
    read = open(f)
    for lines in read:
        line = lines.strip().decode('cp1251').lower()
        if len(line) == 0:
            line = None
            continue
        line_list = line.split()
        if "|" in line_list[0]:
            a = line_list.pop(0)
            a = None
        line_list = rem_p_marks(line_list, stops_l)
        for w in line_list:
            words_l.append(w)
        line = None
    return words_l
###
def rem_p_marks (l_list, s_l):
    words_l = []
    splitter = get_splitter(s_l)
    for elm in l_list:
        elm_list = re.split(splitter,elm)
        for word in elm_list:
            if len(word.strip()) <= 3:
                continue
            words_l.append(word.strip())
    return words_l
###
def get_splitter (spl_list):
    spl = []
    for s in spl_list:      
        spl.append(s)
    return "|".join(spl)
###
def get_aln_blocks (pars_l):
    words_l = pars_l[0]
    text_l = pars_l[1]
    hmm = pars_l[2]
    r = redis.StrictRedis(host='localhost', port=6379, db=pars_l[3])
    min_prob =  pars_l[4]
    print pars_l[5]
    for word in words_l:
        aligned_words = r.keys('*')
        if word.encode('utf-8') in aligned_words:
            continue
        r.set(word, word.encode('cp1251'))
        for w in text_l:
            found_words = r.get(word).decode('cp1251').split(" ")
            if w in found_words:
                continue
            aln = hmm.viterbi({word:word, w:w})
            if aln[1] >= min_prob:
                r.append(word, " "+w.encode('cp1251'))            
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
def min_prob_calc (match, to_match, fine, gap, m_m, st, to_ex, ex, l_b):
    th = log(((match*to_match)**4)*to_match*m_m*(l_b**7)*st*to_ex*ex*to_match*fine*gap, 2)
    return float(int(th)-1)
###
def log_writer (s, si):
    w = open("aln_blocks.txt", 'ab')
    w.write(s.encode('cp1251')+"\t"+si+"\n")
    w.close()
############
###MAIN
#aln_min_prob = min_prob_calc(0.8, match_s_prob, gap_fine, gap_open, mismatch, start_prob, to_exit_prob, exit_prob, len_bonus)
def main_func (*params_list):
###parameters
    num_threads = int(params_list[0])
    in_f = params_list[1]
    draft_matrix_f = params_list[2]
    db_name = params_list[3]
    out_f = params_list[4]
    stops_list = params_list[5]
    aln_min_prob = params_list[6]
    mismatch = params_list[7]
    match_s_prob = params_list[8]
    match_to_match = params_list[9]
    av_match = params_list[10]
    gap_open = params_list[11]
    gap_ext = params_list[12]
    to_exit_prob = params_list[13]
    start_prob = params_list[14]
    ends_bonus = params_list[15]
    gap_fine = params_list[16]
    no_aln_fine = params_list[17]
    exit_prob = params_list[18]
    len_bonus = params_list[19]
    par_list = [match_s_prob, match_to_match, gap_open, gap_ext, to_exit_prob, start_prob, ends_bonus, gap_fine, no_aln_fine, exit_prob, len_bonus]
############
    wr = open("aln_blocks.txt", 'w')
    wr.close()
    r = redis.StrictRedis(host='localhost', port=6379, db=db_name)
    r.flushdb()
    p = subprocess.Popen("redis-server", shell=True)
    words_list = reader(in_f ,stops_list)
    d_m_read = open(draft_matrix_f, 'rb')
    draft_matrix = pickle.load(d_m_read) 
    log_emis = log_matrix(draft_matrix)
    log_params = get_params_log(*par_list)
    hmm = HMM(log_emis, *log_params)
    full_s = len(words_list)
    part_s = full_s/(num_threads-1) + 1
    part_list = []
    n = 1
    m = 0
    run_pars = []
    for word in words_list:
        part_list.append(word)
        if n == part_s:
            m = m + 1
            run_pars.append([part_list, words_list, hmm, db_name, aln_min_prob, m])
            part_list = []
            n = 1
            continue
        n = n + 1
    if part_list != []:
        m = m + 1
        run_pars.append([part_list, words_list, hmm, db_name, aln_min_prob, m])
        part_list = []
    pool = mp.Pool(num_threads)
    result = pool.map(get_aln_blocks, run_pars)
    pool.close()
    pool.join()
    r = redis.StrictRedis(host='localhost', port=6379, db=db_name)
    r_k = r.keys('*')
    for k in r_k:
        log_writer(r.get(k).decode('cp1251'), "")