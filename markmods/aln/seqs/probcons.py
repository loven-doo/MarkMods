import subprocess
from math import *
import multiprocessing as mp
import pickle

import redis

from markmods.models.hmm.hmm import HMM
from markmods.models.hmm.hmm_profile import hmm_profile


class HubsTree:
    def __init__ (self, hubs, levels):
        self.hubs = hubs
        self.levels = levels


class Hub:
    def __init__ (self, content, links):
        self.content = content
        self.links = links


def probcons(sequences, emis, threads, delta, n_iter, match, match_to_match, gap_open, gap_ext, to_exit_prob, start_prob, ends_bonus, gap_fine, no_aln_fine, exit_prob, prof_match_c, len_bonus, dist_step):
    symb_freqs = freqs_calc(sequences)
    log_emis = log_matrix(emis)
    log_params = get_params_log(match, match_to_match, gap_open, gap_ext, to_exit_prob, start_prob, ends_bonus, gap_fine, no_aln_fine, exit_prob)
    l_bonus = log(len_bonus)
    p_match_c = log(prof_match_c)
    temp_f = params_calc(delta, n_iter)
    n = 1
    temp = temp_calc(temp_f, n)
    subprocess.call("echo 1 > /proc/sys/net/ipv4/tcp_tw_reuse", shell=True)
    subprocess.call("echo 1 > /proc/sys/net/ipv4/tcp_tw_recycle", shell=True)
    subprocess.call("redis-cli shutdown 0>log.txt 1>>log.txt 2>err.txt", shell=True)
    p = subprocess.Popen("redis-server 0>>log.txt 1>>log.txt 2>>err.txt", shell=True)
    print("Distances calculating")
    print("Iteration "+str(n))
    hmm = HMM(log_emis, *log_params)
    seqs_pairs = get_seq_pairs(sequences)
    pair_alignments = get_pair_alignments(seqs_pairs, hmm, l_bonus, threads)
    emis_new = emis_calc(pair_alignments[0], emis, symb_freqs, temp)
    while n < n_iter:
        n = n + 1
        print ("Iteration "+str(n))
        temp = temp_calc(temp_f, n)
        emis = emis_new
        log_emis = log_matrix(emis)
        hmm = HMM(log_emis, *log_params)
        pair_alignments = get_pair_alignments(seqs_pairs, hmm, l_bonus, threads)
        emis_new = emis_calc(pair_alignments[0], emis, symb_freqs, temp)
#    matrix_writer("res_matr.txt", log_emis)#
    mult_alignment = mult_from_pairs(pair_alignments, sequences, hmm, p_match_c, l_bonus, dist_step, threads)
    subprocess.call("redis-cli shutdown", shell=True)
    return mult_alignment[0], mult_alignment[1], pair_alignments[1], emis_new
###
def freqs_calc (ss):
    s_nums = {}
    s_num = 0
    s_k = ss.keys()
    for k in s_k:
        for s in ss[k]:
            try:
                s_nums[s] = s_nums[s]+1
            except KeyError:
                s_nums[s] = 1
            s_num = s_num+1
    s_nums = num_to_freq(s_nums, s_num)
    s_nums = pairs_freq(s_nums)
    return s_nums
###
def num_to_freq (nums_d, num):
    n_k = nums_d.keys()
    for k in n_k:
        nums_d[k] = float(nums_d[k])/float(num)
    return nums_d
###
def pairs_freq (freqs_d):
    pairs_fr = {}
    f_k = freqs_d.keys()
    for ki in f_k:
        for kj in f_k:
            pairs_fr[frozenset([ki,kj])] = freqs_d[ki]*freqs_d[kj]
    return pairs_fr
###
def test_print (f, struc):
    w = open(f, 'w')
    w.write(pickle.dumps(struc))
    w.close()    
###
def log_matrix (matr):
    log_matr = {}
    m_k = matr.keys()
    for k in m_k:
        log_matr[k] = log(matr[k])
    return log_matr
###
def get_params_log (*par_l):
    log_par = []
    for p in par_l:
        log_par.append(log(p))
    return log_par
###
def params_calc (d, it):
    k = 1.0
    c = 1.0/d - (k*sqrt(it)+1.0)
    while c > 5.0:
        k = k + 1.0
        c = 1.0/d - (k*sqrt(it)+1.0)
    return c, k
###
def temp_calc (par, i):
    t = par[1]*sqrt(i)+par[0]
    return t
###
def get_seq_pairs (seqs):
    i = 0
    keys_l = seqs.keys()
    l = len(keys_l)
    s_pairs = []
    while i < l-1:
        j = i + 1
        while j < l:
            s_pairs.append({keys_l[i]:seqs[keys_l[i]], keys_l[j]:seqs[keys_l[j]]})
            j = j + 1
        i = i + 1
    return s_pairs
###
def get_pair_alignments (s_pairs, hmm, l_bon, n_thr):
    pair_alns = {}
    probs = {}
    run_list = []
    for pair in s_pairs:
        run_list.append([pair, hmm, l_bon, 0, 1])
    pool = mp.Pool(n_thr)
    result = pool.map(align_pair, run_list)
    pool.close()
    pool.join()
    r = redis.StrictRedis(host='localhost', db=0)
    r_k = r.keys('*')
    for k in r_k:
        pair_alns[pickle.loads(k)] = pickle.loads(r.get(k))
    r.flushdb()
    r = redis.StrictRedis(host='localhost', db=1)
    r_k = r.keys('*')
    for k in r_k:
        probs[pickle.loads(k)] = pickle.loads(r.get(k))
    r.flushdb()    
    return pair_alns, probs
###
def align_pair (run_list):
    pair = run_list[0]
    hmm = run_list[1]
    l_b = run_list[2]
    db1 = run_list[3]
    db2 = run_list[4]
    p_states_seq = hmm.viterbi(pair, l_b)
    pair_aln = build_alignment(p_states_seq[0], pair)
    aln_k = keys_sort(pair_aln.keys())    
    p_aln_k = pickle.dumps(frozenset([aln_k[0],aln_k[1]]))
    r = redis.StrictRedis(host='localhost', db=db1)
    r.append(p_aln_k, pickle.dumps(pair_aln))
    r = redis.StrictRedis(host='localhost', db=db2)
    r.append(p_aln_k, pickle.dumps(float(p_states_seq[1])))
###
def build_alignment (st_seq, s_pair):
    k = s_pair.keys()
    seq1 = []
    seq2 = []
    sts = []
    for pos in st_seq:
        if pos[0] == "start":
            continue
        if pos[0] == "end":
            sts.append(pos[0])
            break
        if pos[0] == "match":
            seq1.append(s_pair[k[0]][pos[2]-1])
            seq2.append(s_pair[k[1]][pos[1]-1])
        if pos[0] == "no_aln":
            seq1.append(s_pair[k[0]][pos[2]-1])
            seq2.append(s_pair[k[1]][pos[1]-1])
        if pos[0] == "insertion_1":
            seq1.append(s_pair[k[0]][pos[2]-1])
            seq2.append("-")        
        if pos[0] == "insertion_2":
            seq1.append("-")
            seq2.append(s_pair[k[1]][pos[1]-1])
        sts.append(pos[0])
    return {k[0]:seq1, k[1]:seq2, "states":sts}
###
def emis_calc (pair_alns, e_d, symb_fr, t):
    num_dict = {}
    subs_dict = {}
    pa_k = pair_alns.keys()
    for k in pa_k:
        aln = pair_alns[k]
        k = keys_sort(aln.keys())
        l = len(aln[k[0]])
        n = 0
        ext = False
        while n < l:
            if aln[k[0]][n] == "-" or aln[k[1]][n] == "-":
                n = n + 1
                continue
            else:
                if aln["states"][n] == "no_aln":
                    if aln["states"][n+1] == "match":
                        ext = True
                    if ext != True:
                        n = n + 1
                        continue
                    ext = False
                else:
                    ext = True
                pair_key = frozenset([aln[k[0]][n],aln[k[1]][n]])
                try:
                    num_dict[pair_key] = num_dict[pair_key] + 2
                except KeyError:
                    num_dict[pair_key] = 2
                try:
                    subs_dict[pair_key] = subs_dict[pair_key] + 2
                except KeyError:
                    subs_dict[pair_key] = 2
                pair_key = None    
            n = n + 1
    emis_probs = subs_calc(subs_dict, symb_fr, num_dict, e_d, t)
    return emis_probs
###
def keys_sort (k_l):
    n = 0
    for k in k_l:
       if k == "states":
           a = k_l.pop(n)
           k_l.append(k)
           return k_l
       n = n + 1
###
def subs_calc (s_dict, fr_dict, n_dict, e, t):
    for key in e.keys():
        try:
            e[key] = (t*e[key]+(float(s_dict[key])/fr_dict[key])/float(n_dict[key]))/(t+1)
        except KeyError:
            e[key] = t*e[key]/(t+1)
    return e
###
def sum_matrix (m):
#the matrix is dictionary!
    s = 0.0
    m_k = m.keys()
    for k in m_k:
        s = s + m[k]
    return s
###
def mult_from_pairs (pair_alns, seqs, hmm, match_c, l_bon, d_step, n_thr):
    s_k = seqs.keys()
    print "Guide tree building" 
#    order_list = komiv(pair_alns[1], s_k[0])#the correct way to go through is NP problem (komiv)
    order_list = s_k#
    hubs_tree = select_clust(order_list, pair_alns[1], d_step)#should be inserted to komiv
    lev_num = len(hubs_tree.levels)
    n = -1
    print "Multiple alignment building"
    while n+lev_num >= 0:
        print "level %s" % abs(n)
        run_list = configure_run(hubs_tree.levels[n], hubs_tree.hubs, seqs, hmm, match_c, l_bon, 0)
        pool = mp.Pool(n_thr)
        result = pool.map(hubs_aln, run_list)
        pool.close()
        pool.join()
        hubs_tree.hubs = read_db(0, hubs_tree.hubs)
        n = n - 1
    aln_key = hubs_tree.levels[0][0]
    mult_aln = hubs_tree.hubs[aln_key].content
    mult_aln_profile = hmm_profile(mult_aln)
    return mult_aln, mult_aln_profile
###
def configure_run (levs, hs, ss, hmm, m_c, l_b, db):
    r_list = []
    for k in levs:
        r_list.append([k, hs, ss, hmm, m_c, l_b, db])
    return r_list
###
def hubs_aln (r_list):
    hub_k = r_list[0]
    hubs = r_list[1]
    seqs = r_list[2]
    hmm = r_list[3]
    match_c = r_list[4]
    len_bon = r_list[5]
    db = r_list[6]
    h_seqs = hubs[hub_k].content
    if len(h_seqs) == 1 or str(type(h_seqs)) == "<type 'dict'>":
        return None
    if str(type(hubs[h_seqs[0]].content)) == "<type 'dict'>":
        aln1 = hubs[h_seqs[0]].content
    else:
        aln1 = {h_seqs[0]:seqs[h_seqs[0]]}
    l_h_seqs = len(h_seqs)
    for i in range (1, l_h_seqs):
        if str(type(hubs[h_seqs[i]].content)) == "<type 'dict'>":
            aln2 = hubs[h_seqs[i]].content
        else:
            aln2 = {h_seqs[i]:seqs[h_seqs[i]]}
        aln1_prof = hmm_profile(aln1)
        aln2_prof = hmm_profile(aln2)
        aln1_states = hmm.profile_viterbi(aln1_prof, len(aln1.keys()), aln2_prof, len(aln2.keys()), match_c, len_bon)
        aln1 = build_mult_aln(aln1_states[0], aln1, aln2)
    r = redis.StrictRedis(host='localhost', db=db)
    r.append(pickle.dumps(hub_k), pickle.dumps(aln1))
###
def select_clust (ord_l, dist, d_st):
    mi_dist = min_dist(dist)
    rmi_dist = round(mi_dist, 0)#this works only with d_st = -0.5
    if rmi_dist > mi_dist:
        m_dist = rmi_dist
    else:
        m_dist = rmi_dist + 0.5
    splits = [m_dist]
    clust_split = m_dist+d_st
    clust_num = len(ord_l)
    cl_num_f = {m_dist:clust_num}
    while clust_num > 1:
        clust_num = cl_num_calc(ord_l, dist, clust_split)
        cl_num_f[clust_split] = clust_num
        splits.append(clust_split)
        clust_split = clust_split + d_st
    split_dists = define_split_dist(cl_num_f, splits)
    print "%s levels in tree" % len(split_dists)
    ord_l = cor_start(ord_l, dist, split_dists[0])
    lev0 = frozenset(ord_l)
    return define_hubs(lev0, dist, split_dists, 1, HubsTree({lev0:Hub(ord_l, (ord_l[0], ord_l[-1]))}, [[lev0]]))
###
def min_dist (d):
    d_k = d.keys()
    min_d = d[d_k[0]]
    for k in d_k:
        if d[k] > min_d:#depends on distances signs
            min_d = d[k]
    return min_d
###
def cl_num_calc (o_l, d, cl_spl):
    cl_n = 0
    n = 0
    l = len(o_l)
    while n < l:
        try:
            if d[frozenset([o_l[n],o_l[n+1]])] < cl_spl: #depends on distances signs
                cl_n = cl_n + 1
        except IndexError:
            if d[frozenset([o_l[n],o_l[0]])] < cl_spl:
                cl_n = cl_n + 1            
        n = n + 1
    if cl_n == 0:
        cl_n = 1
    return cl_n
###
def define_split_dist (f, x):
    spls = []
    f_1 = num_deriv(f, x)
    f_2 = num_deriv(f_1, x)
    loc_max = define_loc_max(f_2, x)
    for xi in x:
        if f_2[xi] > 0.0:
            if xi in loc_max:
                spls.insert(0, xi)
    spls.insert(0, x[-1])
    return spls
###
def define_loc_max (f, x):
    max_l = []
    for n in range(1, len(x)-1):
        if f[x[n]] > f[x[n-1]] and f[x[n]] > f[x[n+1]]:
            max_l.append(x[n])
    return max_l
###
def func_par_wr (f_name, f, x):
    w = open(f_name+"y.txt", 'w')
    w.write(pickle.dumps(f))
    w.close()
    w = open(f_name+"x.txt", 'w')
    w.write(pickle.dumps(x))
    w.close()    
###
def num_deriv (f, x):
    deriv = {x[0]:float(f[x[1]]-f[x[0]])/float(x[1]-x[0])}
    n = 1
    l_x = len(x)
    while n < l_x:
        try:
            deriv[x[n]] = float(f[x[n+1]]-f[x[n-1]])/float(x[n+1]-x[n-1])
        except IndexError:
            deriv[x[n]] = float(f[x[n]]-f[x[n-1]])/float(x[n]-x[n-1])
        n = n + 1
    return deriv
###
def cor_start (o_l, d, spl):
    n = 0
    l = len(o_l)
    while l+n > 0:
        if d[frozenset([o_l[0],o_l[-1]])] >= spl:
            a = o_l.pop(-1)
            o_l.insert(0, a)
            a = None
        else:
            break
        n = n - 1
    return o_l
###
def define_hubs (h_key, d, spls, lev, h_tree):
    o_l = h_tree.hubs[h_key].content
    try:
        cl_spl = spls[lev]
    except IndexError:
        new_hubs = get_hubs_keys(o_l, h_tree.hubs)
        branch_hubs_k = new_hubs[0]
        h_tree.hubs = new_hubs[1]
        return hubs_sort(h_key, h_tree, branch_hubs_k, d)
    branch_hubs = [[o_l[0]]]
    n = 0
    l = len(o_l)
    while n < l-1:
        if d[frozenset([o_l[n],o_l[n+1]])] < cl_spl:
            branch_hubs.append([o_l[n+1]])
        else:
            branch_hubs[-1].append(o_l[n+1])
        n = n + 1
    new_hubs = get_hubs_keys(branch_hubs, h_tree.hubs)
    branch_hubs_k = new_hubs[0]
    h_tree.hubs = new_hubs[1]
    try:
        h_tree.levels[lev] = h_tree.levels[lev]+branch_hubs_k
    except IndexError:
        h_tree.levels.append(branch_hubs_k)
    for i in range(0,len(branch_hubs_k)):
        h_tree = define_hubs(branch_hubs_k[i], d, spls, lev+1, h_tree)
    return hubs_sort(h_key, h_tree, branch_hubs_k, d)
###
def get_hubs_keys (b_hs, hs):
    b_hs_k = []
    for h in b_hs:
        if str(type(h)) == "<type 'list'>":
            if len(h) > 1:
                h_k = frozenset(h)
                hs[h_k] = Hub(h, (h[0], h[-1]))
            else:
                h_k = h[0]
                hs[h_k] = Hub(h, (h[0], h[0]))
        else:
            h_k = h
            hs[h_k] = Hub([h], (h, h))
        b_hs_k.append(h_k)
        h_k = None
    return b_hs_k, hs
###
def hubs_sort (h_k, h_t, b_hk, d):
    h = h_t.hubs
    l = len(b_hk)
    if l == 1:
        return h_t
    min_d_k = frozenset([h[b_hk[0]].links[-1], h[b_hk[1]].links[0]])
    min_d = d[min_d_k]
    min_e = (0, 1)
    for n in range(1, l-1):
        k = frozenset([h[b_hk[n]].links[-1], h[b_hk[n+1]].links[0]]) 
        if d[k] > min_d:            
            min_d_k = None
            min_d_k = k
            min_d = None
            min_d = d[k]
            min_e = None
            min_e = (n, n+1)
        k = None
    h_t.hubs[h_k].content = get_hub_order(h, b_hk, min_e, d)
    return h_t
###
def get_hub_order (h, b_hk, st_pair, d):
    to_prev = None
    to_next = None
    s_h = []
    for p in st_pair:
        s_h.append(b_hk[p])
    n = st_pair[1]
    m = st_pair[0]
    l_h = len(b_hk)-1
    while m > 0 or n < l_h:
        if m > 0:
            to_prev = d[frozenset([h[b_hk[m]].links[0], h[b_hk[m-1]].links[-1]])]
        if n < l_h:
            to_next = d[frozenset([h[b_hk[n]].links[-1], h[b_hk[n+1]].links[0]])]
        if to_prev == None and to_next == None:
            break
        elif to_prev > to_next or to_next == None:
            s_h.append(b_hk[m-1])
            to_prev = None
            m = m - 1
        else:
            s_h.append(b_hk[n+1])
            to_next = None
            n = n + 1
    return s_h        
###
def select_base (p_alns, o_list, ss):
    p_alns_l = []
    n = 0
    l_o_list = len(o_list)
    while n < l_o_list:
        if l_o_list-n == 1:
            p_alns_l.append({o_list[n]:ss[o_list[n]]})
        else:
            p_alns_l.append(p_alns[frozenset([o_list[n], o_list[n+1]])])
        n = n + 2
    return p_alns_l
###
def read_db (db, hubs_d):
    r = redis.StrictRedis(host='localhost', db=db)
    r_k = r.keys('*')
    for k in r_k:
        hubs_d[pickle.loads(k)].content = pickle.loads(r.get(k))
    r.flushdb()
    return hubs_d    
###
def build_mult_aln (st_seq, aln1, aln2):
    aln1_k = aln1.keys()
    aln2_k = aln2.keys()
    aln = aln_prepare(aln1_k, aln2_k)
    for pos in st_seq:
        if pos[0] == "start":
            continue
        if pos[0] == "end":
            break
        if pos[0] == "match":
            aln = append_let(aln, aln1, aln1_k, pos[2]-1)
            aln = append_let(aln, aln2, aln2_k, pos[1]-1)
        if pos[0] == "no_aln":
            aln = append_let(aln, aln1, aln1_k, pos[2]-1)
            aln = append_let(aln, aln2, aln2_k, pos[1]-1)
        if pos[0] == "insertion_1":
            aln = append_let(aln, aln1, aln1_k, pos[2]-1)
            aln = append_gap(aln, aln2, aln2_k)
        if pos[0] == "insertion_2":
            aln = append_gap(aln, aln1, aln1_k)
            aln = append_let(aln, aln2, aln2_k, pos[1]-1)
    return aln
###
def aln_prepare (a1_k, a2_k):
    aln = {}
    for k in a1_k:
         aln[k] = []
    for k in a2_k:
         aln[k] = []
    return aln
###
def append_let (aln, a, a_k, p):
    for k in a_k:
        aln[k].append(a[k][p])
    return aln
###
def append_gap (aln, a, a_k):
    for k in a_k:
        aln[k].append("-")
    return aln
###
def log_writer (s, si):
    w = open("log.txt", 'ab')
    w.write(s.encode('cp1251')+"\t"+si+"\n")
    w.close()
###
def matrix_writer (f, matr):
    m_dict = {}
    m_k = matr.keys()
    for k in m_k:
        f_e = None
        for e in k:
            if len(k) == 1:
                try:
                    m_dict[e][e] = matr[k]
                except KeyError:
                    m_dict[e] = {e:matr[k]}
            elif f_e == None:
                if e not in m_dict.keys():
                    m_dict[e] = {}
                f_e = e
            else:
                m_dict[f_e][e] = matr[k]
                f_e = None
    m_d_k = m_dict.keys()
    m_d_k.sort()
    for k in m_d_k:
        k_k = m_dict[k].keys()
        k_k.sort()
        for ki in k_k:
            writer(f, k+"-"+ki+":"+str(m_dict[k][ki]))
###
def writer (f, s):
    w = open (f, "a")
    w.write(s + "\n")
    w.close()
