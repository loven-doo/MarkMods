###imports
import pickle
import unicodedata
from probcons import probcons
############
###functions
def clust_to_dict (f):
    c_dict = {}
    read = open(f)
    for lines in read:
        line = lines.decode("cp1251").strip()
        if len(line) == 0:
            line = None
            continue
        line_list = line.split()
        c_dict[line_list[0]] = [line_list[0]]
        n = 1
        while n < len(line_list):
            c_dict[line_list[0]].append(line_list[n])
            n = n + 1
        line = None
    return c_dict
###
def join_clust (c_dict, l_size):
    cl_dict = {}
    c_k = c_dict.keys()
    joined_l = []
    for k in c_k:
        if k in joined_l:
            continue
        doub_l = []
        w_l = c_dict[k]
        n = 1
        sm_cl_l = []
        while n < len(w_l):
            links = links_calc(c_dict[w_l[n]], w_l)
            if len(w_l) == l_size:
                if links[0] == l_size:
                    sm_cl_l.append(w_l)
                    n = n + 1
                else:
                    sm_cl_l = []
                    break
            elif links[0] < 1 + l_size:
                a = w_l.pop(n)
                doub_l.append(a)
                a = None
            else:
                joined_l.append(w_l[n])
                for w in links[1]:
                    w_l.append(w)
                n = n + 1
        if sm_cl_l != []:
            for w in sm_cl_l:
                joined_l.append(w)
        if doub_l != []:
            w_l_ap = []
            n = 1
            while n < len(doub_l):
                links = links_calc(c_dict[doub_l[n]], w_l)
                if links[0] >= 1 + l_size:
                    joined_l.append(w_l[n])
                    for w in links[1]:
                        w_l_ap.append(w)
                n = n + 1
            for w in w_l_ap:
                w_l.append(w)
        cl_dict[k] = w_l
    return cl_dict
###
def links_calc (to_add_l, to_join_l):
    link = 0
    ap_l = []
    for w in to_add_l:
        if w in to_join_l:
            link = link + 1
        else:
            ap_l.append(w)
    return link, ap_l
###
def get_clusters (c_dict, matr, d, n_it, par_l, r_f):
    cl_dict = {}
    c_k = c_dict.keys()
    for k in c_k: 
        w_l = c_dict[k]
        if len(w_l) <= 2:#If I write profiles it will be more accurate
            cl_key = w_l[0]
        else:
            ws = read_words(w_l)
            m_aln = probcons(ws, matr, d, n_it, *par_l)
            cl_key = min_dist_sum(m_aln[0])
        cl_dict[cl_key] = w_l
    w_f = open(r_f, 'wb')
    pickle.dump(cl_dict, w_f)
###
def read_words (words):
    s_dict = {}
    for word in words:
        if len(word) == 0:
            continue
        prepared_line = line_prepare(word)
        s_dict[word] = prepared_line
    return s_dict
###
def line_prepare (l):
    s_l = []
    for s in l:
        s_l.append(s)
    return s_l
###
def min_dist_sum (dist_dict):
    dist_sum_d = {}
    d_k = dist_dict.keys()
    for k in d_k:
        ws = "".join(k)
        w_l = ws.split("-")
        if w_l[0] in dist_sum_d:
            dist_sum_d[w_l[0]] = dist_sum_d[w_l[0]] + dist_dict[k]
        else:
            dist_sum_d[w_l[0]] = dist_dict[k]
        if w_l[1] in dist_sum_d:
            dist_sum_d[w_l[1]] = dist_sum_d[w_l[1]] + dist_dict[k]
        else:
            dist_sum_d[w_l[1]] = dist_dict[k]       
        ws = None
#    print dist_sum_d
    ds_k = dist_sum_d.keys()
    key_w = ""
    min_dist = dist_sum_d[ds_k[0]]
    for k in ds_k:
        if dist_sum_d[k] < min_dist:
            key_w = None
            key_w = k
            min_dist = dist_sum_d[k]
    return key_w 
############
###MAIN
def main (params):
###parameters
    in_f = "aln_blocks.txt"
    out_f = "clusters.p"
    link_size = params["link_size"]
    matrix_f = "draft_matrix.p"
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
    par_list = [match_s_prob, match_to_match, gap_open, gap_ext, to_exit_prob, start_prob, ends_bonus, gap_fine, no_aln_fine, exit_prob, len_bonus]
############
    preclust_dict = clust_to_dict(in_f)
    clust_dict = join_clust(preclust_dict, link_size)
    m_read = open(matrix_f, 'rb')
    matrix = pickle.load(m_read)
    get_clusters(clust_dict, matrix, delta, n_iter, par_list, out_f)