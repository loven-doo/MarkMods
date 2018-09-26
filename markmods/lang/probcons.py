###imports
from hmm import HMM
from math import *
############
###functions
def probcons (sequences, emis, delta, n_iter, match, match_to_match, gap_open, gap_ext, to_exit_prob, start_prob, ends_bonus, gap_fine, no_aln_fine, exit_prob, len_bonus):
    log_emis = log_matrix(emis)
    log_params = get_params_log(match, match_to_match, gap_open, gap_ext, to_exit_prob, start_prob, ends_bonus, gap_fine, no_aln_fine, exit_prob, len_bonus)
    temp_f = params_calc(delta, n_iter)
    n = 1
    temp = temp_calc(temp_f, n)
    print ("Iteration "+str(n))
    hmm = HMM(log_emis, *log_params)
    seqs_pairs = get_seq_pairs(sequences)
    pair_alignments = get_pair_alignments(seqs_pairs, hmm)
    emis_new = emis_calc(pair_alignments[0], emis, temp)
    while n < n_iter:
        n = n + 1
        print ("Iteration "+str(n))
        temp = temp_calc(temp_f, n)
        emis = emis_new
        log_emis = log_matrix(emis)
        hmm = HMM(log_emis, *log_params)
        pair_alignments = get_pair_alignments(seqs_pairs, hmm)
        emis_new = emis_calc(pair_alignments[0], emis, temp)
#    test_print(pair_alignments[0])
#    print "hey"
#    mult_alignment = mult_from_pairs(pair_alignments[0], sequences)
#    return mult_alignment, pair_alignments[1], emis_new
    return pair_alignments[1], emis_new
###
def test_print (p_aln):
    w = open("test_alns.txt", 'wb')
    for aln in p_aln:
        for k in aln.keys():
            if k == "states":
                continue
            w.write(k.encode('cp1251')+"\t "+"".join(aln[k]).encode('cp1251')+"\n")
        w.write("\n")
    w.close()    
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
def get_pair_alignments (s_pairs, hmm):
    pair_alns = []
    probs = {}
    for pair in s_pairs:
        p_states_seq = hmm.viterbi(pair)
        pair_aln = build_alignment(p_states_seq[0], pair)
        pair_alns.append(pair_aln)
        aln_k = pair_aln.keys()
        probs[aln_k[0]+"-"+aln_k[1]] = p_states_seq[1]
    return pair_alns, probs
###
def build_alignment (st_seq, s_pair):
    k = s_pair.keys()
    seq1 = []
    seq2 = []
    for pos in st_seq:
        if pos[0] == "start":
            continue
        if pos[0] == "end":
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
    return {k[0]:seq1, k[1]:seq2, "states":st_seq}
###
def emis_calc (pair_alns, e_d, t):
    num_dict = {}
    subs_dict = {}
    for aln in pair_alns:
        k = keys_sort(aln.keys())
        l = len(aln[k[0]])
        n = 0
        ext = False
        while n < l:
            if aln[k[0]][n] == "-" or aln[k[1]][n] == "-":
                n = n + 1
                continue
            elif aln[k[0]][n] == aln[k[1]][n]:
                if aln[k[0]][n] in num_dict.keys():
                    num_dict[aln[k[0]][n]] = num_dict[aln[k[0]][n]] + 2
                else:
                    num_dict[aln[k[0]][n]] = 2
                if aln[k[0]][n] in subs_dict.keys():
                    if aln[k[1]][n] in subs_dict[aln[k[0]][n]].keys():
                        subs_dict[aln[k[0]][n]][aln[k[1]][n]] = subs_dict[aln[k[0]][n]][aln[k[1]][n]] + 2
                    else:
                        subs_dict[aln[k[0]][n]][aln[k[1]][n]] = 2
                else:
                    subs_dict[aln[k[0]][n]] = {aln[k[1]][n]: 2}
            else:
                if aln["states"][n] == "no_aln":
                    if aln["states"][n+1] == "match":
                        ext = True
                    if ext != True:
                        continue
                    ext = False
                else:
                    ext = True
                if aln[k[0]][n] in num_dict.keys():
                    num_dict[aln[k[0]][n]] = num_dict[aln[k[0]][n]] + 1
                else:
                    num_dict[aln[k[0]][n]] = 1
                if aln[k[1]][n] in num_dict.keys():
                    num_dict[aln[k[1]][n]] = num_dict[aln[k[1]][n]] + 1
                else:
                    num_dict[aln[k[1]][n]] = 1
                if aln[k[0]][n] in subs_dict.keys():
                    if aln[k[1]][n] in subs_dict[aln[k[0]][n]].keys():
                        subs_dict[aln[k[0]][n]][aln[k[1]][n]] = subs_dict[aln[k[0]][n]][aln[k[1]][n]] + 1
                    else:
                        subs_dict[aln[k[0]][n]][aln[k[1]][n]] = 1
                else:
                    subs_dict[aln[k[0]][n]] = {aln[k[1]][n]: 1}
                if aln[k[1]][n] in subs_dict.keys():
                    if aln[k[0]][n] in subs_dict[aln[k[1]][n]].keys():
                        subs_dict[aln[k[1]][n]][aln[k[0]][n]] = subs_dict[aln[k[1]][n]][aln[k[0]][n]] + 1                
                    else:
                        subs_dict[aln[k[1]][n]][aln[k[0]][n]] = 1
                else:
                    subs_dict[aln[k[1]][n]] = {aln[k[0]][n]: 1}                  
            n = n + 1
    emis_probs = subs_calc(subs_dict, num_dict, e_d, t)
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
def subs_calc (s_dict, n_dict, e, t):
    for key in e.keys():
        for key_s in e[key].keys():
            try:
                e[key][key_s] = (t*e[key][key_s]+float(s_dict[key][key_s])/float(n_dict[key]))/(t+1)
            except KeyError:
                e[key][key_s] = t*e[key][key_s]/(t+1)
    return e   
###
def mult_from_pairs (pair_alns, seqs):
    alns_dict = get_pairs_dict(pair_alns)
    s_k = seqs.keys()
    first_s = s_k[0]
    log_writer(first_s, "hey")
    mult_aln = prepare_seq(first_s, [0, len(seqs[first_s])+1], seqs[first_s])
    n_seqs = len(seqs)
    for i in range (1, n_seqs):#the correct way to go through is NP problem (komiv)
        log_writer(s_k[i], "hey")
        added_seq = add_seq(mult_aln, alns_dict, s_k[i], 0, len(seqs[s_k[i]]))
        mult_aln = added_seq[0]
    proc_aln = {}
    for k in s_k:
        proc_aln[k] = []
    return aln_processing(mult_aln, s_k, proc_aln)
###
def get_pairs_dict (p_alns):
    a_dict = {}
    for p in p_alns:
        p_k = keys_sort(p.keys())
        if p_k[0] in a_dict.keys():
            a_dict[p_k[0]][p_k[1]] = p
        else:
            a_dict[p_k[0]] = {p_k[1]:p}
        if p_k[1] in a_dict.keys():
            a_dict[p_k[1]][p_k[0]] = p
        else:
            a_dict[p_k[1]] = {p_k[0]:p}
    return a_dict
###
def prepare_seq (f_s, c, s):
    prep_seq = [f_s, c, [[]]]
    for let in s:
        prep_seq[2].append({f_s: let})
        prep_seq[2].append([])
    return prep_seq
###
def add_seq (m_aln, a_dict, s_key, s_l_n, l_s_l_n):#In future this code must be improved
    wrap_aln = [s_key, [s_l_n, s_l_n+1], [[]]]
    let_to_add = []
    lost_let = []
    s_l_n_c = True
    match = True
    s_l_n0 = s_l_n
    aln = a_dict[m_aln[0]][s_key]
    n = define_start(s_l_n, aln[s_key])
    l_n = define_lim(n, aln[m_aln[0]])
    while n < len(aln[m_aln[0]]):
        print l_n
        print wrap_aln
        if aln[s_key][n] != "-":
            s_l_n = s_l_n + 1
        if aln[m_aln[0]][n] != "-":
            l_n = l_n + 1
        if l_n in range (m_aln[1][0], m_aln[1][1]):
            if m_aln[1][0] == 0:
                pos = 2*l_n
            else:
                pos = 2*(l_n-m_aln[1][0]+1)
            if aln[m_aln[0]][n] == "-":
                match = False
                if m_aln[2][pos] == []:
                    wrap_aln[2].append({s_key:aln[s_key][n]})
                    wrap_aln[2].append([])
                else:
                    let_to_add.append(aln[s_key][n])
                    if s_l_n_c == True:
                        s_l_n0 = s_l_n-1
                        s_l_n_c = False
            else:
                if aln[s_key][n] == "-":
                    match = False
                    if let_to_add != []:
                        added_seq = add_seq(m_aln[2][pos-2], a_dict, s_key, s_l_n0, s_l_n)
                        wrap_aln = let_append(added_seq[0][2], wrap_aln)
                        s_l_n = added_seq[1]                 
                        n = define_start(s_l_n, aln[s_key])-1
                        let_to_add = []
                        s_l_n_c = True
                    s_pos = len(wrap_aln[2])-1
                    if wrap_aln[2][s_pos] == []:
                        wrap_aln[2][s_pos] = [m_aln[0], [l_n-1, l_n+1], [m_aln[2][pos-2]]]
                    elif wrap_aln[2][s_pos][0] != m_aln[0]:
                        wrap_aln[2][s_pos] = [m_aln[0], [l_n-1, l_n+1], [wrap_aln[2][s_pos]]]
                    wrap_aln[2][s_pos][2].append(m_aln[2][pos-1])
                    wrap_aln[2][s_pos][2].append(m_aln[2][pos])
                    wrap_aln[2][s_pos][1][1] = l_n+1
                else:
                    if let_to_add != []:
                        added_seq = add_seq(m_aln[2][pos-2], a_dict, s_key, s_l_n0, s_l_n-1)
                        wrap_aln = let_append(added_seq[0][2], wrap_aln)
                        s_l_n = added_seq[1]                 
                        n = define_start(s_l_n, aln[s_key])-1
                        let_to_add = []
                        s_l_n_c = True
                    if l_n == 1 & s_l_n == 1:
                        wrap_aln[2][0] = m_aln[2][pos-2]
                    if match == True:
                        wrap_aln[2][len(wrap_aln[2])-1] = m_aln[2][pos-2] 
                    wrap_aln[2].append(m_aln[2][pos-1])
                    wrap_aln[2][len(wrap_aln[2])-1][s_key] = aln[s_key][n]
                    wrap_aln[2].append([])
                    match = True
                    if n == len(aln[m_aln[0]])-1:
                         wrap_aln[2][len(wrap_aln[2])-1] = m_aln[2][pos]
            wrap_aln[1][1] = s_l_n + 1
            rem_s_l_n = s_l_n
        else:
            if let_to_add != []:
                if aln[s_key][n] != "-":
                    s_l_n = s_l_n - 1
                added_seq = add_seq(m_aln[2][pos], a_dict, s_key, s_l_n0, s_l_n)
                wrap_aln = let_append(added_seq[0][2], wrap_aln)                  
                s_l_n = added_seq[1]                                              
                n = define_start(s_l_n, aln[s_key])-1                             
                let_to_add = []                                                   
                s_l_n_c = True                                                    
            if l_n >= m_aln[1][1] & s_l_n <= l_s_l_n:#
                if lost_let == []:
                    lost_let = [[],{s_key:aln[s_key][n]},[]]
                else:
                    lost_let.append({s_key:aln[s_key][n]})
                    lost_let.append([])                   
                wrap_aln[1][1] = s_l_n + 1
                rem_s_l_n = s_l_n
        n = n + 1
    print wrap_aln
    if let_to_add != []:
        added_seq = add_seq(m_aln[2][pos], a_dict, s_key, s_l_n0, s_l_n)
        wrap_aln = let_append(added_seq[0][2], wrap_aln)
        s_l_n = added_seq[1]                 
        n = define_start(s_l_n, aln[s_key])-1
        let_to_add = []
        s_l_n_c = True
    if lost_let != []:
        wrap_aln = let_append(lost_let, wrap_aln)
        lost_let = []
    return wrap_aln, rem_s_l_n
###
def let_append (let_l, w_aln):
    l = len(w_aln[2])-1
    if w_aln[2][l] == []:
        a = w_aln[2].pop(l)
        a = None
    elif let_l[0] == []:
        a = let_l.pop(0)
        a = None
    else:
         joined_gap = join_gaps(let_l, w_aln[2])
         let_l = joined_gap[0]
         w_aln[2] = joined_gap[1]
    for e in let_l:
        w_aln[2].append(e)
    return w_aln
###
def join_gaps (l_l, w_a):
    l = len(w_a)-1
    if w_a[l] == []:
        a = l_l.pop(0)
        w_a[l].append(a)
        a = None
    else:
        joined_gap = join_gaps(l_l, w_a[l])
        l_l = joined_gap[0]
        w_a[l] = joined_gap[1]
    return l_l, w_a
###
def define_start (lim, seq):
    st = 0
    i = 0
    for s in seq:
        if st == lim:
            break
        if s != "-":
            st = st + 1
        i = i + 1
    return i
###
def define_lim (st, seq):
    lim = 0
    i = 0
    for s in seq:
        if i == st:
            break
        if s != "-":
            lim = lim + 1
        i = i + 1
    return lim 
###
def aln_processing (aln, s_k, proc_aln):
    n = 0
    cur_aln = aln[2]
    for elm in cur_aln:
        if elm == []:
            n = n + 1
            continue
        if int(n)/2 == float(n)/2.0:
            proc_aln = aln_processing(elm, s_k, proc_aln)
        else:
            for k in s_k:
                if k in elm.keys():
                    proc_aln[k].append(elm[k])
                else:
                    proc_aln[k].append("-")
        n = n + 1
    return proc_aln
###
def log_writer (s, si):
    w = open("log.txt", 'ab')
    w.write(s.encode('cp1251')+"\t"+si+"\n")
    w.close()