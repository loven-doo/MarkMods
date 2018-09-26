###imports
import pickle
from math import *
import unicodedata
import re
############
###functions
def freq_calc (f, stops_l, m_match, m_start):
    let_dict = {}
    let_num = 0
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
        next_string = calc_string(line_list, let_num, let_dict, stops_l)
        let_dict = next_string[0]
        let_num = next_string[1]
        line = None
    print let_dict, let_num
    l_freqs = num_to_freq(let_dict, let_num, m_match, m_start)
    return l_freqs
###
def calc_string (l_list, l_num, l_dict, s_l):
    splitter = get_splitter(s_l)
    words_l = []
    for elm in l_list:           
        word_l = re.split(splitter,elm)
        for word in word_l:
            if len(word.strip()) == 0:
                word = None
                continue
            words_l.append(word.strip())
    for word in words_l:
        next_word = calc_word(word, l_num, l_dict)
        l_dict = next_word[0]
        l_num = next_word[1]
        word = None 
    return l_dict, l_num
###
def get_splitter (spl_list):
    spl = []
    for s in spl_list:      
        spl.append(s)
    return "|".join(spl)
###
def calc_word (w, l_n, l_d):
    for let in w:
        l_n = l_n + 1
        if let in l_d.keys():
            l_d[let] = l_d[let] + 1
        else:
            l_d[let] = 1
    return l_d, l_n 
###
def num_to_freq (l_dict, l_num, m_m, m_s):
    l_d_k = l_dict.keys()
    for k in l_d_k:
        l_dict[k] = {k:m_s-(log10(l_dict[k])/log10(l_num))/2}
        for ki in l_d_k:
            if k == ki:
                continue
            l_dict[k][ki] = m_m
    return l_dict
###    
############
###MAIN
def main (params):
###parameters
    in_f = params["input_file"]
    out_f = "draft_matrix.p"
    stops_list = params["spec_symbs"]
    mismatch = params["mismatch"]#start missmatch value for substitutions matrix
    match_start = params["match_start"]#start value of match in matrix
############
    draft_matrix = freq_calc(in_f, stops_list, mismatch, match_start)
    w_f = open(out_f, 'wb')
    pickle.dump(draft_matrix, w_f)
