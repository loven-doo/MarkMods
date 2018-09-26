# This HMMis oly for sequences alignments => needs to be refactored
from math import *


class HMM(object):

    def __init__(self, match_emis, match, match_to_match, gap_open, gap_extantion, to_exit_prob, start_prob, e_bon, gap_fine, no_aln_fine, exit_prob):
        self.start = Start(start_prob, match, gap_open, to_exit_prob, e_bon)
        self.match = Match(match_emis, match+match_to_match, gap_open, to_exit_prob)
        self.insert_1 = Insertion("1", gap_fine, gap_extantion, match, gap_open, to_exit_prob, e_bon)
        self.insert_2 = Insertion("2", gap_fine, gap_extantion, match, gap_open, to_exit_prob, e_bon)
        self.no_aln = NoAlignment(no_aln_fine, gap_extantion, match, gap_open, to_exit_prob, e_bon)
        self.end = End(exit_prob)

    def viterbi(self, seq_pair, len_bonus):
        k = seq_pair.keys()
        seq1 = seq_pair[k[0]]
        seq2 = seq_pair[k[1]]
        l1 = len(seq1)
        l2 = len(seq2)+1
        aln_matrix = [[(self.start.start_prob, self.start, None)]]
        i = 0                         
        while i <= l2:
            aln_matrix.append([])
            j = 0
            while j <= l1:
                if i == l2:
                    j = l1
                    p = aln_matrix[i-1][j][0]+aln_matrix[i-1][j][1].to_exit_prob+self.end.exit_prob
                    aln_matrix[i].append((p, self.end, (i-1, j)))                                                        
                    j = j + 1
                    continue
                elif i == 0:
                    if j != 0:
                        p = aln_matrix[i][j-1][0]+aln_matrix[i][j-1][1].gap_open+self.insert_1.gap_fine+len_bonus
                        aln_matrix[i].append((p, self.insert_1, (i,j-1))) 
                    j = j + 1
                    continue 
                elif j == 0:
                    p = aln_matrix[i-1][j][0]+aln_matrix[i-1][j][1].gap_open+self.insert_2.gap_fine+len_bonus
                    aln_matrix[i].append((p, self.insert_2, (i-1,j)))
                    j = j + 1
                    continue
                aln_matrix[i].append(self.aln_state(i, j, seq1, seq2, len_bonus, aln_matrix))
                j = j + 1
            i = i + 1
        states_seq = self.get_seq(aln_matrix)
        return states_seq

    def aln_state(self, i, j, seq1, seq2, l_bonus, aln_matrix):
        try:
            p_match = aln_matrix[i-1][j-1][0]+aln_matrix[i-1][j-1][1].match+self.match.match_emis[frozenset([seq1[j-1],seq2[i-1]])]+l_bonus
        except KeyError:
            p_match = l_bonus+self.match.min_emis-1
        p_gap2 = aln_matrix[i][j-1][0]+aln_matrix[i][j-1][1].gap_open+self.insert_1.gap_fine+l_bonus
        p_gap1 = aln_matrix[i-1][j][0]+aln_matrix[i-1][j][1].gap_open+self.insert_2.gap_fine+l_bonus
        p_no_aln = aln_matrix[i-1][j-1][0]+aln_matrix[i-1][j-1][1].no_aln_open+self.no_aln.no_aln_fine+l_bonus
        p = max(p_match, p_gap2, p_gap1, p_no_aln)
        if p == p_match:
            aln_s = (p, self.match, (i-1,j-1))
        elif p == p_gap2:
            aln_s = (p, self.insert_1, (i,j-1)) 
        elif p == p_gap1:
            aln_s = (p, self.insert_2, (i-1,j))
        elif p == p_no_aln:
            aln_s = (p, self.no_aln, (i-1,j-1))
        return aln_s    

    def get_seq(self, l):
        i = len(l)-2
        j = 0
        seq = []
        st = l[i][j]
        prob = st[0]
        while st[1] != self.start:
            seq.insert(0, (st[1].name, i, j))
            i = st[2][0]
            j = st[2][1]
            st = None
            st = l[i][j]
        seq.insert(0, (st[1].name, i, j))            
        return seq, prob

    def profile_viterbi(self, prof1, num_s1, prof2, num_s2, match_c, len_bonus):
        l1 = len(prof1)
        l2 = len(prof2)+1
        aln_matrix  = [[(self.start.start_prob, self.start, None)]]
        i = 0                         
        while i <= l2:
            aln_matrix.append([])
            j = 0
            while j <= l1:
                if i == l2:
                    j = l1
                    p = aln_matrix[i-1][j][0]+aln_matrix[i-1][j][1].to_exit_prob+self.end.exit_prob
                    aln_matrix[i].append((p, self.end, (i-1, j)))                                                        
                    j = j + 1
                    continue
                elif i == 0:
                    if j != 0:
                        p = aln_matrix[i][j-1][0]+aln_matrix[i][j-1][1].gap_open+self.insert_1.gap_fine+len_bonus
                        aln_matrix[i].append((p, self.insert_1, (i,j-1))) 
                    j = j + 1
                    continue 
                elif j == 0:
                    p = aln_matrix[i-1][j][0]+aln_matrix[i-1][j][1].gap_open+self.insert_2.gap_fine+len_bonus
                    aln_matrix[i].append((p, self.insert_2, (i-1,j)))
                    j = j + 1
                    continue
                aln_matrix[i].append(self.prof_aln_state(i, j, prof1, num_s1, prof2, num_s2, match_c, len_bonus, aln_matrix))
                j = j + 1
            i = i + 1
        states_seq = self.get_seq(aln_matrix)
        return states_seq

    def prof_aln_state(self, i, j, seq1, n_s1, seq2, n_s2, m_c, l_bonus, aln_matrix):
        p_match = aln_matrix[i-1][j-1][0]+aln_matrix[i-1][j-1][1].match+self.match.profiles_match(seq1[j-1],seq2[i-1])-(m_c+log(n_s1)+log(n_s2))+l_bonus
        p_gap2 = aln_matrix[i][j-1][0]+aln_matrix[i][j-1][1].gap_open+self.insert_1.gap_fine+l_bonus
        p_gap1 = aln_matrix[i-1][j][0]+aln_matrix[i-1][j][1].gap_open+self.insert_2.gap_fine+l_bonus
        p_no_aln = aln_matrix[i-1][j-1][0]+aln_matrix[i-1][j-1][1].no_aln_open+self.no_aln.no_aln_fine+l_bonus
        p = max(p_match, p_gap2, p_gap1, p_no_aln)
        if p == p_match:
            aln_s = (p, self.match, (i-1,j-1))
        elif p == p_gap2:
            aln_s = (p, self.insert_1, (i,j-1)) 
        elif p == p_gap1:
            aln_s = (p, self.insert_2, (i-1,j))
        elif p == p_no_aln:
            aln_s = (p, self.no_aln, (i-1,j-1))
        return aln_s    
#    def forward_backward (self, seq_pair):

#        return states_prob

#    def forward (self, seq_pair):
    
#    def backward (self, seq_pair):


class Start(object):

    def __init__(self, start_prob, match, gap_open, to_exit_prob, e_bon):
        self.name = "start"
        self.start_prob = start_prob
        self.match =  match
        self.gap_open = gap_open+e_bon
        self.no_aln_open = gap_open+e_bon
        self.to_exit_prob = to_exit_prob 


class Match(object):

    def __init__(self, match_emis, match, gap_open, to_exit_prob):
        self.name = "match"
        self.match_emis = match_emis
        self.match = match
        self.gap_open = gap_open
        self.no_aln_open = gap_open
        self.to_exit_prob = to_exit_prob
        self.min_emis = self.find_min(match_emis)

    def find_min(self, matr):
        m_k = matr.keys()
        min_v = matr[m_k[0]]
        for k in m_k:
            v = matr[k]
            if min_v > v:
                min_v = v
        return min_v    

    def profiles_match(self, prof1_col, prof2_col):
        s1_k = prof1_col.keys()
        s2_k = prof2_col.keys()
        score = 0
        for k1 in s1_k:
            for k2 in s2_k:
                try:
                    score = score + prof1_col[k1]*prof2_col[k2]*self.match_emis[frozenset([k1,k2])]
                except KeyError:
                    score = score + self.min_emis-1
        return score


class Insertion(object):

    def __init__(self, name, gap_fine, gap_extantion, match, gap_open, to_exit_prob, e_bon):
        self.name = "insertion_"+name
        self.gap_fine = gap_fine
        self.match = match
        self.gap_open = gap_extantion
        self.no_aln_open = gap_extantion
        self.to_exit_prob = to_exit_prob+e_bon


class NoAlignment(object):

    def __init__(self, no_aln_fine, gap_extantion, match, gap_open, to_exit_prob, e_bon):
        self.name = "no_aln"
        self.no_aln_fine = no_aln_fine
        self.match = match
        self.gap_open = gap_open
        self.no_aln_open = gap_extantion
        self.to_exit_prob = to_exit_prob+e_bon


class End(object):

    def __init__(self, exit_prob):
        self.name = "end"
        self.exit_prob = exit_prob
