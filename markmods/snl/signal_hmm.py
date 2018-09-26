class SignalHMM:
    def __init__ (self, s1_e, s2_e, s3_e, s1_s1, s2_s2, s3_s3, s1_s2, s1_s3, s2_s3, entr_s1, entr_s2, entr_s3,
                  s1_exit, s2_exit, s3_exit):
        self.s1 = State(emis=s1_e, _s1=s1_s1, _s2=s1_s2, _s3=s1_s3, _exit=s1_exit)
        self.s2 = State(emis=s2_e, _s1=s1_s2, _s2=s2_s2, _s3=s2_s3, _exit=s2_exit)
        self.s3 = State(emis=s3_e, _s1=s1_s3, _s2=s2_s3, _s3=s3_s3, _exit=s3_exit)
        self.entrance = Entrance(_s1=entr_s1, _s2=entr_s2, _s3=entr_s3)

    def prepare_signal(self, signal_list):
        prep_signal_list = []
        forward_list = self.forward(signal_list)
        backward_list = self.backward(signal_list)
        for n in range(0,len(signal_list)):
            prep_signal_list.append(self._choose_state(forward_list[n], backward_list[-(n+2)], forward_list[-1]))
        return prep_signal_list

    def forward(self, signal_l):
        forward_l = []
        scores_dict = {}
        scores_dict[self.s1.name] = self.s1.emis.get_emis(signal_l[0])*self.entrance.to_s1
        scores_dict[self.s2.name] = self.s2.emis.get_emis(signal_l[0])*self.entrance.to_s2
        scores_dict[self.s3.name] = self.s3.emis.get_emis(signal_l[0])*self.entrance.to_s3
        forward_l.append(scores_dict)
        for n in range (1, len(signal_l)):
            scores_dict = {}
            scores_dict[self.s1.name] = self.s1.emis.get_emis(signal_l[n])*(forward_l[-1][self.s1.name]*self.s1.to_s1 +
                                                                   forward_l[-1][self.s2.name]*self.s2.to_s1 +
                                                                   forward_l[-1][self.s3.name]*self.s3.to_s1)
            scores_dict[self.s2.name] = self.s2.emis.get_emis(signal_l[n])*(forward_l[-1][self.s1.name]*self.s1.to_s2 +
                                                                   forward_l[-1][self.s2.name]*self.s2.to_s2 +
                                                                   forward_l[-1][self.s3.name]*self.s3.to_s2)
            scores_dict[self.s3.name] = self.s3.emis.get_emis(signal_l[n])*(forward_l[-1][self.s1.name]*self.s1.to_s3 +
                                                                   forward_l[-1][self.s2.name]*self.s2.to_s3 +
                                                                   forward_l[-1][self.s3.name]*self.s3.to_s3)
            forward_l.append(scores_dict)
        forward_l.append(forward_l[-1][self.s1.name]*self.s1.to_exit +
                         forward_l[-1][self.s2.name]*self.s2.to_exit +
                         forward_l[-1][self.s3.name]*self.s3.to_exit)
        return forward_l

    def backward(self, signal_l):
        backward_l = []
        scores_dict = {}
        scores_dict[self.s1.name] = self.s1.to_exit
        scores_dict[self.s2.name] = self.s2.to_exit
        scores_dict[self.s3.name] = self.s3.to_exit
        backward_l.append(scores_dict)
        for n in range (1, len(signal_l)):
            scores_dict = {}
            scores_dict[self.s1.name] = (self.s1.emis.get_emis(signal_l[-n])*backward_l[-1][self.s1.name]*self.s1.to_s1 +
                                         self.s2.emis.get_emis(signal_l[-n])*backward_l[-1][self.s2.name]*self.s1.to_s2 +
                                         self.s3.emis.get_emis(signal_l[-n])*backward_l[-1][self.s3.name]*self.s1.to_s3)
            scores_dict[self.s2.name] = (self.s1.emis.get_emis(signal_l[-n])*backward_l[-1][self.s1.name]*self.s2.to_s1 +
                                         self.s2.emis.get_emis(signal_l[-n])*backward_l[-1][self.s2.name]*self.s2.to_s2 +
                                         self.s3.emis.get_emis(signal_l[-n])*backward_l[-1][self.s3.name]*self.s2.to_s3)
            scores_dict[self.s3.name] = (self.s1.emis.get_emis(signal_l[-n])*backward_l[-1][self.s1.name]*self.s3.to_s1 +
                                         self.s2.emis.get_emis(signal_l[-n])*backward_l[-1][self.s2.name]*self.s3.to_s2 +
                                         self.s3.emis.get_emis(signal_l[-n])*backward_l[-1][self.s3.name]*self.s3.to_s3)
            backward_l.append(scores_dict)
        backward_l.append(self.s1.emis.get_emis(signal_l[-n])*backward_l[-1][self.s1.name]*self.entrance.to_s1 +
                          self.s2.emis.get_emis(signal_l[-n])*backward_l[-1][self.s2.name]*self.entrance.to_s2 +
                          self.s3.emis.get_emis(signal_l[-n])*backward_l[-1][self.s3.name]*self.entrance.to_s3)
        return backward_l

    @staticmethod
    def _choose_state(forw_i, backw_i, full_prob):
        probs_dict = {'all_scores':{}}
        max_score = 0.0
        major_state_name = None
        state_names = forw_i.keys()
        for s_name in state_names:
            score = forw_i[s_name]*backw_i[s_name]/full_prob
            if score > max_score:
                major_state_name = None
                major_state_name = s_name
                max_score = score
            probs_dict['all_scores'][s_name] = score
        probs_dict['max_score'] = major_state_name
        return probs_dict


    @staticmethod
    def norm_trans(_s1, _s2, _s3):
        norm_c = 1.0/(float(_s1) + float(_s2) + float(_s3))
        return float(_s1)/norm_c, float(_s2)/norm_c, float(_s3)/norm_c


class State():
    def __init__(self, emis, _s1, _s2, _s3, _exit):
        self.emis = emis
        self.name = emis.name
        self.to_s1, self.to_s2, self.to_s3 = SignalHMM.norm_trans(_s1=_s1, _s2=_s2, _s3=_s3)
        self.to_exit = _exit


class Entrance:
    def __init__(self, _s1, _s2, _s3):
        self.to_s1 = _s1
        self.to_s2 = _s2
        self.to_s3 = _s3