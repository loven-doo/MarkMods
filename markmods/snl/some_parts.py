class C:

        def get_biggest_freq(self, lines_list):
        biggest_freq_list = [0]
        biggest_freq = 0.0
        minute = 0.0# in seconds units
        msms = 0.0
        ms_level = 0
        rt0 = 0.0
        drt = 1.0
        for lines in lines_list:
            line = lines.strip()
            if line[0:13].lower() == "retentiontime":
                if ms_level == 2:
                    continue
                rt_str = line.split("=")[1].strip('"')
                try:
                    rt = float(rt_str)
                except ValueError:
                    rt = float(rt_str.strip("PTS"))
                if minute == 0.0:
                    minute = float((int(rt) / 60 + 1) + 60)
                    rt0 = rt - 1.0
                drt = rt - rt0
                rt0 = rt
                if rt > minute:
                    biggest_freq_list.append(biggest_freq)
                    biggest_freq = 0.0
                    minute += 60.0
            elif line.lower() == 'mslevel="2"':
                ms_level = 2
                msms += 1.0
            else:
                ms_level = 1
                if msms/drt > biggest_freq:
                    biggest_freq = msms/drt
                msms = 0.0

        signal_hmm = self._build_bf_hmm(max(biggest_freq_list))
        states = signal_hmm.prepare_signal(biggest_freq_list)# states=[{'max_score': <state_name>, 'all_scores':{<state_name>:score, ...}}, ...]
        corrected_bfl = self._correct_bf_drops(biggest_freq_list, states)
        self.biggest_freq_corr = pd.Series(corrected_bfl)
        self.biggest_freq_corr.name = 'max_msms_freq_any_minute'
        return corrected_bfl, self.plot_list(corrected_bfl, "Maximal MS/MS frequency any minute",
                                            "Time, minutes", "MS/MS frequency, Hz", "max_freq_plot.svg")

    @staticmethod
    def _build_bf_hmm(max_freq, delta=0.3, trans_high=0.85, trans_mid=0.5, trans_low=0.15):
        state1_emis = BFEmis(i_max=max_freq, d=1.0, name='high_freq', power=1.0)
        state2_emis = BFEmis(i_max=max_freq, d=delta, name='low_freq', power=1.0)
        state3_emis = BFEmis(i_max=max_freq, d=0.0, name='missed_scans', power=1.0)
        return SignalHMM(s1_e=state1_emis, s2_e=state2_emis, s3_e=state3_emis,
                         s1_s1=trans_high, s2_s2=trans_high, s3_s3=trans_low,
                         s1_s2=trans_mid, s1_s3=trans_high, s2_s3=trans_mid,
                         entr_s1=trans_low, entr_s2=trans_high, entr_s3=0.0,
                         s1_exit=trans_low, s2_exit=trans_high, s3_exit=0.0)


class BFEmis():
    def __init__ (self, i_max, d, name, power=1.0):
        self.emis = None
        self.i_max = float(i_max)
        self.delta = float(d)
        self.name = name
        self.power = power

    def get_emis(self, i, **kwargs):
        return 1.0 - (float(i)/self.i_max - self.delta)**(2.0*self.power)

