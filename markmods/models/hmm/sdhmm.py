from markmods.models.hmm.base import HMMBase, State


class SDHMM(HMMBase):
    # Single Dimension Hidden Markov Model

    def _prepare_trans(self, state_conf, state_name):
        state_trans = self._prepare_pd_trans(state_conf=state_conf, pd_trans=self._scheme["trans"])
        state_trans[state_name] = state_conf["self_trans"]
        return state_trans

    @property
    def scheme(self):
        return

    def viterbi(self, data_array):
        pass

    def fb(self, data_array):
        pass

    def forward(self, data_array):
        pass

    def backward(self, data_array):
        pass

    def bw(self):
        pass
