from markmods.models.hmm.base import HMMBase, State


class MDHMM(HMMBase):
    # Multiple Dimension Hidden Markov Model

    def prepare_trans(self, state_conf, state_name):
        state_trans = dict()
        for pdim in self._scheme["trans"]:
            self.prepare_pd_trans(state_conf=state_conf, pd_trans=self._scheme["trans"][pdim])

        self.nd = 2
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
