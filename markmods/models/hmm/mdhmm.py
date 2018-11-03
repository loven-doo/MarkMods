from markmods.constants import ZERO_PROB
from markmods.models.hmm.base import HMMBase


class MDHMM(HMMBase):
    # Multiple Dimension Hidden Markov Model

    def __init__(self, scheme, zero_prob=ZERO_PROB):
        """
        HMM is recommended to have no less than 2 states
        The usage of zero probabilities is not recommended
        :param scheme: the scheme of the HMM in json format
            (see example in markmods.tests.test_models_hmm.TestMDHMM and
            markmods/tests/test_data/models_hmm/mdhmm_scheme.json)
        """
        pass

    def viterbi(self):
        pass

    def fb(self):
        pass

    def forward(self):
        pass

    def backward(self):
        pass

    def bw(self):
        pass


class State(object):

    def __init__(self, emis_func, name, group, trans):
        pass
