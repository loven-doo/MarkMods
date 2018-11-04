from abc import abstractmethod, ABCMeta

from markmods.constants import ZERO_PROB
from markmods.models.base import ModelBase


class HMMBase(ModelBase, metaclass=ABCMeta):

    def __init__(self, scheme, zero_prob=ZERO_PROB):
        """
        HMM is recommended to have no less than 3 states
        The usage of zero probabilities is not recommended
        :param scheme: the scheme of the HMM in json format
            For Single Dimension HMM see examples in markmods.tests.test_models_hmm.TestSDHMM and
                markmods/tests/test_data/models_hmm/sdhmm_scheme.json
            For Multiple Dimension HMM see examples in markmods.tests.test_models_hmm.TestMDHMM and
                markmods/tests/test_data/models_hmm/mdhmm_scheme.json
        """
        self._scheme = scheme
        self.zero_prob = zero_prob
        self.states = list()
        self.nd = 1
        self.create_hmm()

    def fit(self):
        pass

    def predict(self):
        pass

    @abstractmethod
    def create_hmm(self):
        pass

    @property
    @abstractmethod
    def scheme(self):
        pass

    @abstractmethod
    def viterbi(self):
        pass

    @abstractmethod
    def fb(self):
        pass

    @abstractmethod
    def forward(self):
        pass

    @abstractmethod
    def backward(self):
        pass

    @abstractmethod
    def bw(self):
        pass


class State(object):

    def __init__(self, emis_func, emis_polinomial_c, name, group, self_trans, trans):
        self.name = name
        self.group = group
        self.emis_func = emis_func
        self.emis_polinimial_c = emis_polinomial_c
        self.self_trans = self_trans
        self.trans = trans

    def emis(self, x):
        efr = self.emis_func(x)
        e = 0.0
        for i in range(len(self.emis_polinimial_c)):
            e += self.emis_polinimial_c[i]*(efr**i)
        return e
