import re
from collections import defaultdict
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
        self.groups = defaultdict(set)
        self.states = dict()
        self.nd = 1
        self.create_hmm()

    def fit(self):
        pass

    def predict(self):
        pass

    def dump(self, scheme_path, keep_groups=False):
        pass

    @classmethod
    def load(cls, scheme_path):
        pass

    def create_hmm(self):
        for state_name in self._scheme["states"]:
            self.groups[self._scheme["states"][state_name]["group"]].add(state_name)

        for state_name in self._scheme["states"]:
            state_conf = self._scheme["states"][state_name]
            state_trans = self.prepare_trans(state_conf, state_name)
            self.states[state_name] = State(emis_func=state_conf["emis_func"],
                                            emis_polinomial_c=state_conf["emis_polynomial_c"],
                                            name=state_name,
                                            group=state_conf["group"],
                                            self_trans=state_conf["self_trans"],
                                            trans=state_trans,
                                            nd=self.nd)

    def prepare_pd_trans(self, state_conf, pd_trans):
        state_pd_trans = dict()
        for trans_i in pd_trans:
            group_to = None
            if state_conf["group"] in trans_i:
                trans_i_ori = trans_i.split(">")
                if len(trans_i_ori) == 2:
                    if trans_i_ori[0] == state_conf["group"]:
                        group_to = trans_i_ori[1]
                else:
                    trans_i_unori = trans_i.split("-")
                    if trans_i_unori[0] == state_conf["group"]:
                        group_to = trans_i_unori[1]
                    else:
                        group_to = trans_i_unori[0]
                for state_name in self.groups[group_to]:
                    if state_name in state_pd_trans:
                        if state_pd_trans[state_name] < pd_trans[trans_i]:
                            state_pd_trans[state_name] = pd_trans[trans_i]
                    else:
                        state_pd_trans[state_name] = pd_trans[trans_i]
        return state_pd_trans

    @abstractmethod
    def prepare_trans(self, state_conf, state_name):
        pass

    @abstractmethod
    @property
    def scheme(self):
        pass

    @abstractmethod
    def viterbi(self, data_array):
        pass

    @abstractmethod
    def fb(self, data_array):
        pass

    @abstractmethod
    def forward(self, data_array):
        pass

    @abstractmethod
    def backward(self, data_array):
        pass

    @abstractmethod
    def bw(self):
        pass


class State(object):

    def __init__(self, emis_func, emis_polinomial_c, name, group, self_trans, trans, nd=1):
        self.name = name
        self.group = group
        self.emis_func = emis_func
        self.emis_polinimial_c = emis_polinomial_c
        self.self_trans = self_trans
        self.trans = trans
        self.nd = nd

    def emis(self, x):
        efr = self.emis_func(x)
        e = 0.0
        for i in range(len(self.emis_polinimial_c)):
            e += self.emis_polinimial_c[i]*(efr**i)
        return e
