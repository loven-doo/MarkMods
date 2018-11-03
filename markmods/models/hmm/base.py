from abc import abstractmethod

from markmods.models.base import ModelBase


class HMMBase(ModelBase):

    def fit(cls):
        pass

    def predict(cls):
        pass

    @abstractmethod
    def viterbi(cls):
        pass

    @abstractmethod
    def fb(cls):
        pass

    @abstractmethod
    def forward(cls):
        pass

    @abstractmethod
    def backward(cls):
        pass

    @abstractmethod
    def bw(cls):
        pass
