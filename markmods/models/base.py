from abc import abstractmethod, ABCMeta


class ModelBase(ABCMeta):

    @abstractmethod
    def fit(cls):
        pass

    @abstractmethod
    def predict(cls):
        pass
