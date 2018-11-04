from abc import abstractmethod, ABCMeta


class ModelBase(metaclass=ABCMeta):

    @abstractmethod
    def fit(self):
        pass

    @abstractmethod
    def predict(self):
        pass
