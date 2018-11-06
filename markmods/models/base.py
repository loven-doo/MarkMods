from abc import abstractmethod, ABCMeta


class ModelBase(metaclass=ABCMeta):

    @abstractmethod
    def fit(self):
        pass

    @abstractmethod
    def predict(self):
        pass

    @abstractmethod
    def dump(self, scheme_path, keep_groups):
        pass

    @abstractmethod
    @classmethod
    def load(cls, scheme_path):
        pass
