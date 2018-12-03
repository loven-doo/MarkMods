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

    @classmethod
    @abstractmethod
    def load(cls, scheme_path):
        pass

    @property
    @abstractmethod
    def features(self):
        pass

    @property
    @abstractmethod
    def labels(self):
        pass


def aggregate_array(array, aggr_level, aggr_lims=None):
    if aggr_level < 1:
        return array
    aggr_array = list()
    if aggr_level == 1:
        for array_elm in array:
            aggr_array.extend(array_elm)
    else:
        for array_elm in array:
            aggr_array.extend(aggregate_array(array=array_elm, aggr_level=aggr_level-1))
    return aggr_array


def group_array(aggr_array, aggr_lims):

    return

