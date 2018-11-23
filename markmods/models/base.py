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


def aggregate_array(array, aggr_level, aggr_lims=None):
    if aggr_level < 1:
        return array
    aggr_array = list()
    if aggr_level == 1:
        for array_elm in array:
            for elm in array_elm:
                if type(elm) not in (int, bool, float):###
                    if type(elm) is dict:
                        for feature in elm:
                            if type(elm[feature]) not in (int, bool, float):
                                bill_rec_logger.info(elm)
                    else:
                        bill_rec_logger.info(elm)###
            aggr_array.extend(array_elm)
    else:
        for array_elm in array:
            aggr_array.extend(aggregate_array(array=array_elm, aggr_level=aggr_level-1))
    return aggr_array


def group_array(aggr_array, aggr_lims):

    return

