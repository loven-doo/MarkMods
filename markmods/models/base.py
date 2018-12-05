from abc import abstractmethod, ABCMeta
from collections import defaultdict, OrderedDict

import pandas as pd
from scipy.stats import fisher_exact, chi2_contingency


class ModelBase(metaclass=ABCMeta):

    @abstractmethod
    def fit(self, X, Y):
        pass

    @abstractmethod
    def predict(self, X):
        pass

    def validate(self, X_test, Y_test, labels_to_remove=None):
        tables_2x2 = defaultdict(lambda: defaultdict(int))
        test_pred = self.predict(X=X_test)
        labels = tuple(test_pred[0].keys())
        confusion_dict = OrderedDict()
        for label in labels:
            confusion_dict[label] = [0]*len(labels)
        confusion_matrix = pd.DataFrame.from_dict(data=confusion_dict, orient='index')
        confusion_matrix.columns = labels
        for i in range(len(test_pred)):
            not_TN = [Y_test[i]]
            prob_labels = self._get_prob_labels(test_pred[i], n_prob_states=1, max_d=0.1)
            if Y_test[i] in prob_labels:
                tables_2x2[Y_test[i]]["TP"] += 1
                confusion_matrix[Y_test[i]][Y_test[i]] += 1
            else:
                tables_2x2[Y_test[i]]["FN"] += 1
                tables_2x2[prob_labels[0]]["FP"] += 1
                confusion_matrix[prob_labels[0]][Y_test[i]] += 1
                not_TN.append(prob_labels[0])
            for label in test_pred[i]:
                if label not in not_TN:
                    tables_2x2[label]["TN"] += 1
        if type(labels_to_remove) in (list, tuple, set, frozenset):
            for label in labels_to_remove:
                del tables_2x2[label]
        tables_2x2_array = list()
        for label in tables_2x2.keys():
            tables_2x2_array.append([[tables_2x2[label]["TP"], tables_2x2[label]["FP"]],
                                     [tables_2x2[label]["FN"], tables_2x2[label]["TN"]]])
            tables_2x2[label]["Fisher_P-value"] = fisher_exact(tables_2x2_array[-1],
                                                               alternative='greater')[1]
            try:
                precision = float(tables_2x2[label]["TP"]) / float(tables_2x2[label]["TP"]+tables_2x2[label]["FP"])
                recall = float(tables_2x2[label]["TP"]) / float(tables_2x2[label]["TP"]+tables_2x2[label]["FN"])
                tables_2x2[label]["F1_score"] = 2.0 * (precision*recall) / (precision+recall)
            except ZeroDivisionError:
                tables_2x2[label]["F1_score"] = None
        tables_2x2["Chi-sq_P-value"] = chi2_contingency(tables_2x2_array)[1]
        return tables_2x2, confusion_matrix

    @staticmethod
    def _get_prob_labels(labels_scores, n_prob_states=1, max_d=0.1):
        prob_labels = [(None, 0)]*n_prob_states
        for label in labels_scores:
            for i in range(n_prob_states):
                if i > 0:
                    if prob_labels[i-1][1]-labels_scores[label] > max_d:
                        break
                if prob_labels[i] is None:
                    prob_labels[i] = (label, labels_scores[label])
                    break
                elif labels_scores[label]-prob_labels[i][1] > 0:
                    prob_labels.insert(i, (label, labels_scores[label]))
                    break
        return list(filter(None, map(lambda l: l[0], prob_labels[:n_prob_states])))

    @abstractmethod
    def dump(self, scheme_path, **kwargs):
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


def aggregate_array(array, aggr_level, group_lims=None):
    if aggr_level < 1:
        return array
    if group_lims is None:
        group_lims = [list() for k in range(aggr_level)]
        i = -1
    else:
        try:
            i = group_lims[aggr_level - 1][-1]
        except IndexError:
            i = -1
    aggr_array = list()
    if aggr_level == 1:
        for array_elm in array:
            aggr_array.extend(array_elm)
            i += len(array_elm)
            group_lims[0].append(i)
    else:
        for array_elm in array:
            aggr_array.extend(aggregate_array(array=array_elm, aggr_level=aggr_level-1, group_lims=group_lims)[0])
            i += len(array_elm)
            group_lims[aggr_level-1].append(i)
    return aggr_array, group_lims


def group_array(aggr_array, group_lims):
    if len(group_lims) == 1:
        return _group_level(aggr_array, group_lims.pop(0))
    aggr_array = _group_level(aggr_array, group_lims.pop(0))
    return group_array(aggr_array, group_lims)


def _group_level(aggr_array, level_lims):
    grouped_array = [list()]
    for i in range(len(aggr_array)):
        grouped_array[-1].append(aggr_array[i])
        if i == level_lims[0]:
            level_lims.pop(0)
            grouped_array.append(list())
    return list(filter(None, grouped_array))

