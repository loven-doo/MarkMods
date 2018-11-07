from markmods.general import combine_parts
from markmods.models.hmm.base import HMMBase, State


class MDHMM(HMMBase):
    # Multiple Dimension Hidden Markov Model

    def _prepare_trans(self, state_conf, state_name):
        state_trans = dict()
        max_dim = 1
        for pdim in self._scheme["trans"]:
            state_trans[pdim] = self._prepare_pd_trans(state_conf=state_conf, pd_trans=self._scheme["trans"][pdim])
            pdim_list = pdim.split("-")
            if len(pdim_list) > max_dim:
                max_dim = len(pdim_list)
        self.nd = max_dim
        state_trans = self._expand_dim_keys(state_trans)
        state_self_trans = self._expand_dim_keys(state_conf["self_trans"])
        for pdim in state_self_trans:
            state_trans[pdim][state_name] = state_self_trans[pdim]
        return state_trans

    def _expand_dim_keys(self, dim_dict):
        for dim_key in list(dim_dict.keys()):
            avail_dims = list()
            dims_list = dim_key.split("-")
            for dim in dims_list:
                if len(dim) < 2:
                    avail_dims.append(["d"+str(i+1) for i in range(self.nd)])
                else:
                    avail_dims.append([dim,])
            dim_keys = ["-".join(dkey_list) for dkey_list in combine_parts(avail_dims, without_repeats=True, sort="+")]
            for dkey in dim_keys:
                if dkey not in dim_dict:
                    dim_dict[dkey] = dim_dict[dim_key]
            if "-d-" in dim_key.lower() or dim_key.lower()[:2] == "d-" or dim_key.lower()[-2:] == "-d" or \
                    dim_key.lower() == "d":
                dim_dict.pop(dim_key)
        return dim_dict

    @property
    def scheme(self):
        return

    def viterbi(self, data_array):
        pass

    def fb(self, data_array):
        pass

    def forward(self, data_array):
        pass

    def backward(self, data_array):
        pass

    def bw(self):
        pass
