###imports
import get_draft_matrix
import get_clusters
import get_stop_words
############
###configure
conf_f = "config"
############
###functions
def read_conf (f):
    par_d = {}
    read = open(f)
    for lines in read:
        line = lines.strip()
        if len(line) == 0:
            line = None
            continue
        if line[0] == "#":
            line = None
            continue
        line_list = line.split("##")
        par = line_list[0].strip()
        par_list = par.split("=")
        par_value = define_value(par_list[1].strip())
        par_d[par_list[0].strip()] = par_value
        par_value = None
        par = None
        line = None
    return par_d
###
def define_value (v):
    value = None
    if v[0] == "(":
        value = read_tuple(v)
        print value
    else: 
        try:
            value = float(v)
        except ValueError:
            value = str(v)
    return value
###
def read_tuple (v):
    tup_l = []
    elm_l = []
    tup = v.strip("()").split(",")
    quotes_opened = False
    start_quot = None
    for elms in tup:
        elm = elms.strip()
        if quotes_opened == False:
            if elm[0] == "'" or elm[0] == '"':
                start_quot = elm[0]
                quotes_opened = True
                if elm[len(elm)-1] == elm[0]:
                    quotes_opened = False
                    tup_l.append(elm.strip(elm[0]))
                else:
                    end_space = ""
                    if len(elms) > len(elm):
                        if elm[len(elm)-1] == elms[len(elms)-2]:
                            end_space = None
                            end_space = " "
                    elm_l.append(elm.strip(elm[0])+end_space) 
            else:
                try:
                    tup_l.append(float(elm))
                except ValueError:
                    elm = None
                    print "error while reading configure tuple"
                    continue
        else:
            quotes_opened = quotes_an(elm, start_quot)
            if quotes_opened == True:
                elm_l.append(elms)
            else:
                start_space = ""
                if len(elms) > len(elm):
                    if elm[0] == elms[1]:
                        start_space = None
                        start_space = " "
                start_quot = None                
                elm_l.append(start_space+elm[0:len(elm)-1])
                tup_l.append(",".join(elm_l))
                elm_l = []
        elm = None
    return tuple(tup_l)
###
def quotes_an (e, start_q):
    n = 0
    for s in e:
        if s == start_q:
            n = n + 1
    if n/2 != float(n)/2.0:
        return False
    else:
        return True
############
###MAIN
params = read_conf(conf_f)
get_draft_matrix.main(params)
get_clusters.main(params)
get_stop_words.main(params)