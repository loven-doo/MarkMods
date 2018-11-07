def combine_parts(avail_parts, without_repeats=False, sort=False):
    """
    The function to generate all combination from listed parts
    :param avail_parts: list of parts variants listed (list of lists)
    :param without_repeats: bool, if True combined parts of combined strings will be unique
    :param sort: "+" - combined parts will be in ascending order; "-" - descending order
    :return: list of combined strings
    """
    combined_list = list()
    if len(avail_parts) == 1:
        combined_list = [[part] for part in avail_parts[0]]
    elif len(avail_parts) == 2:
        for part2 in avail_parts[1]:
            for part1 in avail_parts[0]:
                combined_pair = [part1, part2]
                if without_repeats:
                    if len(set(combined_pair)) == 2:
                        if sort == "+":
                            combined_list.append(sorted(combined_pair))
                        elif sort == "-":
                            combined_list.append(sorted(combined_pair, reverse=True))
                        else:
                            combined_list.append(combined_pair)
                else:
                    combined_list.append(combined_pair)
    else:
        for part2 in combine_parts(avail_parts[1:], without_repeats=without_repeats, sort=sort):
            for part1 in avail_parts[0]:
                comb = [part1] + part2
                if without_repeats:
                    if len(set(comb)) == len(comb):
                        if sort == "+":
                            combined_list.append(sorted(comb))
                        elif sort == "-":
                            combined_list.append(sorted(comb, reverse=True))
                        else:
                            combined_list.append(comb)
                else:
                    combined_list.append(comb)
    return combined_list
