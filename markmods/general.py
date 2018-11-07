def combine_str(avail_parts, join_s="", without_repeats=False, sort=False):
    """
    The function to generate all combination from listed parts
    :param avail_parts: list of parts variants listed (list of lists)
    :param join_s: symbol to join parts
    :param without_repeats: bool, if True combined parts of combined strings will be unique
    :param sort: 
    :return: list of combined strings
    """
    combined_list = list()
    if len(avail_parts) == 1:
        combined_list = avail_parts[0]
    elif len(avail_parts) == 2:
        for part2 in avail_parts[1]:
            for part1 in avail_parts[0]:
                combined_list.append(part1+join_s+part2)
    else:
        for part2 in combine_str(avail_parts[1:], join_s):
            for part1 in avail_parts[0]:
                combined_list.append(part1+join_s+part2)
    return combined_list
