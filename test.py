
def determine_REase(sequence, group, default_REase="BtsI", REases=REases):
    """determine the restriction enzyme"""

    overlap = 1
    for REase in REases:
        for i in range(1, len(group) + 1):
            if i == len(group):
                sub = sequence[group[i][0], group[i][1]]
            else:
                sub = sequence[group[i][0], group[i][1] + overlap]
            if 存在酶切位点:
                break
        enzyme = REase
        break

    return enzyme
    # BtsI_num = len(
    #     list(re.finditer(r'(?=GCAGTG|GTGACG|CACTGC|CGTCAC)',
    #                      str(oriSeq))))  # BtsI, oriSeq必须是string，pattern最好是raw string
    # BsrDI_num = len(list(re.finditer(r'(?=GCAATG|GTAACG|CATTGC|CGTTAC)', str(oriSeq))))  # BsrDI
    # print BtsI_num
    # print BsrDI_num
    #
    # for i in range(len(oligoGroups_starts) - 1):
    #     print "i", i
    #     # the first round
    #     if not BtsI_num and not BsrDI_num:  # no restriction sites
    #         print "no restriction sites"
    #     else:
    #         for site in re.finditer(r'(?=GCAGTG|GTGACG|CACTGC|CGTCAC)',
    #                                 str(oriSeq[oligoGroups_starts[i]:oligoGroups_starts[i + 1]])):
    #             BtsI_sites_starts.append(site.start())
    #         for site in re.finditer(r'(?=GCAATG|GTAACG|CATTGC|CGTTAC)',
    #                                 str(oriSeq[oligoGroups_starts[i]:oligoGroups_starts[i + 1]])):
    #             BsrDI_sites_starts.append(site.start())
    #         if BtsI_num >= BsrDI_num:
    #             restrictionSites_starts = BtsI_sites_starts
    #         else:
    #             restrictionSites_starts = BsrDI_sites_starts
    #         print "restrictionSites_starts", restrictionSites_starts
    #         print
    # return

























def determine_rease(dict1, sequence, rease_set):
    for rease in rease_set:
        for i in range(1, len(dict1) + 1):
            subSequence = sequence[dict1[i][0]: dict1[i][1] + 1]
            if rease.find_recognition_site(subSequence):
                break
    return rease
