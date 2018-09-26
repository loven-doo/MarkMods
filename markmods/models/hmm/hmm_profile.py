def hmm_profile(mult_aln):
    profile = []
    aln_k = mult_aln.keys()
    seqs_n = len(aln_k)
    for n in range(0, len(mult_aln[aln_k[0]])):
        colomn = {}
        for k in aln_k:
            try:
                colomn[mult_aln[k][n]] = colomn[mult_aln[k][n]] + 1
            except KeyError:
                colomn[mult_aln[k][n]] = 1   
        profile.append(colomn)
    return profile
