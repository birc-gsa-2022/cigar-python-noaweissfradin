"""A module for translating between alignments and edits sequences."""


from re import M


def get_edits(p: str, q: str) -> tuple[str, str, str]:
    """Extract the edit operations from a pairwise alignment.

    Args:
        p (str): The first row in the pairwise alignment.
        q (str): The second row in the pairwise alignment.

    Returns:
        str: The list of edit operations as a string.

    >>> get_edits('ACCACAGT-CATA', 'A-CAGAGTACAAA')
    ('ACCACAGTCATA', 'ACAGAGTACAAA', 'MDMMMMMMIMMMM')

    """
    assert len(p) == len(q)
    pp = p
    qq = q
    M = ''
    for i in range (len(p)):
        if(p[i]== '-'):
            pp = pp[:i] + pp[i+1:]
            M += 'I'
        elif(q[i]== '-'):
            qq = qq[:i] + qq[i+1:]
            M += 'D'
        else:
            M += 'M'  

    return pp, qq, M


def local_align(p: str, x: str, i: int, edits: str) -> tuple[str, str]:
    """Align two sequences from a sequence of edits.

    Args:
        p (str): The read string we have mapped against x
        x (str): The longer string we have mapped against
        i (int): The location where we have an approximative match
        edits (str): The list of edits to apply, given as a string

    Returns:
        tuple[str, str]: The two rows in the pairwise alignment

    >>> local_align("ACCACAGTCATA", "GTACAGAGTACAAA", 2, "MDMMMMMMIMMMM")
    ('ACCACAGT-CATA', 'A-CAGAGTACAAA')

    """
    
    return align(p,x[i:],edits)


def align(p: str, q: str, edits: str) -> tuple[str, str]:
    """Align two sequences from a sequence of edits.

    Args:
        p (str): The first sequence to align.
        q (str): The second sequence to align
        edits (str): The list of edits to apply, given as a string

    Returns:
        tuple[str, str]: The two rows in the pairwise alignment

    >>> align("ACCACAGTCATA", "ACAGAGTACAAA", "MDMMMMMMIMMMM")
    ('ACCACAGT-CATA', 'A-CAGAGTACAAA')

    """
    pp = p
    qq = q

    for i in range (len(edits)):
        if(edits[i]== 'I'):
            pp = pp[:i] + '-' + pp[i:]
        elif(edits[i]== 'D'):
            qq = qq[:i] + '-' + qq[i:]

    return pp, qq


def edit_dist(p: str, x: str, i: int, edits: str) -> int:
    """Get the distance between p and the string that starts at x[i:]
    using the edits.

    Args:
        p (str): The read string we have mapped against x
        x (str): The longer string we have mapped against
        i (int): The location where we have an approximative match
        edits (str): The list of edits to apply, given as a string

    Returns:
        int: The distance from p to x[i:?] described by edits

    >>> edit_dist("accaaagta", "cgacaaatgtcca", 2, "MDMMIMMMMIIM")
    5
    """
    local = local_align(p, x, i, edits)
    pp = local[0]
    qq = local[1]
    score = 0
    for i in len(pp):
        if (pp[i]!=qq[i]):
            score +=1
    return score
