import scipy as sc
import sys

def read_oneProp_DM(filename):
    """Reads oneProp density matrix in lj presentation as scipy array.
Arguments:
   filename -- path to the density matrix (stored in the OneProp format file.
Output:
   scipy array with density matrix. Basis functions sorted in the following way:
|S,0.5,0.5>  |S,0.5,-0.5>  |P,1.5,1.5>  |P,1.5,0.5>  |P,1.5,-0.5>  |P,1.5,-1.5>
|P,0.5,0.5>  |P,0.5,-0.5>  |D,2.5,2.5>  |D,2.5,1.5>  |D,2.5,0.5>  |D,2.5,-0.5>
|D,2.5,-1.5> |D,2.5,-2.5>  |D,1.5,1.5>  |D,1.5,0.5>  |D,1.5,-0.5>  |D,1.5,-1.5>,
etc.
    """
    res = sc.loadtxt(filename, skiprows=2, converters={0: lambda s: 0})

# simple checks for shape of the obtained array - expected
# number of rows is stored in the first line of the OneProp DM file
# number of columns must equal to the number of rows plus 1
    N, M = sc.shape(res)
    assert(N+1 == M)
    f = open(filename, 'r')
    Nt = int(f.readline())
    assert(N == Nt)
    return res[:, 1:]  # trim first column, where function symbols were placed.


def thr_q_num(N, nnum):
    """ yields basis functions integers from 0 to N-1 along with corresponding
quantum numbers values and basis function numbers as described
in the onePropFormat.
Usage:
    for ii, n, l, j, mj in thr_q_num_(N, [l0n, l1n, l2n]):
        ...

Arguments:
    N -- maximal index value
    nnum -- sequence of basis function numbers for given l values i.e.
    nnum[0] -- number of s basis functions, nnum[1] number of p basis functions,
    etc.
Result:
    iterator for using in loops, see description above

    """
    ii = 0
    for l in xrange(len(nnum)):  # l loop
        for n in xrange(nnum[l]):  # n loop
            j2 = 2*l + 1
            while j2 >= max(0, 2*l-1):  # j loop, max function for case l=0
                for mj2 in range(j2, -j2 - 1, -2):
                    yield (ii, n, l, j2*0.5, mj2*0.5)
                    ii += 1
                    if ii >= N:
                        raise StopIteration  # stop iteration
                j2 -= 2


def diagonal_masks(N, nnum):
    """ makes scipy masks, that give submatrices of density matrix,
corresponding to the each possible triple of the l, j, mj values.
Arguments:
    N -- size of the density matrix.
    nnum -- sequence of basis function numbers for given l values i.e.
    nnum[0] -- number of s basis functions, nnum[1] number of p basis functions,

Output:
    dictionary res[(l,j)] = mask, where l and j is the quantum numbers
    and mask is the mask"""

    res = {}
    for i, n, l, j, mj in thr_q_num(N, nnum):
        if (l, j, mj) not in res.keys():
            res[(l, j, mj)] = sc.zeros((N,), dtype=bool)
        for ii, nn, ll, jj, mmj in thr_q_num(N, nnum):
            if l != ll or j != jj or mj != mmj:
                continue
            res[(l, j, mj)][ii] = True
    return res


def test():
    # oneProp test on correct DM
    dm = read_oneProp_DM('./tests/JRed_DM_Re')
    # oneProp test on empty DM
    try:
        print read_oneProp_DM('./tests/empty_DM')
        raise ValueError("passed on empty!!!")
    except Exception:
        pass
    # oneProp test on non-square DM
    try:
        print read_oneProp_DM('./tests/rect')
        raise ValueError("passed on rectangular!!!")
    except AssertionError:
        pass
    try:
        print read_oneProp_DM('./tests/non_consist')
        raise ValueError("passed on consistent!!!")
    except AssertionError:
        pass
    for ii, n, l, j, mj in thr_q_num(36, [2, 2, 2]):
        print "{0} {1}{2}{3}/2, {4}".format(ii, n+1, 'spdfgh'[l], int(2*j), mj)
    # diagonal matrices tests ( one basis function per l value)
    for ljm, mask in diagonal_masks(18, [2, 2, 2, 2]).iteritems():
        print "-"*20
        print ljm
        print "-"*20
        print mask
        print dm[mask][:, mask]  # black magic here
    # diagonal matrices tests ( two basis function per l value)
    dm2 = read_oneProp_DM('JRed_DM_Re')
    for ljm, mask in diagonal_masks(36, [2, 2, 2, 2]).iteritems():
        print "-"*20
        print ljm
        print "-"*20
        print dm2[mask][:, mask]  # black magic here

def get_dm_blocks(filename, N, nnum):
    dm2 = read_oneProp_DM(filename)
    res = {}
    for ljm, mask in diagonal_masks(N, nnum).iteritems():
        l, j, m = ljm
        l = "spdfgh"[l]
        j = int(j*2)
        tmp =dm2[mask][:, mask]
        if (l, j) not in res.keys():
            res[l,j] = tmp*1e0
        else:
            res[l,j] += tmp
    return res
if __name__ == '__main__':
    #test()
    l, j, N, nnum = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), map(int, sys.argv[4:])
    print get_dm_blocks('JRed_Dm_Re', N, nnum)[(l, j)]
