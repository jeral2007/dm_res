import source_hfd_python.hfd_dat as hfd
import source_hfd_python.hfd_in_out as hinout
import scipy as sc
from scipy.linalg import inv


def eval_dm():
    import subprocess as sp
    sp.call(['hfj-dft.exe'])


class Evaluator(object):
    def __init__(self, rc, hfj):
        self.nc = sc.argwhere(hfj.grid < rc)[-1]
        self.wh = hfj.weights[:self.nc]*hfj.h
        self.hfj = hfj

    def orb_dot(self, a, b):
            return sc.trapz(a[:self.nc]*b[:self.nc]*self.wh)

    def get_dm_block(self, fis, bfis, qfracs):
        fvals = [q*self.hfj.getorb_by_number(fi)[0] for q, fi in zip(qfracs,
                                                                     fis)]
        bfvals = [self.hfj.getorb_by_number(fi)[0] for fi in bfis]
        norms = [self.orb_dot(orb, orb) for orb in bfvals]
        smati = inv([[self.orb_dot(f2, f1) for f1 in bfvals] for f2 in bfvals])
        dots_vec = [sum(self.orb_dot(bf, f) for f in fvals) for bf in bfvals]
        coefs = sc.dot(smati, dots_vec)
        return sc.outer(coefs, coefs)
test = hfd.Hfj('HFJ.DAT')


r_and_orbs = sc.zeros((491, 13))

r_and_orbs[:, 0] = test.grid

for ni in xrange(1, 13):
    r_and_orbs[:, ni] = test.getorb_by_number(ni)[0]
    orb = r_and_orbs[:, ni]
    print(sc.trapz(orb**2*test.weights*test.h))
    print(sc.trapz(orb**2*test.grid*test.weights*test.h))


sc.savetxt('test.out', r_and_orbs)
sc.savetxt('test2.out', test.weights)

e = Evaluator(1.0, test)

orbs_from_hfj_res = hinout.make_orbitals_parser(fmt=['no', 'nl', 'j', 'occ',
                                                     'kp', 'en', 'd', 'delt',
                                                     'Rm2'])
fil = open('HFJ.RES', 'r')

orbs = orbs_from_hfj_res(fil)[0]
print orbs
inds = [i for i, o in enumerate(orbs) if 'P' in o['nl'] and o['j'] == '1/2']

b_inds = inds[1:]

ci = 0
smat = sc.zeros((2,2))
for i in b_inds:
    tmp = []
    cj = 0
    for j in b_inds:
        smat[ci, cj]= e.orb_dot(r_and_orbs[:,i],r_and_orbs[:,j])
        cj +=1
    ci += 1
print smat
print inv(smat)
print sc.dot(inv(smat), smat)
print e.get_dm_block(inds, b_inds, [orbs[i]['occ'] for i in inds])
