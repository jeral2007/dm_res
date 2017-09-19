import source_hfd_python.hfd_dat as hfd
import source_hfd_python.hfd_in_out as hinout
import scipy as sc
from scipy.linalg import inv, sqrtm


def eval_dm():
    import subprocess as sp
    sp.call(['hfj-dft.exe'])
    basis = hfd.Hfj('HFJ_basis.DAT')

class Evaluator(object):
    def __init__(self, rc, hfj):
        self.nc = sc.argwhere(hfj.grid < rc)[-1]
        self.wh = hfj.weights[:self.nc]*hfj.h
        self.hfj = hfj

    def orb_dot(self, a, b):
            return sc.trapz(a[:self.nc]*b[:self.nc]*self.wh)

    def get_dm_block(self, fis, bfis, qfracs):
        fvals = [self.hfj.getorb_by_number(fi)[0] for fi in fis]
        bfvals = [self.hfj.getorb_by_number(fi)[0] for fi in bfis]
        norms = [self.h
        smat = sc.array([[self.orb_dot(f2, f1)/self.orb_dot( for f1 in bfvals] for f2 in bfvals])
        smat
        smati =
        #smati = sqrtm(inv(smat))
        #smati = sc.eye(smat.shape[0],smat.shape[1])
        def aux(orb, basis):
            #  prepare matrix for orb
            cf = sc.dot(smati, [self.orb_dot(orb, b) for b in basis])
            return sc.outer(cf, cf)

        return sum(q*aux(f, bfvals) for q, f in zip(qfracs, fvals))

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

e = Evaluator(0.7, test)

orbs_from_hfj_res = hinout.make_orbitals_parser(fmt=['no', 'nl', 'j', 'occ',
                                                     'kp', 'en', 'd', 'delt',
                                                     'Rm2'])
fil = open('HFJ.RES', 'r')

orbs = orbs_from_hfj_res(fil)[0]
print orbs
inds = [i for i, o in enumerate(orbs) if 'S' in o['nl'] and o['j'] == '1/2']
fracs = [orbs[i]['occ'] for i in inds]
b_inds = [inds[1], inds[2]]

print inds, b_inds
coefs = []
ci = 0
for i in inds:
    orb = test.getorb_by_number(i)[0]
    coefs += [sc.array([e.orb_dot(orb, test.getorb_by_number(j)[0])
                        for j in b_inds])]

print coefs
res =  e.get_dm_block(inds, b_inds, [orbs[i]['occ'] for i in inds])
print res
print inv(res)
print orbs[0]['occ']
mt = sc.array([[3.617, 8.340], [8.340, 19.393]])
print sc.dot(inv(mt), res)

