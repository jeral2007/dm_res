import source_hfd_python.hfd_dat as hfd
import source_hfd_python.hfd_in_out as hinout
import scipy as sc
from scipy.linalg import inv, sqrtm, eig
from scipy import sqrt
import sys
import readOneProp as rop
def eval_dm():
    import subprocess as sp
    sp.call(['hfj-dft.exe'])
    basis = hfd.Hfj('HFJ_basis.DAT')

class Evaluator(object):
    def __init__(self, rc, hfj, hfj_bas):
        self.nc = sc.argwhere(hfj.grid < rc)[-1]
        self.wh = hfj.weights[:self.nc]*hfj.h
        self.hfj = hfj
        self.hfj_bas = hfj_bas
    def orb_dot(self, a, b):
        return sc.trapz(a[:self.nc]*b[:self.nc]*self.wh)
    def norm(self, a):
        return sc.trapz(a[:self.nc]**2*self.wh)

    def get_dm_block(self, fis, bfis, qfracs):
        fvals = [self.hfj.getorb_by_number(fi)[0] for fi in fis]
        bfvals = [self.hfj_bas.getorb_by_number(fi)[0] for fi in bfis]
#        smat = sc.array([[self.orb_dot(f2, f1) for f1 in bfvals] for f2 in bfvals])
#        w, cbs = eig(smat)
#        bas2  = [sum(cbs[i,j]*f for i, f in enumerate(bfvals))
#                 for j in xrange(cbs.shape[0])]

        #smati = (inv(smat))
        def aux(orb, basis):
            #  prepare matrix for orb
            #cf = sc.dot(smati, [self.orb_dot(orb, b)/self.norm(b) for b in basis])
            cf = sc.array([self.orb_dot(orb, b)/self.norm(b) for b in basis])*sqrt(0.7)
            orb2 = orb*1e0
            for b, c in zip(basis, cf):
                orb2 -= b*c
#            print "||orb2||^2 = {} ".format(self.orb_dot(orb2, orb2))
 #           print "||orb||^2 = {} ".format(self.orb_dot(orb, orb))
            tmp_res =  sc.outer(cf, cf)
            #mult = sc.array([[1./sqrt(self.norm(b1)*self.norm(b2)) for b1 in basis] for b2 in basis])
            return tmp_res

        return (sum(q*aux(f, bfvals) for q, f in zip(qfracs, fvals)))

def ljs(ls):
    for l in ls:
        if l == 's':
            yield l, '1/2'
        else:
            lv = "spdfgh".index(l)
            yield l, str(2*lv-1)+'/2'
            yield l, str(2*lv+1)+'/2'

def tot_error(hfj_fname, bas_hfj_fname, ls):
    hfj = hfd.Hfj(hfj_fname+'.DAT')
    hfj_bas = hfd.Hfj(bas_hfj_fname+'.DAT')
    e = Evaluator(0.7, hfj, hfj_bas)
    orbs_from_hfj_res = hinout.make_orbitals_parser(fmt=['no', 'nl', 'j',
                                                         'occ', 'kp', 'en',
                                                         'd', 'delt', 'Rm2'])
    fil = open(hfj_fname+'.RES', 'r')
    orbs = orbs_from_hfj_res(fil)[0]
    fil.close()
    bfil = open(hfj_fname +'.res', 'r')
    borbs = orbs_from_hfj_res(bfil)[0]
    bfil.close()
    err = 0e0
    for larg, jarg in ljs(ls):
        inds = [i for i, o in enumerate(orbs) if larg.upper() in o['nl'] and
                o['j'] == jarg]
        b_inds = [i for i, o in enumerate(borbs) if larg.upper() in o['nl'] and
                o['j'] == jarg]
        b_inds = b_inds[-2:]
        tm1 =  e.get_dm_block(inds, b_inds, [orbs[i]['occ'] for i in inds])
        tm = op_dm_blocks[larg.lower(), int(jarg[0])]
        err += sc.sum((tm1/tm - sc.ones_like(res))**2)
    return err

if __name__ == '__main__':
    test = hfd.Hfj('HFJ_basis.DAT')

    e = Evaluator(0.7, test, test)

    orbs_from_hfj_res = hinout.make_orbitals_parser(fmt=['no', 'nl', 'j', 'occ',
                                                        'kp', 'en', 'd', 'delt',
                                                        'Rm2'])
    fil = open('HFJ_basis.RES', 'r')

    orbs = orbs_from_hfj_res(fil)[0]
    larg, jarg, ls = sys.argv[1], sys.argv[2], sys.argv[3]
    op_dm_blocks = rop.get_dm_blocks('JRed_DM_Re', 36, [2,2,2])

    inds = [i for i, o in enumerate(orbs) if larg.upper() in o['nl'] and o['j'] == jarg]
    print larg.upper(), jarg, inds
    b_inds = inds[-2:]
    res =  e.get_dm_block(inds, b_inds, [orbs[i]['occ'] for i in inds])
    tm = op_dm_blocks[larg.lower(), int(jarg[0])]
    print res
    print tm
    print res/tm - sc.ones_like(res)
    print orbs[0]['occ']
    fil.close()
    print tot_error('HFJ_basis', 'HFJ_basis', ls)
#mt = sc.array([[3.617 + 3.565, 8.340+8.201], [8.340+8.201, 19.393+19.021]])
#print  mt/res


