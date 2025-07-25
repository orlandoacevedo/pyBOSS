from .utils import (
    zmat2cor,  find_all_rings_dfs, ClosedForm, align_onto_z_axis
)
from .constants import (
    ATOMIC_RADIUS, AMBER_TO_OPLSAA_SYNONYM_ATOM_TYPE, BONDED_ATOMIC_RADIUS,
    ATOMIC_ELECTRONEGATIVITY, ATOMIC_WEIGHTS, AMBER_SURFACE_AREA
)

import math


class InitSolutes:
    """Initialize Solute perturbation parameters
    
    Methods Execution Sequence after `cls.__init__`:

                get_xyzs                get_pars
                   |                       |
                    \                     /
                     \                   /   calc_nonbond_pars
                      \                 /   /   \
                       \               /  /     calc_sasa
                        \             / /
                         calc_neighbors   ---  calc_atoms_symmetry_list
                        /         |    \         \
                       /          |     \         symmetrize_charges  --- calc_coulombic_pars
        calc_pert_dhiedhral_pars  |      \
                   calc_pert_bond_pars    calc_pert_angle_pars
                                  \      /
                                   calc_pert_xyzs
                                  /      |       \
                calc_maximum_overlap     |      recentroid
                                  calc_onto_z_axis
    """
    def __init__(
        self,rc0=None,rc1=None,rc2=None,scllj=None,
        nonebn=None,tk=None,
        izlong=None,maxovl=None,
        ncent1=None,ncent2=None,
        qmscale=None,
        bond_pars=None,dihedral_pars=None,                  # file operation
        oplsaa_bond_pars=None, oplsaa_angle_pars=None,      # file operation
        solutezmat=None,                                    # file operation
        *args, **kws
    ):
        self.nice = True
        self.info = ''
        self.esq = 332.06
        self.dsqesq = pow(self.esq,0.5)
        self.rc0 = 0.0 if rc0 is None else rc0
        self.rc1 = 0.5 if rc1 is None else rc1
        self.rc2 = 1.0 if rc2 is None else rc2
        self.scllj = scllj if scllj else 1.0
        self.nonebn = True if nonebn else False
        self.tk = tk if tk else 298.15
        self.ncent1 = ncent1
        self.ncent2 = ncent2
        self.izlong = True if izlong else False
        self.maxovl = True if maxovl else False
        self.qmscale = 1.0 if qmscale is None else qmscale

        self.solutezmat = solutezmat
        self.bond_pars = bond_pars
        self.dihedral_pars = dihedral_pars
        self.oplsaa_bond_pars = oplsaa_bond_pars
        self.oplsaa_angle_pars = oplsaa_angle_pars

        if self.ncent1 and isinstance(self.ncent1,int):
            self._bo_recenter = True
            n = len(self.solutezmat['data'])
            if self.ncent1 <= 0 or self.ncent1 > n:
                self.info = f'Fatal: not valid recenter atom index: range: [1,{n}]: ncent1: {self.ncent1}'
                self.nice = False
                return
            if self.ncent2:
                if isinstance(self.ncent2,int) and self.ncent2 > 0 or self.ncent2 <= n:
                    pass
                else:
                    self.info = f'Fatal: not valid recenter atom index: range: [1,{n}]: ncent2: {self.ncent2}'
                    self.nice = False
                    return
            else:
                self.ncent2 = self.ncent1
        else:
            self._bo_recenter = False
        self.total_number_atoms = len(self.solutezmat['data'])

        self.solutesdata = {'reference':{}, 'first':{}, 'second':{}}

        self._full_bond_pars = [*self.solutezmat['bond_pars'], *self.bond_pars]
        self.amber_to_oplsaa_atom_type = {
            k.lower():v for k,v in AMBER_TO_OPLSAA_SYNONYM_ATOM_TYPE.items()
        }
        self.bonded_atomic_radius = {
            k.lower():(k,v) for k,v in BONDED_ATOMIC_RADIUS.items()
        }
        self.amber_surface_area = {
            k.lower():v for k,v in AMBER_SURFACE_AREA.items()
        }
        # number of atom offset
        self.atoms_number_offset = []
        b = 0
        for i in range(self.solutezmat['number_of_entries']):
            e = self.solutezmat['number_of_atoms_for_entry_'+str(i+1)] + b
            self.atoms_number_offset.append((b,e))
            b += e
        self.solutesdata['atoms_number_offset'] = self.atoms_number_offset

    def get_xyzs(self):
        if not self.nice: return
        mat = [i[3:9] for i in self.solutezmat['data']]
        self.solutesdata['xyz_initial'] = zmat2cor(mat)

        # deep copy
        final = [[j for j in i] for i in mat]
        refer = [[j for j in i] for i in mat]
        for v in self.solutezmat['bonds_geometry']:
            if v[1] == 1:           # for bond
                final[v[0]][1] = v[2]
                refer[v[0]][1] += self.rc0 * (v[2]-refer[v[0]][1])
            elif v[1] == 2:         # for angle
                final[v[0]][3] = v[2]
                refer[v[0]][3] += self.rc0 * (v[2]-refer[v[0]][3])
            else:                   # for dihedral
                final[v[0]][5] = v[2]
                refer[v[0]][5] += self.rc0 * (v[2]-refer[v[0]][5])
        self.solutesdata['reference'] = {'xyz': zmat2cor(refer) }
        self.solutesdata['xyz_final'] = zmat2cor(final)

    def get_pars(self):
        if not self.nice: return
        self.solutesdata['initial_bond_pars'] = []
        self.solutesdata['final_bond_pars'] = []
        self.solutesdata['reference']['atomtypes'] = []
        self.solutesdata['first']['atomtypes'] = []
        self.solutesdata['second']['atomtypes'] = []
        for v in self.solutezmat['data']:
            if v[1] == -1:
                pi = [0, 0, ' ', 0.0, 0.0, 0.0]    # for dummy atom
            else:
                # custom ff will always be in first
                bo = False
                for pi in self._full_bond_pars:
                    if v[1] == pi[0]:
                        bo = True
                        break
                if not bo:
                    self.info = f'Fatal: parameter not found: atom: {v[:8]}'
                    self.nice = False
                    return
            self.solutesdata['initial_bond_pars'].append(pi)
            if v[1] == -1 or v[2] == 0 or v[2] == v[1]:
                pf = [t for t in pi]  # deep copy
            else:
                bo = False
                for pf in self._full_bond_pars:
                    if v[2] == pf[0]:
                        bo = True
                        break
                if not bo:
                    self.info = f'Fatal: parameter not found: atom index: {v[:8]}'
                    self.nice = False
                    return
            self.solutesdata['final_bond_pars'].append(pf)
            
            # now global determination
            v1 = pi[1] if self.rc0 <= 0.5 else pf[1]
            v2 = pi[1] if self.rc1 <= 0.5 else pf[1]
            v3 = pi[1] if self.rc2 <= 0.5 else pf[1]
            self.solutesdata['reference']['atomtypes'].append(v1)
            self.solutesdata['first']['atomtypes'].append(v2)
            self.solutesdata['second']['atomtypes'].append(v3)
        # for SASA
        lt = []
        for i,p in enumerate(self.solutesdata['initial_bond_pars']):
            if p[1] > 0:
                if self.rc0 > 0.5:
                    p = self.solutesdata['final_bond_pars'][i]
                v = self.amber_surface_area[p[2].lower()]
            else:
                v = 1
            lt.append(v)
        self.solutesdata['amber_sasa_atomtypes'] = lt

    def calc_nonbond_pars(self):
        if not self.nice: return
        sqrtesq = pow(self.esq,0.5)
        rc = [self.rc0, self.rc1, self.rc2]
        for i,k in enumerate(['reference', 'first', 'second']):
            self.solutesdata[k]['Q'] = []
            self.solutesdata[k]['A'] = []
            self.solutesdata[k]['B'] = []
            self.solutesdata[k]['pert_sigma'] = []
            self.solutesdata[k]['pert_epsilon'] = []
            for j,pi in enumerate(self.solutesdata['initial_bond_pars']):
                pf = self.solutesdata['final_bond_pars'][j]
                sig = (pi[4]+rc[i]*(pf[4]-pi[4])) * self.scllj
                eps = (pf[5]+rc[i]*(pf[5]-pi[5])) * self.scllj
                self.solutesdata[k]['Q'].append(pi[3]*sqrtesq)
                self.solutesdata[k]['A'].append(math.sqrt(4.0*eps*sig**12))
                self.solutesdata[k]['B'].append(math.sqrt(4.0*eps*sig**6))
                self.solutesdata[k]['pert_sigma'].append(sig)
                self.solutesdata[k]['pert_epsilon'].append(eps)

        # look for H on N, O, F, P, CL, Br, and I, then set LJ parameters to zero
        for i in range(len(self.solutezmat['data'])):
            if self.solutezmat['data'][i][1] != 1: continue
            for j in self.nbor[i]:
                if self.solutezmat['data'][j][1] in [7,8,9,15,16,17,35,53]:
                    self.solutesdata['reference']['A'][i] = 0.0
                    self.solutesdata['reference']['B'][i] = 0.0
                    if not self.nonebn:
                        self.solutesdata['first']['A'][i] = 0.0
                        self.solutesdata['first']['B'][i] = 0.0
                        self.solutesdata['second']['A'][i] = 0.0
                        self.solutesdata['second']['B'][i] = 0.0
                    break

    def calc_sasa(self):
        if not self.nice: return
        self.solutesdata['radius_sa'] = [
            i*0.561231024 for i in self.solutesdata['reference']['pert_sigma']
        ]

    def calc_neighbors(self):
        if not self.nice: return
        # calculate neighbor list
        n = len(self.solutezmat['data'])
        nbor = {i:set() for i in range(n)}
        unique_cob = 0
        for a in range(n):
            if self.solutezmat['data'][a][1] == -1 or self.solutezmat['data'][a][1] == 100 \
                or self.solutezmat['data'][a][2] == 100:
                continue
            ia = self.solutesdata['xyz_initial'][a]
            fa = self.solutesdata['xyz_final'][a]
            ita = self.solutesdata['reference']['atomtypes'][a]
            fta = self.solutesdata['first']['atomtypes'][a]
            ira = ATOMIC_RADIUS[ita]
            fra = ATOMIC_RADIUS[fta]
            for b in range(n):
                if a == b: continue
                if self.solutezmat['data'][b][1] == -1 or self.solutezmat['data'][b][1] == 100 \
                    or self.solutezmat['data'][b][2] == 100:
                    continue
                bo = False
                for p in [*self.solutezmat['bonds_equal'],*self.solutezmat['bonds_additional']]:
                    if (a in p) or (b in p):
                        unique_cob += 1
                        nbor[a].add(b)
                        nbor[b].add(a)
                        bo = True
                        break
                if bo: continue
                if self.nonebn: continue
                itb = self.solutesdata['reference']['atomtypes'][b]
                ftb = self.solutesdata['first']['atomtypes'][b]
                # exclude H-H and Cl-Cl pairs
                if ita+fta+itb+ftb == 4: continue
                if ita == 17 and fta == 17 and itb == 17 and ftb == 17: continue
                ib = self.solutesdata['xyz_initial'][b]
                fb = self.solutesdata['xyz_final'][b]
                di = sum([u*u for u in [ia[s]-ib[s] for s in range(3)]])
                df = sum([u*u for u in [fa[s]-fb[s] for s in range(3)]])
                if di > 4.0 and df > 4.0: continue
                irb = ATOMIC_RADIUS[itb]
                frb = ATOMIC_RADIUS[ftb]
                qi = (ira+irb+0.1)**2
                qf = (fra+frb+0.1)**2
                if di <= qi or df <= qf:
                    nbor[a].add(b)
                    nbor[b].add(a)
                    if b > a:
                        unique_cob += 1
        self.nbor = nbor

    def calc_pert_dihedral_pars(self):
        if not self.nice: return
        #TODO for `dihedrals_"other"`
        dpars = []
        limp = []
        for i in self.solutezmat['dihedrals_variable']:
            p1 = self.solutezmat['data'][i]
            p2 = self.solutezmat['data'][p1[3]-1]
            p3 = self.solutezmat['data'][p1[5]-1]
            p4 = self.solutezmat['data'][p1[7]-1]

            if p1[1] < 0 or p2[1] < 0 or p3[1] < 0 or p4[1] < 0:
                dpars.append(([],'Ignored'))
                continue
            
            if p1[0] == 'DUM' or p2[0] == 'DUM' or p3[0] == 'DUM' or p4[0] == 'DUM':
                dpars.append(([],'Ignored'))
                continue

            # check whether they are in the same solute
            for t in self.atoms_number_offset:
                if i >= t[0] and i < t[1]:
                    break
            if p1[7]-1 < t[0] or p1[7]-1 >= t[1]:       # `data` starts from 1
                print(f'Note: dihedral not in the same solute: ignored: {p1}')
                dpars.append(([],'Ignored'))
                continue
            
            # for improper dihedrals
            if p1[5]-1 in self.nbor[i]:              # 1-3 improper
                asy = '{:<2}'.format(self.solutesdata['initial_bond_pars'][p1[5]-1][2])
                nnb = len(self.nbor[i])
            elif p1[7]-1 in self.nbor[p1[3]-1]:      # 2-4 improper
                asy = '{:<2}'.format(self.solutesdata['initial_bond_pars'][p1[7]-1][2])
                nnb = len(self.nbor[p1[3]-1])
            else:
                asy = '  '
            if asy == '  ':
                limp.append(False)
                # construct dihedrals   --   reversed order, but it does not matter
                da = [
                    self.solutesdata['initial_bond_pars'][p1[7]-1][2],
                    self.solutesdata['initial_bond_pars'][p1[5]-1][2],
                    self.solutesdata['initial_bond_pars'][p1[3]-1][2],
                    self.solutesdata['initial_bond_pars'][i][2],
                ]
                dstr = '{:<2}-{:<2}-{:<2}-{:<2}'.format(*da)
                da = dstr.split('-')        # reproduce
                knq = 0
                for s in da:
                    if s[0] == '?': knq += 1
                    if s[1] == '?': knq += 1

                kmax = 0
                pf = None
                for p in self.dihedral_pars:
                    knti = kntf = 0
                    for t in range(2):
                        e = 2 - t
                        if da[0][t] == p[9][t]: knti += e
                        if da[1][t] == p[10][t]: knti += e + 1
                        if da[2][t] == p[11][t]: knti += e + 1
                        if da[3][t] == p[12][t]: knti += e

                        if da[3][t] == p[9][t]: kntf += e
                        if da[2][t] == p[10][t]: kntf += e + 1
                        if da[1][t] == p[11][t]: kntf += e + 1
                        if da[0][t] == p[12][t]: kntf += e
                    k = 2 * max(knti,kntf) + knq
                    if k == 32:     # means exact match
                        dpars.append((p, 'Exactly match'))
                        kmax = 0    # reset
                        break
                    if k > kmax:
                        kmax = k
                        pf = p      # alias
                if kmax > 0:        # means guessing
                    dpars.append((pf, 'Guess'))
            else:
                limp.append(True)
                index = None
                if asy[0] == 'C':
                    if asy[1] == ' ':
                        index = 160
                    elif nnb == 3:
                        index = 162
                elif asy[0] == 'N':
                    if asy[1] in [' ', 'A', '2', 'O']:
                        index = 161
                if index:
                    bo = False
                    for p in self.dihedral_pars:
                        if p[0] == index:
                            dpars.append((p, 'Improper'))
                            bo = True
                            break
                    if not bo:
                        print(f'Fatal: dihedral parameter not found: index: {index}')
                else:
                    dpars.append(([],'Ignored'))

        # works only on initial and final dihedrals are the same
        # except `dihedrals_variable`, not other types are allowed
        self.solutesdata['initial_dihedral_pars'] = dpars
        self.solutesdata['final_dihedral_pars'] = dpars

        #TODO restrained dihedral `index=500`

        lrv = []
        lrp = []
        lfv = []
        lfp = []
        lsv = []
        lsp = []
        for gi,gf in zip(
            self.solutesdata['initial_dihedral_pars'],
            self.solutesdata['final_dihedral_pars']
        ):
            gi = gi[0]
            gf = gf[0]
            dv = [gf[t]-gi[t] for t in range(1,5)]
            dp = [gf[t]-gi[t] for t in range(5,8)]
            lrv.append([gi[t]+u*self.rc0 for u,t in zip(dv,range(1,5))])
            lrp.append([gf[t]+u*self.rc0 for u,t in zip(dp,range(5,8))])
            lfv.append([gi[t]+u*self.rc1 for u,t in zip(dv,range(1,5))])
            lfp.append([gf[t]+u*self.rc1 for u,t in zip(dp,range(5,8))])
            lsv.append([gi[t]+u*self.rc2 for u,t in zip(dv,range(1,5))])
            lsp.append([gf[t]+u*self.rc2 for u,t in zip(dp,range(5,8))])
        self.solutesdata['reference']['pert_dihedrals_v'] = lrv
        self.solutesdata['reference']['pert_dihedrals_p'] = lrp
        self.solutesdata['first']['pert_dihedrals_v'] = lfv
        self.solutesdata['first']['pert_dihedrals_p'] = lfp
        self.solutesdata['second']['pert_dihedrals_v'] = lsv
        self.solutesdata['second']['pert_dihedrals_p'] = lsp

        ldt = []
        rings = find_all_rings_dfs(self.nbor)
        for bo,i,g in zip(
            limp,
            self.solutezmat['dihedrals_variable'],
            self.solutesdata['reference']['pert_dihedrals_v']
        ):
            if bo:
                ldt.append(5.0)
                continue
            p = self.solutezmat['data'][i]
            d1 = p[3] - 1       # index starts from 1
            d2 = p[5] - 1
            n2 = len([1 for t in rings if d1 in t])     # be aware of the sequence
            n1 = len([1 for t in rings if d2 in t])     # dihedrals are got from right to left
            v = 15.0
            if n1 != 0: v = 10.0
            if n2 != 0: v = 5.0
            if g[1] > 4.0:
                v = 2.0
            elif g[1] > 3.0:
                v = 5.0
            ldt.append(v)
        self.solutesdata['pert_delta_dihedrals'] = ldt

    def calc_pert_bond_pars(self):
        if not self.nice: return
        dpars = []
        for i in self.solutezmat['bonds_variable']:
            p1 = self.solutezmat['data'][i]
            p2 = self.solutezmat['data'][p1[3]-1]

            if p1[1] < 0 or p2[1] < 0:
                dpars.append(([],'Ignored'))
                continue
            
            if p1[0] == 'DUM' or p2[0] == 'DUM':
                dpars.append(([],'Ignored'))
                continue

            # check whether they are in the same solute
            for t in self.atoms_number_offset:
                if i >= t[0] and i < t[1]:
                    break
            if p1[3]-1 < t[0] or p1[3]-1 >= t[1]:       # `data` starts from 1
                print(f'Note: bond not in the same solute: ignored: {p1}')
                dpars.append(([],'Ignored'))
                continue
            
            a = '{:<2}'.format(p1[0])
            b = '{:<2}'.format(p2[0])
            bo = False
            for p in self.oplsaa_bond_pars:
                if (a == p[0] and b == p[1]) or (a == p[1] and b == p[0]):
                    bo = True
                    dpars.append((p,'Exactly match'))
                    break
            if bo: continue

            # `a` and `b` have the same weight, they cannot be simultaneously changed
            sa = a
            sb = b
            if a.lower() in self.amber_to_oplsaa_atom_type:
                # try `new a` and `original b`
                a = self.amber_to_oplsaa_atom_type[a.lower()]
                for p in self.oplsaa_bond_pars:
                    if (a == p[0] and b == p[1]) or (a == p[1] and b == p[0]):
                        bo = True
                        dpars.append((p,f'Match type-AmberSynonym: {sa}-{a}'))
                        break
                if bo: continue

                if b.lower() in self.amber_to_oplsaa_atom_type:
                    # try `original a` and `new b`
                    a = sa
                    b = self.amber_to_oplsaa_atom_type[b.lower()]
                    for p in self.oplsaa_bond_pars:
                        if (a == p[0] and b == p[1]) or (a == p[1] and b == p[0]):
                            bo = True
                            dpars.append((p,f'Match type-AmberSynonym: {sb}-{b}'))
                            break
                    if bo: continue

                    # try `new a` and `new b`
                    a = self.amber_to_oplsaa_atom_type[a.lower()]
                    b = self.amber_to_oplsaa_atom_type[b.lower()]
                    for p in self.oplsaa_bond_pars:
                        if (a == p[0] and b == p[1]) or (a == p[1] and b == p[0]):
                            bo = True
                            dpars.append((p,f'Match type-AmberSynonym: {sa}-{a}, {sb}-{b}'))
                            break
                    if bo: continue
            else:
                if b.lower() in self.amber_to_oplsaa_atom_type:
                    # try `original a` and `new b`
                    b = self.amber_to_oplsaa_atom_type[b.lower()]
                    for p in self.oplsaa_bond_pars:
                        if (a == p[0] and b == p[1]) or (a == p[1] and b == p[0]):
                            bo = True
                            dpars.append((p,f'Match type-AmberSynonym: {sb}-{b}'))
                            break
                    if bo: continue
            # reset
            a = sa
            b = sb

            n1 = self.solutesdata['initial_bond_pars'][i][1]
            n2 = self.solutesdata['initial_bond_pars'][p1[3]-1][1]
            if n1 <= 0 or n1 > 54 or n2 <= 0 or n2 > 54:
                print(f'Warning: bond parameters: atomic number out of the range')
                dpars.append(([],'Ignored'))
                continue

            # estimate
            r = 0.0
            if a.lower() in self.bonded_atomic_radius:
                sa,l = self.bonded_atomic_radius[a.lower()]
                r += l
            else:
                sa = ' x'
                print(f'Warning: OPLS-AA bond parameter: not found: {a}')
            if b.lower() in self.bonded_atomic_radius:
                sb,l = self.bonded_atomic_radius[b.lower()]
                r += l
            else:
                sb = ' x'
                print(f'Warning: OPLS-AA bond parameter: not found: {b}')
            
            r1 = ATOMIC_ELECTRONEGATIVITY[n1-1]     # index starts from 1
            r2 = ATOMIC_ELECTRONEGATIVITY[n2-1]
            f = 120.0 * pow((r1*r2)/(r*r),0.75)
            dpars.append(((sa,sb,f,r),'Estimation'))

        # works only on initial and final bond are the same
        # except `bonds_variable`, not other types are allowed
        self.solutesdata['initial_oplsaa_bond_pars'] = dpars
        self.solutesdata['final_oplsaa_bond_pars'] = dpars
        lrv = []
        lrp = []
        lfv = []
        lfp = []
        lsv = []
        lsp = []
        for gi,gf in zip(
            self.solutesdata['initial_oplsaa_bond_pars'],
            self.solutesdata['final_oplsaa_bond_pars']
        ):
            gi = gi[0]
            gf = gf[0]
            if len(gi) or len(gf):
                vi = gi[2]
                vf = gf[3]
                dv = gf[2] - gi[2]
                dr = gf[3] - gi[3]
            else:           # care: empty
                vi = vf = dv = dr = 0.0
            lrv.append(vi+self.rc0*dv)
            lrp.append(vf+self.rc0*dr)
            lfv.append(vi+self.rc1*dv)
            lfp.append(vf+self.rc1*dr)
            lsv.append(vi+self.rc2*dv)
            lsp.append(vf+self.rc2*dr)
        self.solutesdata['reference']['pert_bonds_k'] = lrv
        self.solutesdata['reference']['pert_bonds_r'] = lrp
        self.solutesdata['first']['pert_bonds_k'] = lfv
        self.solutesdata['first']['pert_bonds_r'] = lfp
        self.solutesdata['second']['pert_bonds_k'] = lsv
        self.solutesdata['second']['pert_bonds_r'] = lsp

        ldt = []
        for i,v in zip(
            self.solutezmat['bonds_variable'],
            self.solutesdata['reference']['pert_bonds_k']
        ):
            g = self.solutezmat['data'][i]
            p = self.solutezmat['data'][g[3]]
            if g[1] < 0 or p[1] < 0:
                ldt.append(0.05)
                continue
            if v > 0.0:
                t = 0.029*(self.tk+50.0)/300.0*pow(300.0/v,0.5)
                ldt.append(t)
            else:
                ldt.append(0.1)
        self.solutesdata['pert_delta_bonds'] = ldt

    def calc_pert_angle_pars(self):
        """
        Logical:
            -> a, b, c
            if na:
                -> na, b, c
                if nc:
                    -> a,  b, nc
                    -> na, b, nc
                    if nb:
                        -> a,  nb, c
                        -> na, nb, c
                        -> a,  nb, nc
                        -> na, nb, nc
                if nb:
                    -> a,  nb, c
                    -> na, nb, c
                    -> a,  nb, nc
                    -> na, nb, nc
            if nc:
                -> a,  b, nc
                -> na, b, nc
                if nb:
                    -> a,  nb, c
                    -> na, nb, c
                    -> a,  nb, nc
                    -> na, nb, nc
            if nb:
                -> a,  nb, c
                -> na, nb, c
                -> a,  nb, nc
                -> na, nb, nc
        """
        if not self.nice: return
        dpars = []
        for i in self.solutezmat['angles_variable']:
            p1 = self.solutezmat['data'][i]
            p2 = self.solutezmat['data'][p1[3]-1]
            p3 = self.solutezmat['data'][p1[5]-1]

            if p1[1] < 0 or p2[1] < 0 or p3[1] < 0:
                dpars.append(([],'Ignored'))
                continue
            
            if p1[0] == 'DUM' or p2[0] == 'DUM' or p3[0] == 'DUM':
                dpars.append(([],'Ignored'))
                continue

            # check whether they are in the same solute
            for t in self.atoms_number_offset:
                if i >= t[0] and i < t[1]:
                    break
            if p1[5]-1 < t[0] or p1[5]-1 >= t[1]:       # `data` starts from 1
                print(f'Note: bond not in the same solute: ignored: {p1}')
                dpars.append(([],'Ignored'))
                continue

            a = '{:<2}'.format(p1[0])
            b = '{:<2}'.format(p2[0])
            c = '{:<2}'.format(p3[0])
            bo = False
            for p in self.oplsaa_angle_pars:
                if b != p[1]: continue
                if (a == p[0] and c == p[2]) or (a == p[2] and c == p[0]):
                    bo = True
                    dpars.append((p,'Exactly match'))
                    break
            if bo: continue

            # `a` and `c` have the same weight, they cannot be simultaneously changed
            sa = a
            sb = b
            sc = c
            if a.lower() in self.amber_to_oplsaa_atom_type:
                # try `new a`, `original b` and `original c`
                a = self.amber_to_oplsaa_atom_type[a.lower()]
                for p in self.oplsaa_angle_pars:
                    if b != p[1]: continue
                    if (a == p[0] and c == p[2]) or (a == p[2] and c == p[0]):
                        bo = True
                        dpars.append((p,f'Match type-AmberSynonym: {sa}-{a}'))
                        break
                if bo: continue

                if c.lower() in self.amber_to_oplsaa_atom_type:
                    # try `original a`, `original b` and `new c`
                    a = sa
                    c = self.amber_to_oplsaa_atom_type[c.lower()]
                    for p in self.oplsaa_angle_pars:
                        if b != p[1]: continue
                        if (a == p[0] and c == p[2]) or (a == p[2] and c == p[0]):
                            bo = True
                            dpars.append((p,f'Match type-AmberSynonym: {sc}-{c}'))
                            break
                    if bo: continue

                    # try `new a`, `original b` and `new c`
                    a = self.amber_to_oplsaa_atom_type[a.lower()]
                    c = self.amber_to_oplsaa_atom_type[c.lower()]
                    for p in self.oplsaa_angle_pars:
                        if b != p[1]: continue
                        if (a == p[0] and c == p[2]) or (a == p[2] and c == p[0]):
                            bo = True
                            dpars.append((p,f'Match type-AmberSynonym: {sa}-{a}, {sc}-{c}'))
                            break
                    if bo: continue

                    if b.lower() in self.amber_to_oplsaa_atom_type:
                        # try `original a`, `new b` and `original c`
                        a = sa
                        b = self.amber_to_oplsaa_atom_type[b.lower()]
                        c = sc
                        for p in self.oplsaa_angle_pars:
                            if b != p[1]: continue
                            if (a == p[0] and c == p[2]) or (a == p[2] and c == p[0]):
                                bo = True
                                dpars.append((p,f'Match type-AmberSynonym: {sb}-{b}'))
                                break
                        if bo: continue

                        # try `new a`, `new b` and `original c`
                        a = self.amber_to_oplsaa_atom_type[a.lower()]
                        b = self.amber_to_oplsaa_atom_type[b.lower()]
                        c = sc
                        for p in self.oplsaa_angle_pars:
                            if b != p[1]: continue
                            if (a == p[0] and c == p[2]) or (a == p[2] and c == p[0]):
                                bo = True
                                dpars.append((p,f'Match type-AmberSynonym: {sa}-{a}, {sb}-{b}'))
                                break
                        if bo: continue

                        # try `original a`, `new b` and `new c`
                        a = sa
                        b = self.amber_to_oplsaa_atom_type[b.lower()]
                        c = self.amber_to_oplsaa_atom_type[c.lower()]
                        for p in self.oplsaa_angle_pars:
                            if b != p[1]: continue
                            if (a == p[0] and c == p[2]) or (a == p[2] and c == p[0]):
                                bo = True
                                dpars.append((p,f'Match type-AmberSynonym: {sb}-{b}, {sc}-{c}'))
                                break
                        if bo: continue

                        # try `new a`, `new b` and `new c`
                        a = self.amber_to_oplsaa_atom_type[a.lower()]
                        b = self.amber_to_oplsaa_atom_type[b.lower()]
                        c = self.amber_to_oplsaa_atom_type[c.lower()]
                        for p in self.oplsaa_angle_pars:
                            if b != p[1]: continue
                            if (a == p[0] and c == p[2]) or (a == p[2] and c == p[0]):
                                bo = True
                                dpars.append((p,f'Match type-AmberSynonym: {sa}-{a}, {sb}-{b}, {sc}-{c}'))
                                break
                        if bo: continue

                if b.lower() in self.amber_to_oplsaa_atom_type:
                    # try `original a`, `new b` and `original c`
                    a = sa
                    b = self.amber_to_oplsaa_atom_type[b.lower()]
                    c = sc
                    for p in self.oplsaa_angle_pars:
                        if b != p[1]: continue
                        if (a == p[0] and c == p[2]) or (a == p[2] and c == p[0]):
                            bo = True
                            dpars.append((p,f'Match type-AmberSynonym: {sb}-{b}'))
                            break
                    if bo: continue

                    # try `new a`, `new b` and `original c`
                    a = self.amber_to_oplsaa_atom_type[a.lower()]
                    b = self.amber_to_oplsaa_atom_type[b.lower()]
                    c = sc
                    for p in self.oplsaa_angle_pars:
                        if b != p[1]: continue
                        if (a == p[0] and c == p[2]) or (a == p[2] and c == p[0]):
                            bo = True
                            dpars.append((p,f'Match type-AmberSynonym: {sa}-{a}, {sb}-{b}'))
                            break
                    if bo: continue

                    # try `original a`, `new b` and `new c`
                    a = sa
                    b = self.amber_to_oplsaa_atom_type[b.lower()]
                    c = self.amber_to_oplsaa_atom_type[c.lower()]
                    for p in self.oplsaa_angle_pars:
                        if b != p[1]: continue
                        if (a == p[0] and c == p[2]) or (a == p[2] and c == p[0]):
                            bo = True
                            dpars.append((p,f'Match type-AmberSynonym: {sb}-{b}, {sc}-{c}'))
                            break
                    if bo: continue

                    # try `new a`, `new b` and `new c`
                    a = self.amber_to_oplsaa_atom_type[a.lower()]
                    b = self.amber_to_oplsaa_atom_type[b.lower()]
                    c = self.amber_to_oplsaa_atom_type[c.lower()]
                    for p in self.oplsaa_angle_pars:
                        if b != p[1]: continue
                        if (a == p[0] and c == p[2]) or (a == p[2] and c == p[0]):
                            bo = True
                            dpars.append((p,f'Match type-AmberSynonym: {sa}-{a}, {sb}-{b}, {sc}-{c}'))
                            break
                    if bo: continue

            if c.lower() in self.amber_to_oplsaa_atom_type:
                # try `original a`, `original b` and `new c`
                a = sa
                c = self.amber_to_oplsaa_atom_type[c.lower()]
                for p in self.oplsaa_angle_pars:
                    if b != p[1]: continue
                    if (a == p[0] and c == p[2]) or (a == p[2] and c == p[0]):
                        bo = True
                        dpars.append((p,f'Match type-AmberSynonym: {sc}-{c}'))
                        break
                if bo: continue

                # try `new a`, `original b` and `new c`
                a = self.amber_to_oplsaa_atom_type[a.lower()]
                c = self.amber_to_oplsaa_atom_type[c.lower()]
                for p in self.oplsaa_angle_pars:
                    if b != p[1]: continue
                    if (a == p[0] and c == p[2]) or (a == p[2] and c == p[0]):
                        bo = True
                        dpars.append((p,f'Match type-AmberSynonym: {sa}-{a}, {sc}-{c}'))
                        break
                if bo: continue

                if b.lower() in self.amber_to_oplsaa_atom_type:
                    # try `original a`, `new b` and `original c`
                    a = sa
                    b = self.amber_to_oplsaa_atom_type[b.lower()]
                    c = sc
                    for p in self.oplsaa_angle_pars:
                        if b != p[1]: continue
                        if (a == p[0] and c == p[2]) or (a == p[2] and c == p[0]):
                            bo = True
                            dpars.append((p,f'Match type-AmberSynonym: {sb}-{b}'))
                            break
                    if bo: continue

                    # try `new a`, `new b` and `original c`
                    a = self.amber_to_oplsaa_atom_type[a.lower()]
                    b = self.amber_to_oplsaa_atom_type[b.lower()]
                    c = sc
                    for p in self.oplsaa_angle_pars:
                        if b != p[1]: continue
                        if (a == p[0] and c == p[2]) or (a == p[2] and c == p[0]):
                            bo = True
                            dpars.append((p,f'Match type-AmberSynonym: {sa}-{a}, {sb}-{b}'))
                            break
                    if bo: continue

                    # try `original a`, `new b` and `new c`
                    a = sa
                    b = self.amber_to_oplsaa_atom_type[b.lower()]
                    c = self.amber_to_oplsaa_atom_type[c.lower()]
                    for p in self.oplsaa_angle_pars:
                        if b != p[1]: continue
                        if (a == p[0] and c == p[2]) or (a == p[2] and c == p[0]):
                            bo = True
                            dpars.append((p,f'Match type-AmberSynonym: {sb}-{b}, {sc}-{c}'))
                            break
                    if bo: continue

                    # try `new a`, `new b` and `new c`
                    a = self.amber_to_oplsaa_atom_type[a.lower()]
                    b = self.amber_to_oplsaa_atom_type[b.lower()]
                    c = self.amber_to_oplsaa_atom_type[c.lower()]
                    for p in self.oplsaa_angle_pars:
                        if b != p[1]: continue
                        if (a == p[0] and c == p[2]) or (a == p[2] and c == p[0]):
                            bo = True
                            dpars.append((p,f'Match type-AmberSynonym: {sa}-{a}, {sb}-{b}, {sc}-{c}'))
                            break
                    if bo: continue
            if b.lower() in self.amber_to_oplsaa_atom_type:
                # try `original a`, `new b` and `original c`
                a = sa
                b = self.amber_to_oplsaa_atom_type[b.lower()]
                c = sc
                for p in self.oplsaa_angle_pars:
                    if b != p[1]: continue
                    if (a == p[0] and c == p[2]) or (a == p[2] and c == p[0]):
                        bo = True
                        dpars.append((p,f'Match type-AmberSynonym: {sb}-{b}'))
                        break
                if bo: continue

                # try `new a`, `new b` and `original c`
                a = self.amber_to_oplsaa_atom_type[a.lower()]
                b = self.amber_to_oplsaa_atom_type[b.lower()]
                c = sc
                for p in self.oplsaa_angle_pars:
                    if b != p[1]: continue
                    if (a == p[0] and c == p[2]) or (a == p[2] and c == p[0]):
                        bo = True
                        dpars.append((p,f'Match type-AmberSynonym: {sa}-{a}, {sb}-{b}'))
                        break
                if bo: continue

                # try `original a`, `new b` and `new c`
                a = sa
                b = self.amber_to_oplsaa_atom_type[b.lower()]
                c = self.amber_to_oplsaa_atom_type[c.lower()]
                for p in self.oplsaa_angle_pars:
                    if b != p[1]: continue
                    if (a == p[0] and c == p[2]) or (a == p[2] and c == p[0]):
                        bo = True
                        dpars.append((p,f'Match type-AmberSynonym: {sb}-{b}, {sc}-{c}'))
                        break
                if bo: continue

                # try `new a`, `new b` and `new c`
                a = self.amber_to_oplsaa_atom_type[a.lower()]
                b = self.amber_to_oplsaa_atom_type[b.lower()]
                c = self.amber_to_oplsaa_atom_type[c.lower()]
                for p in self.oplsaa_angle_pars:
                    if b != p[1]: continue
                    if (a == p[0] and c == p[2]) or (a == p[2] and c == p[0]):
                        bo = True
                        dpars.append((p,f'Match type-AmberSynonym: {sa}-{a}, {sb}-{b}, {sc}-{c}'))
                        break
                if bo: continue
            # reset
            a = sa
            b = sb
            c = sc

            n = 0
            vt = 0.0
            rt = 0.0
            for p in self.oplsaa_angle_pars:
                if p[1] == b:
                    n += 1
                    vt += p[3]
                    rt += p[4]
            if n == 0:
                v = 63.0
                r = 112.4
            else:
                v = vt / n
                r = rt / n
            dpars.append(((' x', ' x', ' x', v, r), 'Estimation'))

        # works only on initial and final bond are the same
        # except `bonds_variable`, not other types are allowed
        self.solutesdata['initial_oplsaa_angle_pars'] = dpars
        self.solutesdata['final_oplsaa_angle_pars'] = dpars

        lrv = []
        lrp = []
        lfv = []
        lfp = []
        lsv = []
        lsp = []
        for gi,gf in zip(
            self.solutesdata['initial_oplsaa_angle_pars'],
            self.solutesdata['final_oplsaa_angle_pars']
        ):
            gi = gi[0]
            gf = gf[0]
            if len(gi) or len(gf):
                vi = gi[3]
                vf = gf[4]
                dv = gf[3] - gi[3]
                dr = gf[4] - gi[4]
            else:           # care: empty
                vi = vf = dv = dr = 0.0
            lrv.append(vi+self.rc0*dv)
            lrp.append(vf+self.rc0*dr)
            lfv.append(vi+self.rc1*dv)
            lfp.append(vf+self.rc1*dr)
            lsv.append(vi+self.rc2*dv)
            lsp.append(vf+self.rc2*dr)
        self.solutesdata['reference']['pert_angles_k'] = lrv
        self.solutesdata['reference']['pert_angles_a'] = lrp
        self.solutesdata['first']['pert_angles_k'] = lfv
        self.solutesdata['first']['pert_angles_a'] = lfp
        self.solutesdata['second']['pert_angles_k'] = lsv
        self.solutesdata['second']['pert_angles_a'] = lsp

        ldt = []
        for i,v in zip(
            self.solutezmat['angles_variable'],
            self.solutesdata['reference']['pert_angles_k']
        ):
            g = self.solutezmat['data'][i]
            p = self.solutezmat['data'][g[3]]
            u = self.solutezmat['data'][g[5]]
            if g[1] < 0 or p[1] < 0 or u[1] < 0:
                ldt.append(2.0)
                continue
            if v > 0.0:
                t = 2.87*pow(80.0/v,0.5)*(self.tk+50.0)/300
                ldt.append(t)
            else:
                ldt.append(2.0)
        self.solutesdata['pert_delta_angles'] = ldt

    def calc_pert_xyzs(self):
        if not self.nice: return
        # deep copy
        mat = [i[3:9] for i in self.solutezmat['data']]
        first = [[j for j in i] for i in mat]
        second = [[j for j in i] for i in mat]

        for v in self.solutezmat['bonds_geometry']:
            if v[1] == 1:
                first[v[0]][1] += self.rc1 * (v[2]-first[v[0]][1])
                second[v[0]][1] += self.rc2 * (v[2]-second[v[0]][1])
            elif v[1] == 2:
                first[v[0]][3] += self.rc1 * (v[2]-first[v[0]][1])
                second[v[0]][3] += self.rc2 * (v[2]-second[v[0]][3])
            else:
                first[v[0]][5] += self.rc1 * (v[2]-first[v[0]][1])
                second[v[0]][5] += self.rc2 * (v[2]-second[v[0]][5])

        for i,v in enumerate(self.solutezmat['bonds_variable']):
            br = self.solutesdata['reference']['pert_bonds_r'][i]
            bf = self.solutesdata['first']['pert_bonds_r'][i]
            bs = self.solutesdata['second']['pert_bonds_r'][i]
            first[v][1] += bf - br
            second[v][1] += bs - br

        for i,v in enumerate(self.solutezmat['angles_variable']):
            br = self.solutesdata['reference']['pert_angles_a'][i]
            bf = self.solutesdata['first']['pert_angles_a'][i]
            bs = self.solutesdata['second']['pert_angles_a'][i]
            first[v][3] += bf - br
            second[v][3] += bs - br

        self.solutesdata['first']['xyz'] = zmat2cor(first)
        self.solutesdata['second']['xyz'] = zmat2cor(second)

    def calc_maximum_overlap(self):
        if not self.nice: return
        if not self.maxovl: return

        cf = ClosedForm()
        rcor = self.solutesdata['reference']['xyz']

        fit = cf.calc_bestfit(vl=rcor, vr=self.solutesdata['first']['xyz'], centroid=False)
        self.solutesdata['first']['xyz'] = fit

        fit = cf.calc_bestfit(vl=rcor, vr=self.solutesdata['second']['xyz'], centroid=False)
        self.solutesdata['second']['xyz'] = fit

    def recentroid(self):
        if not self.nice: return
        if not self._bo_recenter: return
        rcor = self.solutesdata['reference']['xyz']
        a = self.ncent1 - 1     # index minus 1 to start at 0
        b = self.ncent2 - 1
        offset = [(rcor[a][i]+rcor[b][i])/2.0 for i in range(3)]
        for lr,lf,ls in zip(
            self.solutesdata['reference']['xyz'],
            self.solutesdata['first']['xyz'],
            self.solutesdata['second']['xyz']
        ):
            for i in range(3):
                lr[i] -= offset[i]
                lf[i] -= offset[i]
                ls[i] -= offset[i]

    def run(self):
        self.get_xyzs()
        self.get_pars()
        self.calc_neighbors()
        self.calc_nonbond_pars()
        self.calc_sasa()
        self.calc_atoms_symmetry_list()
        self.calc_pert_dihedral_pars()
        self.calc_pert_bond_pars()
        self.calc_pert_angle_pars()
        self.calc_pert_xyzs()
        self.calc_maximum_overlap()
        self.recentroid()
        self.calc_onto_z_axis()
        if self.nice:
            self.qmidx = [
                (i,v)
                for i,v in enumerate(self.solutesdata['reference']['atomtypes'])
                if v not in [-1,99,0,2]
            ]

    def calc_onto_z_axis(self):
        if not self.nice: return
        if not self.izlong: return
        cor = self.solutesdata['reference']['xyz']
        sol1 = self.solutesdata['first']['xyz']
        sol2 = self.solutesdata['second']['xyz']
        fin0, fin1, fin2 = align_onto_z_axis(cor,sol1,sol2)  # align is only based on reference
        self.solutesdata['reference']['xyz'] = fin0
        self.solutesdata['first']['xyz'] = fin1
        self.solutesdata['second']['xyz'] = fin2

    def get_qm_xyzs(self):
        refer = [[v,*self.solutesdata['reference']['xyz'][i]] for i,v in self.qmidx]
        first = [[v,*self.solutesdata['first']['xyz'][i]] for i,v in self.qmidx]
        second = [[v,*self.solutesdata['second']['xyz'][i]] for i,v in self.qmidx]
        return refer, first, second

    def calc_atoms_symmetry_list(self):
        self.atoms_symmetry_list = []
        for i,v in enumerate(self.solutesdata['reference']['atomtypes']):
            if v in [-1,99,0,2]: continue
            if v in [1,8,9,17]: continue
            for a in [1,8,9,17]:
                ls = [
                    j for j in self.nbor[i] if self.solutesdata['reference']['atomtypes'][j] == a
                ]
                if len(ls) > 1:
                    self.atoms_symmetry_list.append(ls)

    def symmetrize_charges(self,charges,key='reference',fullpars=True):
        if fullpars:
            newcharges = [i for i in charges]
        else:
            newcharges = [0.0 for i in range(len(self.solutezmat['data']))]
            i = 0
            for j,p in enumerate(self.solutesdata['initial_bond_pars']):
                if p[1] > 0:
                    newcharges[j] = charges[i]
                    i += 1
        for t in self.solutesdata[key]['atomtypes']:
            if t in [-1,99,0,2]: continue
            for v in self.atoms_symmetry_list:
                tot = [newcharges[j] for j in v]
                a = sum(tot)/len(v)
                for j in v:
                    newcharges[j] = a
        return newcharges

    def calc_coulombic_pars(self,charges,key='reference',fullpars=True):
        new = self.symmetrize_charges(charges,key,fullpars)
        self.solutesdata[key]['Q'] = [i*self.qmscale*self.dsqesq for i in new]

    def set_pert_charges(self,charges,fullpars=True):
        news = charges if len(charges) == 3 else [charges, charges, charges]
        for i,key in enumerate(['reference','first','second']):
            self.calc_coulombic_pars(news[i],key,fullpars)

    def set_pert_energies(self,energies):
        news = energies if len(energies) == 3 else [energies, energies, energies]
        for i,key in enumerate(['reference','first','second']):
            self.solutesdata[key]['energy'] = news[i]

    def calc_total_molweights(self):
        self.total_molecular_weights = sum([
            ATOMIC_WEIGHTS[i-1] for i in self.solutesdata['reference']['atomtypes']
            if i > 0
        ])


def calc_dipoles(xyz,ql):
    qp = [0.0 for t in range(6)]
    r0 = xyz[0]
    dxx = dyy = dzz = 0.0
    n = len(xyz)
    for i in range(2,n):
        ri = xyz[i]
        q = ql[i] * 0.0548771714
        ri0 = [ri[t]-r0[t] for t in range(3)]
        dxx += ri0[0] * q
        dyy += ri0[1] * q
        dzz += ri0[2] * q
        qp[0] += ri0[0] * ri0[0] * q
        qp[1] += ri0[1] * ri0[1] * q
        qp[2] += ri0[2] * ri0[2] * q
        qp[3] += ri0[0] * ri0[1] * q
        qp[4] += ri0[0] * ri0[2] * q
        qp[5] += ri0[1] * ri0[2] * q
    dipole = pow(dxx*dxx+dyy*dyy+dzz*dzz, 0.5) * 4.802813198
    # xx, yy, zz, xy, xz, yz
    quadrupole = [i*4.802813198 for i in qp]
    return dipole, quadrupole


