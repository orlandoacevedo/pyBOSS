from PyBOSS.utils import randint, genatom
from .constants import SOLVENT_MODE, SOLVENT_PARS
from .inter import wxpot

import os
import math


def get_solvents_data(svmod1,svmod2,bond_pars):
    g_c_esq = 332.06
    solvent1_pars = SOLVENT_PARS[svmod1.lower()]
    solvent2_pars = SOLVENT_PARS[svmod2.lower()]
    s_mod1 = solvent1_pars[1]
    s_mod2 = solvent2_pars[1]

    spars = {}
    if s_mod1 == 1 or s_mod2 == 1:
        spars = SOLVENT_MODE[1]
    elif s_mod1 == 2 or s_mod2 == 2:
        spars = SOLVENT_MODE[2]

    solventpars = {
        '1': {
            'name': solvent1_pars[0],
            'modenum': solvent1_pars[1],
            'atom_names': solvent1_pars[2],
            'atom_indexes': solvent1_pars[3],
            'natoms': len(solvent1_pars[2])
        },
        '2': {
            'name': solvent2_pars[0],
            'modenum': solvent2_pars[1],
            'atom_names': solvent2_pars[2],
            'atom_indexes': solvent2_pars[3],
            'natoms': len(solvent2_pars[2])
        },
    }
    solventpars.update(spars)
    solventpars['natmx'] = max(len(solvent1_pars[2]),len(solvent2_pars[2]))
    solventpars['isvaty'] = [solvent1_pars[3], solvent2_pars[3]]

    sqrtesq = pow(g_c_esq,0.5)
    for s in ['1','2']:
        p = solventpars[s]
        bp = []
        for i in p['atom_indexes']:
            for t in bond_pars:
                if t[0] == i:
                    bp.append(t)
                    break
        p['bond_pars'] = bp
        p['QW'] = [l[3]*sqrtesq for l in bp]
        p['AW'] = [pow(4.0*l[5]*(l[4]**12),0.5) for l in bp]
        p['BW'] = [pow(4.0*l[5]*(l[4]**6), 0.5) for l in bp]

    inters = {}
    n = len(solventpars['1']['QW'])
    for i in range(n):
        q = solventpars['1']['QW'][i]
        a = solventpars['1']['AW'][i]
        b = solventpars['1']['BW'][i]
        for j in range(i,n):
            qq = solventpars['1']['QW'][j] * q
            aa = solventpars['1']['AW'][j] * a
            bb = solventpars['1']['BW'][j] * b
            ka = f'{i}-{j}'
            kb = f'{j}-{i}'
            inters[ka] = inters[kb] = [qq,aa,bb]
    solventpars['1-1'] = inters

    return solventpars


def read_waterboxes_tip4p(filename):
    """
        format:
            header: (I4,3F12.8)
            data:   (6F13.8)
    """
    if not os.path.isfile(filename):
        print(f'Warning: not a file: {filename}')
        return {}
    waterboxes = {}
    with open(filename,'rt') as f:
        while True:
            line = f.readline()
            if not line: break
            try:
                nmol = int(line[:4])
                bx = float(line[4:16])
                by = float(line[16:28])
                bz = float(line[28:40])
            except (ValueError,TypeError):      # TypeError when list empty
                print(f'Fatal: TIP4P waterbox: wrong input: {line}')
                return waterboxes
            data = []
            for l in range(nmol*4*3//6):
                line = f.readline()
                try:
                    a = float(line[:13])
                    b = float(line[13:26])
                    c = float(line[26:39])
                    d = float(line[39:52])
                    e = float(line[52:65])
                    t = float(line[65:78])
                except (ValueError,TypeError):
                    print(f'Fatal: TIP4P waterbox: wrong input: {line}')
                    return {}
                else:
                    data.extend([a,b,c,d,e,t])
            box = []
            s = 0
            for i in range(4):
                rl = []
                s += 3 * nmol
                for j in range(3):
                    u = s - (3-j) * nmol
                    rl.append([data[u+k] for k in range(nmol)])
                box.append(rl)
            mols = []
            for i in range(nmol):
                mols.append([
                    [box[k][0][i], box[k][1][i], box[k][2][i]] for k in range(4)
                ])
            waterboxes[nmol] = {
                'edge' : [bx,by,bz],
                'xyz' : mols
            }
    return waterboxes


class InitSolvents:
    def __init__(
        self, svmod1=None, svmod2=None, bond_pars=None,
        waterboxfile=None, ibox=None, irn=None, nmol=None, nmol2=None,
        *args, **kws
    ):
        self.nice = True
        self.info = ''

        waterboxes = read_waterboxes_tip4p(waterboxfile)

        self.ibox = ibox if ibox else 216
        self.waterbox = waterboxes[self.ibox]
        self.irn = irn if irn else 786123
        self.nmol = nmol
        self.nmol2 = nmol2
        self.solventsdata = get_solvents_data(svmod1, svmod2, bond_pars)

        #!!kws
        self.kws = kws
        self.kws['ac'] = self.waterbox['xyz']
        self.kws['edg2'] = e = self.waterbox['edge']
        self.kws['edge'] = [e[0]*2, e[1]*2, e[2]*2]
        self.kws['aw'] = []
        self.kws['bw'] = []
        self.kws['qw'] = []
        for k in ['1','2']:
            self.kws['aw'].append(self.solventsdata[k]['AW'])
            self.kws['bw'].append(self.solventsdata[k]['BW'])
            self.kws['qw'].append(self.solventsdata[k]['QW'])
        self.kws.update(self.solventsdata)
        self.kws['movtyp'] = 1
        self.kws['ncutat'] = None
        self.kws['nmov'] = None
        self.kws['asol'] = None
        self.kws['asol1'] = None
        self.kws['asol2'] = None
        self.kws['icut'] = 2        # tmp

    def run(self,key=None):
        self.buildup()
        self.boss_type_filter()
        self.replace_solvent_atoms()
        if self.solventsdata['1']['modenum'] <= 2:
            self.fix_solvent_waterbox()
        self.solventsdata.update(self.waterbox)

        #self.kws.pop('icut')                # keep it
        for i in self.nmr:
            self.kws['nstyp'][i] = 1
        self.kws['ac'] = self.waterbox['xyz']
        self.kws['nmr'] = self.nmr
        self.kws['irn'] = self.irn

    def buildup(self,key=None):
        if key is None: key = 'reference'
        if key != 'reference': assert False
        if key:
            self.energies = []
            for nm in range(self.ibox):
                self.kws['nm'] = nm
                gkw = wxpot(**self.kws)
                self.energies.append(gkw['elj'][0])
        else:
            assert False

    def boss_type_filter(self):
        # to maximumly keep consistence
        esw = self.energies
        sortedndx = list(range(self.ibox))
        for i in range(self.ibox):
            for j in range(i):
                if esw[j] >= esw[i]:
                    e = esw.pop(i)
                    esw.insert(j,e)
                    n = sortedndx.pop(i)
                    sortedndx.insert(j,n)
        # to keep consistent
        newndx = sorted([sortedndx[i] for i in range(self.nmol)])
        # update
        self.waterbox['xyz'] = [self.waterbox['xyz'][i] for i in newndx]
        return sortedndx

    def replace_solvent_atoms(self,irn=None):
        if not irn: irn = self.irn
        natoms2 = len(self.solventsdata['2']['atom_names'])
        irn, nmr = randint(self.nmol,self.nmol2,irn)
        anew = [[0.0, 0.0, 0.0] for i in range(50)]
        anew[1][0] = 0.945
        anew[2][0] = -0.45374566
        anew[2][1] = -1.35610283
        for i in nmr:
            xyz = self.waterbox['xyz'][i][0]
            for j in range(1,natoms2):
                for k in range(3):
                    self.waterbox['xyz'][i][j][k] = anew[j][k] + xyz[k]
        self.irn = irn
        self.nmr = self.waterbox['nmr'] = nmr

    def fix_solvent_waterbox(self):
        rcd = 0.9572
        thbcd = 104.52 / 180.0 * math.pi
        rcm = 0.15
        th2 = thbcd / 2.0
        wbox = self.waterbox['xyz']
        for i in range(self.nmol):
            if i in self.nmr: continue
            a = wbox[i][1]
            b = wbox[i][0]
            c = wbox[i][2]
            d = genatom(a,b,c,rcd,thbcd,0.0)
            a = genatom(d,c,b,rcd,thbcd,0.0)
            c = genatom(d,a,b,rcd,thbcd,0.0)
            for j in range(3):
                wbox[i][1][j] = a[j]
                wbox[i][2][j] = c[j]
                if len(wbox[i]) == 4:
                    e = genatom(d,a,b,rcm,th2,0.0)
                    wbox[i][3][j] = e[j]
        self.waterbox['xyz'] = wbox






