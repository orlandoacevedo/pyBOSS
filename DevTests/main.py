from .utils import zmat2cor
from .constants import SOLVENT_MODE, SOLVENT_PARS
from .init_pars import init_pars
from .get_pars import read_config_file
from .init_solutes import InitSolutes
from .init_solvents import InitSolvents
from .interactions import Interactions
from .read_zmat import read_zmat
from .read_ff_pars import read_ff_file, read_oplsaa_file
from .qm import qm

import os
import math
import collections
import json
import copy

SRCPATH    = os.path.dirname(__file__)        # this will be udpated
OPLSAA_PAR  = os.path.join(SRCPATH,'data','oplsaa.par')
OPLSAA_SB   = os.path.join(SRCPATH,'data','oplsaa.sb')
ORG1BOX     = os.path.join(SRCPATH,'data','org1box')
ORG2BOX     = os.path.join(SRCPATH,'data','org2box')
WATBOX      = os.path.join(SRCPATH,'data','watbox')


def allclose(a, b, tol=None) -> bool:
    """recursively compare a and b on given tolerance"""
    tol = tol if tol else 0.000001
    if a is None and b is None:
        return True
    if isinstance(a, bool) and isinstance(b, bool):
        return a == b
    if isinstance(a, (int,float)) and isinstance(b, (int,float)):
        if abs(a-b) > tol:
            return False
        return True
    if isinstance(a, (tuple,list)) and isinstance(b, (tuple,list)):
        if len(a) != len(b):
            return False
        if not a:           # for empty object
            return True
        for vi,vj in zip(a,b):
            if not allclose(vi, vj, tol):
                return False
        return True
    if isinstance(a, dict) and isinstance(b, dict):
        if set(list(a.keys())) != set(list(b.keys())):
            return False
        for k in a.keys():
            if not allclose(a[k], b[k], tol):
                return False
        return True
    return False


def my_print_tld(dic,offset=0):
    """powerful print for `tuple,list,dictionary`"""
    skip = ' '*offset
    if isinstance(dic,(int,float,str)):
        print(f'{skip}{dic}')
    elif isinstance(dic,(list,tuple,set)):
        for i in dic:
            print(f'{skip}{i}')
    elif isinstance(dic,dict):
        for k,v in dic.items():
            print(skip + '='*4)
            print(f'{skip}{k}:')
            my_print_tld(v,offset+4)
    else:
        print(f'{skip} "Unknown Type PASS"')
    print()


def my_dump_json(dic,file='checkrst.json'):
    print(f'Note: Dumped to {file}')
    f = open(file,'wt')
    json.dump(dic,f,indent=4)
    f.close()


CWD = os.path.dirname(__file__)
CONFIG = os.path.join(CWD,'configs.txt')
ZMAT = os.path.join(CWD,'pmfzmat')


dpars = read_config_file(CONFIG)
dpars = init_pars(dpars)


solutezmat = read_zmat(ZMAT)

bond_pars, dihedral_pars = read_ff_file(OPLSAA_PAR,skip=0,msg=False)
oplsaa_bond_pars, oplsaa_angle_pars = read_oplsaa_file(OPLSAA_SB,msg=False)


gkws = copy.deepcopy(dpars)

Xsolute = InitSolutes(
    bond_pars=bond_pars, dihedral_pars=dihedral_pars,
    oplsaa_bond_pars=oplsaa_bond_pars, oplsaa_angle_pars=oplsaa_angle_pars,
    solutezmat=solutezmat,
    **gkws
)
Xsolute.run()


from .interface import Interface
Gsolutesdata = g = Xsolute.solutesdata
Xif = Interface(
    refer_xyz=g['reference']['xyz'],refer_atnum=g['reference']['atomtypes'],
    first_xyz=g['first']['xyz'],first_atnum=g['first']['atomtypes'],
    second_xyz=g['second']['xyz'],second_atnum=g['second']['atomtypes'],
)
Xif.run()
charges = [Xif.refer_charges, Xif.first_charges, Xif.second_charges]
energies = [Xif.refer_energy, Xif.first_energy, Xif.second_energy]
Xsolute.set_pert_charges(charges,fullpars=False)
Xsolute.set_pert_energies(energies)




#!!gkws
gkws['iztyp'] = [g[1] for g in Xsolute.solutezmat['data']]
gkws['nsatm'] = Xsolute.total_number_atoms
gkws['ityp'] = [i for i in range(gkws['nsatm'])]      # one-on-one correspond with `a`, `b`, `q`
gkws['ityp1'] = [i for i in range(gkws['nsatm'])]     # 
gkws['ityp2'] = [i for i in range(gkws['nsatm'])]     # 
gkws['a'] = [Xsolute.solutesdata['reference']['A'], Xsolute.solutesdata['first']['A'], Xsolute.solutesdata['second']['A']] 
gkws['b'] = [Xsolute.solutesdata['reference']['B'], Xsolute.solutesdata['first']['B'], Xsolute.solutesdata['second']['B']] 
gkws['q'] = [Xsolute.solutesdata['reference']['Q'], Xsolute.solutesdata['first']['Q'], Xsolute.solutesdata['second']['Q']] 
gkws['anew'] = Xsolute.solutesdata['reference']['xyz']
gkws['anew1'] = Xsolute.solutesdata['first']['xyz']
gkws['anew2'] = Xsolute.solutesdata['second']['xyz']
gkws['nstyp'] = [0 for i in range(1024)]             # index of solvent, `MXS=2500`
gkws['isolute'] = p = []
for c,g in enumerate(Xsolute.solutesdata['atoms_number_offset']):
    p.extend([c for i in range(*g)])


Xsolvent = InitSolvents(
    movetype=1, waterboxfile=WATBOX,
    bond_pars=bond_pars,
    solutezmat=solutezmat,
    solutesdata=Gsolutesdata,
    **gkws
)
Xsolvent.run('reference')
Gsolventsdata = Xsolvent.solventsdata
gkws.update(Xsolvent.kws)




nmol = dpars['nmol']

nvchg = 5
nschg = 3
nconsv = 0

vdel = 0.0
vdel2 = 2.0 * vdel
nsolute = 1
isolec = 1
if isolec > nsolute: isolec = 1
gkws['isolec'] = isolec

icut = dpars['icut']
if icut == 3:
    ncutat = 2
    gkws['icutat'][0], gkws['icutat'][1] = dpars['ncent1'], dpars['ncent2']
elif icut == 4:
    ls = Gsolutesdata['reference']['atomtypes']
    ncutat = len(ls) - ls.count(0)
elif icut == 5:
    ls = Gsolutesdata['reference']['atomtypes']
    ncutat = len(ls) - ls.count(0) - ls.count(-1) - ls.count(1)
else:
    ncutat = 0
gkws['ncutat'] = ncutat

# added keys
gkws['modsv1'] = 2
gkws['modsv2'] = 3
if gkws['vdel'] <= 0.0:
    if gkws['modsv1'] > 2 and gkws['modsv1'] != 3:
        gkws['vdel'] = float(gkws['nmol']//10*15)



from .interface import Interface

it = Interface(
    refer_xyz=Gsolutesdata['reference']['xyz'],
    refer_atnum=Gsolutesdata['reference']['atomtypes'],
    first_xyz=Gsolutesdata['first']['xyz'],
    first_atnum=Gsolutesdata['first']['atomtypes'],
    second_xyz=Gsolutesdata['second']['xyz'],
    second_atnum=Gsolutesdata['second']['atomtypes'],
)
it.run()
Gsolutesdata['reference']['Q'] = it.refer_charges
Gsolutesdata['first']['Q'] = it.first_charges
Gsolutesdata['second']['Q'] = it.second_charges



dipole = []
quadrupole = []
n = len(Gsolutesdata['reference']['xyz'])
for k in ['reference','first','second']:
    qp = [0.0 for t in range(6)]
    r0 = Gsolutesdata[k]['xyz'][0]
    dxx = dyy = dzz = 0.0
    j = -1
    for i in range(2,n):
        if Gsolutesdata[k]['atomtypes'][i] <= 0: continue
        j += 1
        ri = Gsolutesdata[k]['xyz'][i]
        q = Gsolutesdata[k]['Q'][j] * 0.0548771714
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
    quadrupole.append(qp)
    dipole.append(pow(dxx*dxx+dyy*dyy+dzz*dzz, 0.5))
dipole = [i*4.802813198 for i in dipole]
# xx, yy, zz, xy, xz, yz
quadrupole = [[i*4.802813198 for i in j] for j in quadrupole]


# need more debug -- hydrogen donor & acceptor
xyz = Gsolutesdata['reference']['xyz']
atomtypes = Gsolutesdata['reference']['atomtypes']
a = Gsolutesdata['reference']['A']
numdonor = numacceptor = 0
for rxyz,st,a in zip(xyz,atomtypes,a):
    #                                 heteroatom only if hydrogen
    if st in [7,8,16] or (st == 1 and a <= 0.0):
        for i,wxyz in enumerate(Gsolventsdata['xyz']):
            if i in Gsolventsdata['nmr']:
                bp = Gsolventsdata['2']['bond_pars']
            else:
                bp = Gsolventsdata['1']['bond_pars']
            sn = len(bp)
            for j in range(len(bp)):
                sp = bp[j][2]
                mxyz = xyz[j]
                if a == 1:
                    if sp in [7,8,16]:
                        d = sum([rxyz[t]-mxyz][t] for t in range(3))
                        if d < 6.25:
                            numdonor += 1
                elif sp == 1:
                    d = sum([rxyz[t]-mxyz][t] for t in range(3))
                    if d < 6.25:
                        numacceptor += 1


from .inter import enmtx

ckws = enmtx(**gkws)








wiks = [1.0/(i+wkc) for i in wiks]
total = sum(wiks)
wiks = [i/total for i in wiks]
maxwik = max(wiks)


ecut = inter.ecut()

radsa = Gsolutesdata['radius_sa']
asol = Gsolutesdata['reference']['xyz']

from boss_calc import PyBOSSCalc
pbc = PyBOSSCalc(radsa=radsa, asol=asol, rsolv=dpars['rsolv'])

sasa = pbc.calc_amber_sasa(Gsolutesdata['amber_sasa_atomtypes'])

Gsolutesdata['sasa'] = sasa

from utils import randu

x,y = randu(dpars['irn'])
dpars['irn'] = x





last = 0
nfl = 0
nsolmv = 0
nflptr = 0
nflprj = 0
iflip = 0
xmol = float(nmol)
xmlbet = (xmol + float(nsolut)) / beta
epin = 1.0 / eprinc
dince = 1.0 / edinc
epins = 1.0 / essinc
dinces = 1.0 / ebsinc
ihist = 0
mhist = mxcon // (120 * nschg) + 1

ndihmx = max(ntdih, 8)

tk20 = 20.0 / beta
tk25 = 25.0 / beta
nschg5 = 5 * nschg
nsolup = 20 * nschg
nsa = 0
nfrqsa = min((1 + mxcon // 10), nsolup)
elow = 1e20

if nbuse == 1:
    rncsq = (rcut + 1.5)**2
    #call neibor()


# 10
ncon += 1
movvol = 0
movsol = 0
movtyp = 0
iflip = 0

if iewald == 1:
    raise NotImplementedError()

if nsolut > 2:
    if ncon % nsolup == 0:
        for l in range(3, nsolut):
            if nsoltr[l] > 0:
                fac = float(nsolac[l]) / float(nsoltr[l])
                if fac > 0.6:
                    rdelas[l] = min(rdelas[l] * 1.1, 100.0)
                    adelas[l] = min(adelas[l] * 1.1, twopi)
                elif fac < 0.4:
                    rdelas[l] = rdelas[l] * 0.9
                    adelas[l] = adelas[l] * 0.9


if ncon % nschg5 == 0:
    #call dipole(dip, qdr)

    ndicnt += 1
    dipsum[0] += dip[0]
    dipsum[1] += dip[1]
    dipsum[2] += dip[2]

    if nmol != 0:
        #call hbond(nacp, ndon)
        nhbcnt += 1
        nhbnda += nacp
        nhbndd += ndon


if ncon % nfrqsa == 0:
    nsa += 1
    if nosolu == 0:
        #call savol2(0)
        avsa += sa
        avvl += vl
        for i in range(4):
            avsatp[i] += satp[i]
    #call ssljco()
    eslj += esljol
    esco += escool


bo40 = True
bo50 = False
bo80 = False
if not (nvchg != 999999 and ncon % nvchg == 0):
    # !! goto 40
    bo40 = False
    if ncon % nschg == 0:
        # goto 50
        bo50 = True
    else:
        if nopref == 0:
            while True:
                xmov = xmol * ranu() + 1.0
                nmov = int(xmov)
                x = ranu()
                wr = wi[nmov] / wimax
                if wr > x: break
        else:
            xmov = xmol * ranu() + 1.0
            nmov = int(xmov)

        movtry[nmov] += 1
        #call movsvn()
        # goto 80
        bo80 = True

bo70 = False
if bo40 and not bo50 and not bo80:
    # 40
    movvol = 1
    nmov = int(xmol * ranu() + 1.0)
    # goto 70
    bo70 = True

if bo50 and not bo70 and not bo80:
    # 50
    movtyp = 1
    nmov = 0
    movsol = 1
    nsolmv += 1
    if ntdih != 0:
        if nsolmv % mhist == 0:
            ihist += 1
            for i in range(1, ndihmx + 1):
                histdi[ihist][i] = 1 + int(48.0 * phi[i] / twopi)


if bo70 and not bo80:
    # 70
    #call movmol()

    if movsol == 1:
        nsotr[isol] += 1
        if abs(enbne-999999.0) < 0.1 or abs(enbne1-999999.0) < 0.1 or abs(enbne2-999999.0) < 0.1:
            print('Error: outside scripts failed')
            enbne = 999999.0
            enbne1 = 999999.0
            enbne2 = 999999.0

        if nflptr < 0:
            print('SUBROUTINE MONTE NFLPTR < 0')
            nflptr = -nflptr
            iflip = 1


# 80
delh = enew - eold
headfile = '~+`#@!'
xztype = '~+%#@!'
if headfile == 'HEAD  ':
    with open(headfile, 'r') as f:
        xztemp, xztype = f.readline().split()

fac = 1.0
if movtyp != 0:
    if movsol != 1:
        delh += pvcon * (vnew - vold) - xmlbet * math.log(vnew/vold)

    if delh <= 0.0:
        if movsol == 1:
            if xztype == 'FALSE' or xztype == 'TRUE ':
                with open(headfile, 'a') as f:
                    if qmname == 'G09U':
                        if xztype == 'FALSE':
                            f.write('Error: MONTE Bug, -5-1-\n')
                        else:
                            f.write('GAUSS_NORM_WIN\n')
                    elif qmname == 'G091':
                        if xztype == 'FALSE':
                            f.write('Error: MONTE Bug, -5-2-\n')
                        else:
                            f.write('AEGAU_NWNW\n')
                    elif qmname == 'G092':
                        if xztype == 'FALSE':
                            f.write('Error: MONTE Bug, -5-3-\n')
                        else:
                            f.write('AEGAU_TWNW\n')
                    elif qmname == 'G09D':
                        if xztype == 'TRUE ':
                            f.write('Error: MONTE Bug, -5-5-\n')
                        else:
                            f.write('AENET_TLNW, small train problem\n')
                        qmname = 'G09U'
                        # goto 550
                    elif qmname == 'G09L':
                        if xztype == 'TRUE ':
                            f.write('Error: MONTE Bug, -5-6-\n')
                        else:
                            f.write('AENET_NLNW, big train problem\n')
                        qmname = 'G09U'
                        # goto 550
                    else:
                        f.write('Error: MONTE Bug, -5-4-\n')
            qmname = 'G09U'
        # goto 90
else:
    if nopref == 0:
        fac = wnew[nmov] / wi[nmov]


if delh > tk20:
    if movsol == 1:
        if xztype == 'FALSE' or xztype == 'TRUE ':
            with open(headfile, 'a') as f:
                if qmname == 'G09U':
                    if xztype == 'FALSE':
                        f.write('AENET_NLNL\n')
                    else:
                        f.write('GAUSS_NORM_LOSE\n')
                elif qmname == 'G091':
                    if xztype == 'FALSE':
                        f.write('Error: MONTE Bug, -1-1-\n')
                    else:
                        f.write('AEGAU_NWNL, should not happen\n')
                elif qmname == 'G092':
                    if xztype == 'FALSE':
                        f.write('Error: MONTE Bug, -1-2-\n')
                    else:
                        f.write('AEGAU_TWNL\n')
                elif qmname == 'G09D':
                    if xztype == 'TRUE ':
                        f.write('Error: MONTE Bug, -1-4-\n')
                    else:
                        f.write('AENET_TLNL\n')
                elif qmname == 'G09L':
                    if xztype == 'TRUE ':
                        f.write('Error: MONTE Bug, -1-6-\n')
                    else:
                        f.write('AENET_NLNL\n')
                else:
                    f.write('Error: MONTE Bug, -1-3-\n')
        qmname = 'G09U'
    # goto 550
elif delh < -tk20:
    if movsol == 1:
        if xztype == 'FALSE' or xztype == 'TRUE ':
            with open(headfile, 'a') as f:
                if qmname == 'G09U':
                    if xztype == 'FALSE':
                        f.write('Error: MONTE Bug, -2-1-\n')
                    else:
                        f.write('GAUSS_NORM_WIN\n')
                elif qmname == 'G091':
                    if xztype == 'FALSE':
                        f.write('Error: MONTE Bug, -2-2-\n')
                    else:
                        f.write('AEGAU_NWNW\n')
                elif qmname == 'G092':
                    if xztype == 'FALSE':
                        f.write('Error: MONTE Bug, -2-3-\n')
                    else:
                        f.write('AEGAU_TWNW\n')
                elif qmname == 'G09D':
                    if xztype == 'TRUE ':
                        f.write('Error: MONTE Bug, -2-5-\n')
                    else:
                        f.write('AENET_TLNW, small training problem\n')
                    qmname = 'G09U'
                elif qmname == 'G09L':
                    if xztype == 'TRUE ':
                        f.write('Error: MONTE Bug, -2-6-\n')
                    else:
                        f.write('AENET_NLNW, big training problem\n')
                    qmname = 'G09U'
                else:
                    f.write('Error: MONTE Bug, -2-4-\n')
        qmname = 'G09U'
    # goto 90


x = ranu()
xbeta = beta

if movsol == 1:
    if isol == lhtsol:
        xbeta = betlht
    if lhtsol == 9999:
        xbeta = betlht

emet = fac * math.exp(-xbeta * delh)
if emet < x:
    if movsol == 1:
        if xztype == 'FALSE' or xztype == 'TRUE ':
            with open(headfile, 'a') as file:
                if qmname == 'G09U':
                    if xztype == 'FALSE':
                        file.write('Error: MONTE Bug, -3-7-\n')
                    else:
                        file.write('GAUSS_TOSS_LOSE\n')
                elif qmname == 'G091':
                    if xztype == 'FALSE':
                        file.write('Error: MONTE Bug, -3-1-\n')
                    else:
                        file.write('AEGAU_NWTL\n')
                elif qmname == 'G092':
                    if xztype == 'FALSE':
                        file.write('Error: MONTE Bug, -3-2-\n')
                    else:
                        file.write('AEGAU_TWTL\n')
                elif qmname == 'G09D':
                    if xztype == 'TRUE ':
                        file.write('Error: MONTE Bug, -3-4-\n')
                    else:
                        file.write('AENET_TLTL, lucky\n')
                elif qmname == 'G09L':
                    if xztype == 'TRUE ':
                        file.write('Error: MONTE Bug, -3-6-\n')
                    else:
                        file.write('AENET_NLTL, small training problem\n')
                else:
                    file.write('Error: MONTE Bug, -3-3-\n')
            qmname = 'G09U'
    # goto 550


if movsol == 1:
    if xztype == 'FALSE' or xztype == 'TRUE ':
        with open(headfile, 'a') as file:
            if qmname == 'G09U':
                if xztype == 'FALSE':
                    file.write('Error: MONTE Bug, -4-5-\n')
                    # goto 550
                else:
                    file.write('GAUSS_TOSS_WIN\n')
            elif qmname == 'G091':
                if xztype == 'FALSE':
                    file.write('Error: MONTE Bug, -4-1-\n')
                else:
                    file.write('AEGAU_NWTW\n')
            elif qmname == 'G092':
                if xztype == 'FALSE':
                    file.write('Error: MONTE Bug, -4-2-\n')
                else:
                    file.write('AEGAU_TWTW\n')
            elif qmname == 'G09D':
                if xztype == 'TRUE ':
                    file.write('Error: MONTE Bug, -4-4-\n')
                else:
                    file.write('AENET_TLTW, really unlucky\n')
                qmname = 'G09U'
                # goto 550
            elif qmname == 'G09L':
                if xztype == 'TRUE ':
                    file.write('Error: MONTE Bug, -4-6-\n')
                else:
                    file.write('AENET_NLTW, little training problem\n')
                qmname = 'G09U'
                # goto 550
            else:
                file.write('Error: MONTE Bug, -4-3-\n')
        qmname = 'G09U'


# 90
if ncon != 1:
    if movsol == 1:
        headfile = '~+`#@!'
        headfile = os.getenv('HEAD', headfile)
        if headfile == 'HEAD  ':
            with open(headfile, 'a') as file:
                file.write(f'Accepted Energy Solute Ref: {enbne}\n')
                file.write(f'Accepted Energy Solute 1st: {enbne1}\n')
                file.write(f'Accepted Energy Solute 2nd: {enbne2}\n')
    
    xrep = float(nrepet)
    etot += xrep * eold
    x = xrep * vold
    v += x
    vsq += x * vold
    y = eold + pvcon * vold
    x = xrep * y
    h += x
    vh += x * vold
    
    if igbsa == 1:
        esx += (egb + esasa) * xrep
    else:
        esx += xrep * esonol
        esxco += xrep * esonco
        esxlj += xrep * esonlo
        esx1 += xrep * esol1
        esx2 += xrep * esol2

    einter = eold - edihol - enbol - ebndol - eangol - esbnol - esanol - esnbol - esdiol - ebcold
    x = einter + pvcon * vold
    y1 = xrep * x
    hinter += y1
    hinsq += y1 * x

    if nofep == 0:
        x = exxold + esonol + edihol + enbol + ebndol + eangol + epoold
        del1 = exxol1 + esol1 + ediol1 + enbol1 + ebnol1 + eanol1 + epool1 - x
        del2 = exxol2 + esol2 + ediol2 + enbol2 + ebnol2 + eanol2 + epool2 - x
        if iewald == 1:
            raise NotImplementedError()

        x = math.exp(-beta * del1) * xrep
        edelg += x
        udel1 += (y + del1) * x
        x = math.exp(-beta * del2) * xrep
        edelg2 += x
        udel2 += (y + del2) * x
        x = math.exp(-beta * del1 * 0.5) * xrep
        sos1 += x
        x = math.exp(-beta * del2 * 0.5) * xrep
        sos2 += x

    if nsatm > 1:
        exx += xrep * exxold
        exx1 += xrep * exxol1
        exx2 += xrep * exxol2
        epol += xrep * epoold
        epol1 += xrep * epool1
        epol2 += xrep * epool2
        exxco += xrep * exxolc
        exxlj += xrep * exxoll
        ebc += xrep * ebcold
        ebc1 += xrep * ebcol1
        ebc2 += xrep * ebcol2

    if ntdih != 0:
        edih += xrep * edihol
        edih1 += xrep * ediol1
        edih2 += xrep * ediol2

    enb += xrep * enbol
    enb1 += xrep * enbol1
    enb2 += xrep * enbol2

    if ntbnd != 0:
        ebnd += xrep * ebndol
        ebnd1 += xrep * ebnol1
        ebnd2 += xrep * ebnol2

    if ntang != 0:
        eang += xrep * eangol
        eang1 += xrep * eanol1
        eang2 += xrep * eanol2

    if icussl != 0:
        esbnd += xrep * esbnol
        esang += xrep * esanol
        esdih += xrep * esdiol
        esnb += xrep * esnbol
        esint += xrep * esinol


    if natmx != 0:
        # !! goto 140

        for i in range(nrdl):
            if nrdf != 0:
                for j in range(nrdf):
                    xr[i][j] += xrep * float(ior[i][j])

            if nrdfs != 0:
                for j in range(nrdfs):
                    xrs[i][j] += xrep * float(iors[i][j])

        i = int(dince * (eonold - edmin)) + 1
        if i >= 1:
            if i <= 50:
                xedist[i] += xrep

        for i in range(50):
            xepar[i] += xrep * float(ioepar[i])
            xess[i] += xrep * float(ioess[i])

        i = int(dinces * (esonol - ebsemin)) + 1
        if i >= 1:
            if i <= 50:
                xbesol[i] += xrep


# 140
if ntdih != 0:
    for i in range(ntdih):
        j = int(phi[i] / (6.0 * dtorad)) + 1
        if j > 60:
            j -= 60
        iphidi[i, j] += nrepet


naccpt += nrepet
if last == 1:
    # goto 690
    pass

eold = enew
vold = vnew
epoold = eponew
epool1 = epone1
epool2 = epone2
nrepet = 1

if igbsa == 1:
    raise NotImplementedError()

if iewald == 1:
    raise NotImplementedError()


if movsol != 1:
    if movvol != 1:
        nmr = int(xmol * ranu() + 1.0)
        movacp[nmov] += 1
        esonol = esone
        esonco = esonc
        esonlo = esonl
        esol1 = eson1
        ess1[nmov] = es1
        esol2 = eson2
        ess2[nmov] = es2
        
        if icussl != 0:
            raise NotImplementedError()

        j = int(epins * (ess[nmov] - essmin)) + 1
        k = int(epins * (esol - essmin)) + 1
        ess[nmov] = esol
        essc[nmov] = ecsx
        essl[nmov] = elsx

        if j <= 50 and j >= 1:
            ioess[j] -= 1
        if k <= 50 and k >= 1:
            ioess[k] += 1

        ioepar = [0 for i in range(50)]
        eonold = 0.0
        for i in range(1, nmov):
            n = nmol * (i - 1) + nmov - i - (i * (i - 1)) // 2
            eij[n] = emov[i]

        j = nmol * (nmov - 1) - nmov - (nmov * (nmov - 1)) // 2
        for i in range(nmov + 1, nmol + 1):
            n = j + i
            eij[n] = emov[i]

        for i in range(nmr):
            n = nmol * (i - 1) + nmr - i - (i * (i - 1)) // 2
            e = eij[n]
            eonold += e
            j = int(epin * (e - eprmin)) + 1
            if j <= 50 and j >= 1:
                ioepar[j] += 1

        k = nmol * (nmr - 1) - nmr - (nmr * (nmr - 1)) // 2
        for i in range(nmr + 1, nmol + 1):
            e = eij[i + k]
            eonold += e
            j = int(epin * (e - eprmin)) + 1
            if j <= 50 and j >= 1:
                ioepar[j] += 1

        ntyp = nstyp[nmov]
        natoms = nsvat[ntyp]
        for j in range(3):
            for k in range(natoms):
                x = anew[k][j]
                anew[k][j] = ac[nmov][k][j]
                ac[nmov][k][j] = x

        if not (ntyp == 2 or nrdfs == 0):
            # !! goto 260

            movtyp = 2
            if modsv1 <= 2:
                x = wxpot(0)
            elif modsv1 > 2:
                x = sxpot(0)

            movtyp = 0
            for i in range(nrdfs):
                m = knew[i]
                if m <= nrdl and m >= 1:
                    iors[m][i] += 1
                m = int(rinc * (rnew[i] - rdlmin)) + 1
                if m <= nrdl and m >= 1:
                    iors[m][i] -= 1

        # 260
        if nbuse == 1:
            nbmov[nmov] += 1
            if nbmov[nmov] == nblim:
                #call neighbor()
                pass

        if nrdf == 0:
            # goto 530
            pass

        ntyp = nstyp[nmr]
        if ntyp == 2:
            while True:
                nmr = int(xmol * ranu() + 1.0)
                ntyp = nstyp[nmr]
                if ntyp != 2: break

        natoms = nsvat[1]
        ncen = 1
        if modcus == 1:
            ncen = ncents

        for i in range(3):
            for j in range(natoms):
                anew[j][i] = ac[nmr][j][i]


        for i in range(nrdf):
            for j in range(nrdl):
                ior[j][i] = 0

        for i in range(nmol):
            if nstyp[i] == 2 or i == nmr: continue

            for j in range(3):
                for k in range(natoms):
                    anm[k][j] = ac[i][k][j]

            x2 = 0.0
            for k in range(3):
                x = anew[ncen][k] - anm[ncen][k]
                if not (islab == 1 and k == 3):
                    # !! goto 370

                    if not (x <= edg2[k]):
                        # !! goto 360
                        y = edge[k]

                        # 340
                        for j in range(natoms):
                            anm[j][k] += y
                        x = x - y
                        # goto 370
                    else:
                        # 360
                        if x < -edg2[k]:
                            y = -edge[k]
                            # goto 340
                            # 340
                            for j in range(natoms):
                                anm[j][k] += y
                            x = x - y
                            # goto 370
                # 370
                x2 += x ** 2

            if x2 <= rcutsq:
                for j in range(nrdf):
                    k = nrdfs1[j]
                    l = nrdfs2[j]
                    x2 = 0.0
                    for m in range(3):
                        x2 += (anew[k][m] - anm[l][m]) ** 2
                    x = x2 ** 0.5
                    m = int(rinc * (x - rdlnmin)) + 1
                    if m <= nrdl and m >= 1:
                        ior[m][j] += 1
        # goto 530
else:
    for i in range(nmol):
        ess[i] = emov[i]
        essc[i] = emovc[i]
        essl[i] = emovl[i]
        ess1[i] = esmov[i]
        ess2[i] = esmov2[i]


if natmx != 0:
    # !! goto 460

    esonol = esone
    esonco = esonc
    esonlo = esonl
    esol1 = eson1
    esol2 = eson2
    ioess = [0 for i in range(50)]
    for i in range(1, nmol + 1):
        j = int(epins * (ess[i] - essmin)) + 1
        if j <= 50 and j >= 1:
            ioess[j] += 1

    for j in range(nrdfs):
        for i in range(nrdl):
            iors[i][j] = idists[i][j]

# 460
if movvol != 1:
    for j in range(3):
        for i in range(nsatm):
            asol[i][j] = anew[i][j]
            asol1[i][j] = anew1[i][j]
            asol2[i][j] = anew2[i][j]

exxold = exxnew
exxol1 = exxne1
exxol2 = exxne2
exxolc = exxnec
exxoll = exxnel
ebcold = ebcnew
ebcol1 = ebcne1
ebcol2 = ebcne2
if ntbnd != 0:
    ebndol = ebndne
    ebnol1 = ebnne1
    ebnol2 = ebnne2
    for i in range(ntbnd):
        bnd[i] = bndnew[i]

if ntang != 0:
    eangol = eangne
    eanol1 = eanne1
    eanol2 = eanne2
    for i in range(ntang):
        ang[i] = angnew[i]

if ntdih != 0:
    edihol = edihne
    ediol1 = edine1
    edioll2 = edine2
    for i in range(ntdih):
        phi[i] = phinew[i]

if nbp != 0 or isqm == 1:
    enbol = enbne
    enbol1 = enbne1
    enbol2 = enbne2

nsolac[isol] += 1
if nconsv == -1:
    if enew < elow:
        # rewind(idsks)
        # writzm(idsks, enew)
        elow = enew
else:
    eonold = eone
    ioepar = [0 for i in range(50)]
    for i in range(nmol):
        if i != nmov:
            j = int(epin * (emov[i] - eprmin)) + 1
            if j <= 50 and j >= 1:
                ioepar[j] += 1

if nopref == 0:
    for i in range(nmol):
        wi[i] = wnew[i]
    wimax = wnmax
    wisum = wnsum
# goto 680


nrject += 1
nrepet += 1

if iflip == 1:
    nflprj += 1

if movsol == 1 and ichmod != 0:
    for i in range(nats):
        j = irenum[i]
        q[ityp[j] ] = qold[i][0]
        q[ityp1[j]] = qold[i][1]
        q[ityp2[j]] = qold[i][2]

if movvol == 1:
    for i in range(3):
        edge[i] = oldedge[i]
        edg2[i] = 0.5 * edge[i]

    for j in range(3):
        if ivxyz == 1 and j != ivax: continue

        for i in range(nmol):
            tmp = ac[i][0][j] * slvfac / (slvfac + 1.0)
            natoms = nsvat[nstyp[i]]
            for k in range(1, natoms + 1):
                ac[i][k][j] -= tmp

        for i in range(nsatm):
            asol [i][j] -= soldel[j]
            asol1[i][j] -= soldel[j]
            asol2[i][j] -= soldel[j]

    if noss != 1:
        # !! goto 645
        movtyp = 0
        knt = 1
        n = nmol - 1
        for i in range(n):
            k = i + 1
            jtyp = nstyp[i]
            natoms = nsvat[jtyp]
            ntyp = modsv1
            if jtyp == 2:
                ntyp = modsv2
            nmov = i
            for m in range(3):
                for l in range(natoms):
                    anew[l][m] = ac[i][l][m]

            for j in range(k, nmol):
                if nsvat[2] == 0:
                    if modsv1 > 2:
                        eij[knt] = sspot(j)
                    else:
                        eij[knt] = wwpot(j)
                else:
                    jtyp = modsv1
                    if nstyp[j] == 2:
                        jtyp = modsv2
                    if ntyp <= 2 and jtyp <= 2:
                        eij[knt] = wwpot(j)
                    else:
                        eij[knt] = sspot(j)
                knt += 1

    # 645
    movtyp = 1
    for j in range(3):
        for i in range(nsatm):
            anew[i][j] = asol[i][j]
            anew1[i][j] = asol1[i][j]
            anew2[i][j] = asol2[i][j]

    for i in range(nmol):
        modsv = modsv1
        if nstyp[i] == 2:
            modsv = modsv2
        if modsv > 2:
            ess[i] = sxpot(i)
        else:
            ess[i] = wxpot(i)
        essc[i] = ecsx
        essl[i] = elsx
        ess1[i] = es1
        ess2[i] = es2

    nfl += 1
    vnew = vold

if nconsv > 0 and nconsv != 999999:
    print('subroutine monte after 680, inside condition')
    if ncon % nconsv == 0:
        print(f'{nmol} {nsatm} {nsolut} {nvdih} {nvbnd} {nvang} {naccpt} {nrject} {irn} {iversn} {nr} {nx}')
        print(f'{t} {p} {edg2}')
        print(f'{eold} {eonold} {esonol} {esol1} {esol2} {exxold} {exxol1}')
        print(f'{exxol2} {edihol} {ediol1} {ediol2} {enbol} {enbol1} {enbol2}')
        print(f'{ebcold} {ebcol1} {ebcol2} {ebndol} {ebnol1} {ebnol2} {eangol} {eanol1} {eanol2}')
        print(f'{esonco} {esonlo} {exxolc} {exxoll}')
        print(f'{esbnol} {esanol} {esdiol} {esnbol} {esinol} {esljol}')

if ncon != mxcon:
    # goto 10
    pass

last = 1
# goto 90


# 690
if iewald == 1:
    raise NotImplementedError()

i = mxcon // nvchg
k = i - nfl
j = 0
if i != 0:
    j = (100 * k) // i
print(f"Accepted Volume Moves = {k} Attempted = {i} ({j} %)")

if nflptr != 0:
    i = nflptr - nflprj
    j = (100 * i) // nflptr
print(f"Accepted Flip Moves = {i} Attempted = {nflptr} ({j} %)")

if ntdih != 0:
    #call plthis
    pass






