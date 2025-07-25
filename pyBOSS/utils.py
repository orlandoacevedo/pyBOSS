import os
import math
import numpy as np


def norcross(a,b):
    """calculate normalization for cross product a and b, 3x1"""
    x = a[1]*b[2] - b[1]*a[2]
    y = b[0]*a[2] - a[0]*b[2]
    z = a[0]*b[1] - b[0]*a[1]
    p = x*x + y*y + z*z
    s = 1.0/math.sqrt(p) if p > 0.0 else 0.0
    return [x*s, y*s, z*s]


def hessian(a,b,c):
    """calculate internal Hessian matrix, 3x3"""
    vca = [a[i]-c[i] for i in range(3)]
    vcb = [b[i]-c[i] for i in range(3)]
    y = norcross(vca,vcb)
    x = norcross(vcb,y)
    z = norcross(x,y)
    return [x,y,z]


def genatom(a,b,c,r,angle,dihedral):
    """a is the reference, b and c are anchors, angle/dihedral should be radian"""
    st = math.sin(angle)
    x = r * st * math.cos(dihedral)
    y = r * st * math.sin(dihedral)
    z = -r * math.cos(angle)
    h = hessian(a,b,c)
    xx = x*h[0][0] + y*h[1][0] + z*h[2][0] + c[0]
    yy = x*h[0][1] + y*h[1][1] + z*h[2][1] + c[1]
    zz = x*h[0][2] + y*h[1][2] + z*h[2][2] + c[2]
    return [xx,yy,zz]


def zmat2cor(zmat,radian=None):
    """convert zmatrix to coordinates
    
    Rule:
        firat atom will be (0.0, 0.0, 0.0)
        second atom will at positive x-axis
        third atom will in positive xy plane
    
    Return:
        cor: List[(x,y,z),   ]
    """
    if len(zmat) == 0: return []
    if len(zmat) == 1: return [0.0, 0.0, 0.0]

    if radian is not True:
        # avoid change the initial values
        nzmat = [[j for j in i] for i in zmat]
        for v in nzmat[1:]:
            v[3] = v[3] * math.pi / 180.0
            v[5] = v[5] * math.pi / 180.0
        zmat = nzmat

    cor = [[0.0,0.0,0.0], [zmat[1][1],0.0,0.0], ]
    if len(zmat) == 2: return cor

    y = zmat[2][1] * math.sin(zmat[2][3])
    x = zmat[2][1] * math.cos(zmat[2][3])
    if zmat[2][0] == 2:
        x = cor[1][0] - x
    cor.append([x,y,0.0])

    for v in zmat[3:]:
        xyz = genatom(cor[v[4]-1],cor[v[2]-1],cor[v[0]-1],v[1],v[3],v[5])
        cor.append(xyz)
    return cor



def find_all_rings_dfs(nbor):
    visited = [False for i in range(len(nbor))]
    dfs = []
    rings = []
    while False in visited:
        t = visited.index(False)
        visited[t] = True
        for v in nbor[t]:
            dfs.append([t,v])
        while len(dfs):
            g = dfs.pop()
            t = g[-1]
            if visited[t]:
                for v in nbor[t]:
                    if v == g[0]:
                        if len(g) > 2:
                            rings.append(g)
            else:
                visited[t] = True
                for v in nbor[t]:
                    if v == g[0]:
                        if len(g) > 2:
                            rings.append(g)
                    else:
                        for i in range(len(dfs)):
                            if dfs[i][-1] == t:
                                dfs[i].append(v)
                        for i in range(len(g)):
                            dfs.append([*g[i:],v])
    return rings


class ClosedForm:
    """Calculate RMSD for obmols by using closedForm algorithm

    Args:
        vl : 2D n*3f : List[[float, float, float]] : Left input matrix
        vr : 2D n*3f : List[[float, float, float]] : Right input matrix
        Left matrix is used for reference
            -> new_vr_fit_on_vl = BestFit * vr

    Methods:
        centroid
        calc_N
        calc_left_rotM
        calc_bestfit

    Reference:
        Closed-form solution of absolute orientation using unit quaternions
        Berthold K. P. Horn
        J. Opt. Soc. Am. A 4, 629-642 (1987)
    """
    def __init__(self,vl=None,vr=None,*args,**kws):
        self.vl = vl
        self.vr = vr

    def calc_bestfit(self,vl=None,vr=None,centroid=None,centroid_vl=True,centroid_vr=True):
        """calc best fit structure
        both vl(reference) and vr(target) can be reused
        centroid: whether centroid results: None(Yes), False(on vl), True(on vr)
        """
        if vl is None: vl = self.vl
        if vr is None: vr = self.vr
        if centroid_vl:
            cvl,tl = self.centroid(vl)
        else:
            cvl = vl
            tl = [0.0, 0.0, 0.0]
        if centroid_vr:
            cvr,tr = self.centroid(vr)
        else:
            cvr = vr
            tr = [0.0, 0.0, 0.0]
        M = self.calc_left_rotM(cvl,cvr)
        fit = []
        for v in cvr:
            rx = v[0]*M[0][0] + v[1]*M[1][0] + v[2]*M[2][0]
            ry = v[0]*M[0][1] + v[1]*M[1][1] + v[2]*M[2][1]
            rz = v[0]*M[0][2] + v[1]*M[1][2] + v[2]*M[2][2]
            fit.append([rx,ry,rz])
        if centroid is True:
            tk = tr
        elif centroid is False:
            tk = tl
        else:
            tk = [0.0, 0.0, 0.0]
        for i in range(len(fit)):
            fit[i][0] += tk[0]
            fit[i][1] += tk[1]
            fit[i][2] += tk[2]
        return fit

    def centroid(self,v):
        """calc centroid vector and translation vector"""
        x = [i[0] for i in v]
        y = [i[1] for i in v]
        z = [i[2] for i in v]
        t = len(v)
        ax = sum(x) / t
        ay = sum(y) / t
        az = sum(z) / t
        # care, do not use in-place operation
        cv = [[0.0 for i in range(3)] for j in range(len(v))]
        for i in range(t):
            cv[i][0] = v[i][0] - ax
            cv[i][1] = v[i][1] - ay
            cv[i][2] = v[i][2] - az
        return cv,(ax,ay,az)

    def calc_N(self,vl,vr):
        """calc 4x4 real symmetric N matrix"""
        XxYx = 0.0
        XxYy = 0.0
        XxYz = 0.0
        XyYx = 0.0
        XyYy = 0.0
        XyYz = 0.0
        XzYx = 0.0
        XzYy = 0.0
        XzYz = 0.0
        # careful of the sequence: X-r, Y-l
        # for rotation l = Rr
        for i,p in enumerate(vl):
            XxYx += p[0] * vr[i][0]
            XxYy += p[0] * vr[i][1]
            XxYz += p[0] * vr[i][2]

            XyYx += p[1] * vr[i][0]
            XyYy += p[1] * vr[i][1]
            XyYz += p[1] * vr[i][2]

            XzYx += p[2] * vr[i][0]
            XzYy += p[2] * vr[i][1]
            XzYz += p[2] * vr[i][2]

        N = [[0.0, 0.0, 0.0, 0.0] for i in range(4)]

        N[0][0] = XxYx + XyYy + XzYz
        N[0][1] = XyYz - XzYy
        N[0][2] = XzYx - XxYz
        N[0][3] = XxYy - XyYx

        N[1][0] = N[0][1]
        N[1][1] = XxYx - XyYy - XzYz
        N[1][2] = XxYy + XyYx
        N[1][3] = XzYx + XxYz

        N[2][0] = N[0][2]
        N[2][1] = N[1][2]
        N[2][2] = -XxYx + XyYy - XzYz
        N[2][3] = XyYz + XzYy

        N[3][0] = N[0][3]
        N[3][1] = N[1][3]
        N[3][2] = N[2][3]
        N[3][3] = -XxYx - XyYy + XzYz

        return N

    def calc_left_rotM(self,vl,vr):
        """calc rotation matrix for vr*M,
        quaternion is got from the vector which is corresponding to
        largest positive eigenvalue

        M  : 2D 3*3f : rotation matrix for vr, vr*M, in element-wise operation
        """
        N = self.calc_N(vl,vr)
        values,vectors = np.linalg.eig(N)
        ndx = np.where(values == max(values))
        ndx = ndx[0][0]
        # For numpy, eigenvectors are correspondingly put in column
        # note, this vector has already been normalized
        V = vectors[:,ndx]
        M = [[0.0, 0.0, 0.0] for i in range(3)]

        M[0][0] = 1 - 2 * (V[2]*V[2] + V[3]*V[3])
        M[0][1] = 2 * (V[1]*V[2] - V[3]*V[0])
        M[0][2] = 2 * (V[1]*V[3] + V[2]*V[0])

        M[1][0] = 2 * (V[1]*V[2] + V[3]*V[0])
        M[1][1] = 1 - 2 * (V[1]*V[1] + V[3]*V[3])
        M[1][2] = 2 * (V[2]*V[3] - V[1]*V[0])

        M[2][0] = 2 * (V[1]*V[3] - V[2]*V[0])
        M[2][1] = 2 * (V[2]*V[3] + V[1]*V[0])
        M[2][2] = 1 - 2 * (V[1]*V[1] + V[2]*V[2])

        return M


def align_onto_z_axis(cor,sol1=None,sol2=None):
    n = 0
    rmax = 0.0
    for i,l in enumerate(cor):
        r = l[0]*l[0] + l[1]*l[1] + l[2]*l[2]
        if r > rmax:
            rmax = r
            n = i
    v = cor[n]

    # first, rotate y axis to yz-plane
    t = pow(v[0]*v[0]+v[2]*v[2],0.5)
    ct = v[2] / t
    a = math.acos(ct)
    st = math.sin(a)
    if v[0] < 0.0: st = -st
    ry = [                  # counter clockwise
        [ct,  0.0,   st],
        [0.0, 1.0,  0.0],
        [-st, 0.0,   ct]
    ]

    # then, rotate x axis to xz-plane
    t = pow(rmax,0.5)
    st = v[1] / t
    a = math.asin(st)
    ct = math.cos(a)
    if v[1] < 0.0:
        st = -st
        ct = -ct
    rx = [                  # counter clockwise
        [1.0,   0.0, 0.0],
        [0.0,   ct,  st],
        [0.0,   -st,   ct]
    ]

    rmat = [
        [sum([ry[i][k]*rx[k][j] for k in range(3)]) for j in range(3)] for i in range(3)
    ]

    fin = []
    for c in cor:
        t = [sum([c[j]*rmat[j][i] for j in range(3)]) for i in range(3)]
        fin.append(t)

    finsol1 = []
    if sol1:
        for c in sol1:
            t = [sum([c[j]*rmat[j][i] for j in range(3)]) for i in range(3)]
            finsol1.append(t)

    finsol2 = []
    if sol2:
        for c in sol2:
            t = [sum([c[j]*rmat[j][i] for j in range(3)]) for i in range(3)]
            finsol2.append(t)

    return fin, finsol1, finsol2


def randu(irn,ichg=0):
    icnt = 0
    ichg = 1167
    ichg0 = 1167
    imul = 1173
    icon = 458753759
    imod = 1048573
    icn0 = 458753759
    jmul = 1161
    jcon = 458716759
    jmod = 1048573
    jrn = 12469
    icnt += 1
    if icnt == ichg:
        jrn = jrn * jmul + jcon
        jrn = jrn % jmod
        rnj = float(jrn) / float(jmod)
        fac = 1.0+0.5*rnj if rnj > 0.5 else 1.0-0.5*rnj
        fac = float(icn0) * fac
        icon = int(fac)
        jrn = jrn * jmul + jcon
        jrn = jrn % jmod
        rnj = float(jrn) / float(jmod)
        fac = 1.0+0.5*rnj if rnj > 0.5 else 1.0-0.5*rnj
        fac = float(ichg0) * fac
        ichg = int(fac)
        icnt = 0
    irn = irn * imul + icon
    irn = irn % imod
    return irn, float(irn)/float(imod)


def randint(n1, n2, x=None):
    """number of `n2` replaced in `n1`, thus `n1>=n2`"""
    if n1 < n2: n1, n2 = n2, n1
    if not x: x = 786123
    nmr = [0 for i in range(n2)]
    for v in range(n2):
        while True:
            x, y = randu(x)
            t = int(n1*y) + 1
            bo = True
            for i in range(v):
                if nmr[i] == t:
                    bo = False
                    break
            if bo:
                nmr[v] = t
                break
    return x,[i-1 for i in sorted(nmr)]



