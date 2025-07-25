"""
notes:

    BOSS: AW, BW, QW: format: [solvent_atoms, solvent_type]
    This:   [solvent_type, solvent_atoms]

    BOSS: ityp, ityp1, ityp2,   format: [sol1-aidx1, sol2-aidx1, sol3-aidx1, sol1-aidx2, sol2-aidx2, sol3-aidx2, ...]
      ->  a, b, q:              format: [sol1      , sol2      , sol3,      ...]
    This: a, b, q:  [[sol1], [sol2], [sol3]]

    BOSS: nstyp[i] starts at 1
"""


def hbond(nsatm, iztyp, asol, iatno, nmol, nsvat, nstyp, isatno, ac, aw, a, ityp, **kws):
    rhb2 = 6.25
    naccp = 0
    ndon = 0
    for i in range(nsatm):
        if iztyp[i] <= 0:
            continue
        j = iatno[i]
        x = asol[i][0]
        y = asol[i][1]
        z = asol[i][2]
        bo = False
        if j == 1:
            #if a[ityp[i]] == 0:
            if a[0][i] == 0:
                bo = True
        if j in [7, 8]:
            bo = True
        if not (bo or j==16):
            continue
        for l in range(nmol):
            ltyp = nstyp[l]
            for k in range(nsvat[ltyp]):
                m = isatno[k][ltyp]
                if j == 1:
                    if m in [7, 8, 16]:
                        r = (x - ac[l][k][0])**2 + (y - ac[l][k][1])**2 + (z - ac[l][k][2])**2
                        if r <= rhb2:
                            ndon += 1
                elif m == 1:
                    if aw[ltyp][k] != 0:
                        continue
                    r = (x - ac[l][k][0])**2 + (y - ac[l][k][1])**2 + (z - ac[l][k][2])**2
                    if r <= rhb2:
                        naccp += 1
    return naccp, ndon


def ssljco(nmol, noss, nstyp, modsv1, modsv2, anew, nsvat, ac, **kws):
    esljol = 0.0
    escool = 0.0
    sslj = 0.0
    if nmol == 0:
        return esljol, sslj, escool
    if noss == 1:
        return esljol, sslj, escool
    movtyp = 0
    n = nmol - 1
    for i in range(n):
        k = i + 1
        nmov = i
        jtyp = nstyp[i]
        natoms = nsvat[jtyp]
        ntyp = modsv1
        if jtyp == 2:
            ntyp = modsv2
        for l in range(natoms):
            for m in range(3):
                anew[l][m] = ac[i][l][m]
        e = 0.0
        for j in range(k,nmol):
            if nsvat[1] == 0:
                if modsv1 > 2:
                    e = sspot(j)
                else:
                    e = wwpot(j)
            else:
                jtyp = modsv1
                if nstyp[j] == 1:
                    jtyp = modsv2
                if ntyp <= 2 and jtyp <= 2:
                    e = wwpot(j)
                else:
                    e = sspot(j)
            esljol += sslj
            escool += e - sslj
    esljol += ecut()
    return esljol, sslj, escool


def sspot(
    nm, edge, islab, nstyp, nsvat, ac, nmov, ncutas, icutas, anew, edg2, rcutsq, nmol2,
    qw, aw, bw, iewald, irfon, rffac, natmx, rl2, rul, ru2,
    isvaty, #aas, bbs, qqs,         # only for solvent-solvent interactions
                                    # i-j valid when isvaty[i][j] != 0
    **kws
):
    sspot = 0.0
    sslj = 0.0
    xedg = edge[0]
    yedg = edge[1]
    zedg = edge[2]
    
    if islab == 1:
        zedg = 0.0

    n1 = nstyp[nmov]
    n2 = nstyp[nm]
    nat2 = nsvat[n2]
    
    c = [[0.0,0.0,0.0] for i in range(nat2)]
    for i in range(3):
        for j in range(nat2):
            c[j][i] = ac[nm][j][i]

    wik = 100000000
    for i in range(ncutas[n1]):
        nw = icutas[n1][i]
        for j in range(ncutas[n2]):
            nc = icutas[n2][j]
            xiz = 0.0
            zz = abs(anew[nw][2] - c[nc][2])
            if zz > edg2[2]:
                xiz = zedg
            twik = (zz - xiz) ** 2
            if twik > rcutsq:
                continue
            xix = 0.0
            xx = abs(anew[nw][0] - c[nc][0])
            if xx > edg2[0]:
                xix = xedg
            twik += (xx - xix) ** 2
            if twik > rcutsq:
                continue
            xiy = 0.0
            yy = abs(anew[nw][1] - c[nc][1])
            if yy > edg2[1]:
                xiy = yedg
            twik += (yy - xiy) ** 2
            wik = min(wik, twik)

    if wik > rcutsq:
        return sspot

    if nmol2 != 0:  # done
        for i in range(nsvat[n1]):
            if isvaty[n1][i] <= 0:      # change `isvaty[i][n1]` to `isvaty[n1][i]`
                continue

            x = anew[i][0]
            y = anew[i][1]
            z = anew[i][2]
            qi = qw[n1][i]
            ai = aw[n1][i]
            iflag = 1 if ai == 0.0 else 0
            bi = bw[n1][i]

            for j in range(nat2):
                if isvaty[n2][j] <= 0:
                    continue
                
                xx = abs(x - c[j][0])
                yy = abs(y - c[j][1])
                zz = abs(z - c[j][2])
                xix = xedg if xx > edg2[0] else 0.0
                xiy = yedg if yy > edg2[1] else 0.0
                xiz = zedg if zz > edg2[2] else 0.0
                
                rr = (xx - xix) ** 2 + (yy - xiy) ** 2 + (zz - xiz) ** 2
                if iflag != 1:
                    r6 = 1.0 / (rr ** 3)
                    aaij = ai * aw[n2][j]
                    bbij = bi * bw[n2][j]
                    xlj = (aaij*r6 - bbij) * r6
                    sspot += xlj
                    sslj += xlj
                
                qqij = qi * qw[n2][j]
                if qqij != 0.0:
                    if iewald == 1:
                        raise NotImplementedError('Ewald not supported')
                    else:
                        r1 = 1.0 / pow(rr,0.5)
                        if irfon == 0:
                            sspot += qqij * r1
                        else:
                            sspot += qqij * (r1 + rr * rffac)
    else:
        # future, sovlent self interaction
        k = 0
        for i in range(natmx):
            if isvaty[0][i] <= 0:
                continue

            x = anew[i][0]
            y = anew[i][1]
            z = anew[i][2]

            for j in range(natmx):
                if isvaty[0][j] <= 0:
                    continue
                
                xx = abs(x - c[j][0])
                yy = abs(y - c[j][1])
                zz = abs(z - c[j][2])
                xix = xedg if xx > edg2[0] else 0.0
                xiy = yedg if yy > edg2[1] else 0.0
                xiz = zedg if zz > edg2[2] else 0.0
                
                rr = (xx - xix) ** 2 + (yy - xiy) ** 2 + (zz - xiz) ** 2
                if aas[k] != 0.0:
                    r6 = 1.0 / (rr ** 3)
                    xlj = (aas[k] * r6 - bbs[k]) * r6
                    sspot += xlj
                    sslj += xlj

                if qqs[k] != 0.0:
                    if iewald == 1:
                        raise NotImplementedError('Ewald not supported')
                    else:
                        r1 = 1.0 / pow(rr,0.5)
                        if irfon == 0:
                            sspot += qqs[k] * r1
                        else:
                            sspot += qqs[k] * (r1 + rr*rffac)
                k += 1

    if wik < rl2:
        return sspot
    
    # not done
    feat = rul * (ru2 - wik)
    sspot *= feat
    sslj *= feat

    return sspot


def wwpot(
    nm, anew, ac, islab, edg2, edge, rcutsq, aoo, coo, natoms, irfon, rffac, qq, rl2, ru2, rul,
    kq, iewald, **kws
):      # done
    xiz = 0.0
    zz = abs(anew[0][2] - ac[nm][0][2])
    if islab != 1:
        if zz > edg2[2]:
            xiz = edge[2]
    wik = (zz - xiz) ** 2

    ckws = {'sslj':0.0, 'ssall':0.0}

    if wik > rcutsq:
        return ckws

    xix = 0.0
    xx = abs(anew[0][0] - ac[nm][0][0])
    if xx > edg2[0]:
        xix = edge[0]
    wik += (xx - xix) ** 2

    if wik > rcutsq:
        return ckws

    xiy = 0.0
    yy = abs(anew[0][1] - ac[nm][0][1])
    if yy > edg2[1]:
        xiy = edge[1]
    wik += (yy - xiy) ** 2

    if wik > rcutsq:
        return ckws

    r6 = 1.0 / (wik * wik * wik)
    val = r6 * (aoo * r6 - coo)
    ckws['sslj'] = val

    anm = [[0.0,0.0,0.0] for i in range(natoms)]
    for i in range(3):
        for j in range(natoms):
            anm[j][i] = ac[nm][j][i]

    knt = 0
    for i in range(kq, natoms):
        for j in range(kq, natoms):
            xx = abs(anew[i][0] - anm[j][0])
            yy = abs(anew[i][1] - anm[j][1])
            zz = abs(anew[i][2] - anm[j][2])
            rr = (xx - xix) ** 2 + (yy - xiy) ** 2 + (zz - xiz) ** 2
            if iewald == 1:
                raise NotImplementedError('Ewald not supported')
            else:
                r1 = 1.0 / pow(rr,0.5)
                if irfon == 0:
                    val += qq[knt] * r1
                else:
                    val += qq[knt] * (r1 + rr * rffac)
            knt += 1
    if wik < rl2:
        ckws['ssall'] = val
        return ckws

    val *= rul * (ru2 - wik)
    ckws['ssall'] = val
    return ckws


def ecut(natmx, nocut, noss, rcut, nmol, nsvat, modsv1, nmol2, modsv2, aw, bw, vnew, twopi, **kws):
    val = 0.0

    if natmx == 0 or nocut == 1 or noss == 1:
        return val

    rc3 = 1.0 / (rcut ** 3)
    rc9 = rc3 ** 3
    xn = float(nmol)
    xm = 0.0

    if nsvat[1] == 0:
        if modsv1 <= 2 or modsv1 == 13:
            return val
    else:
        xn = float(nmol - nmol2)
        xm = float(nmol2)
        if modsv1 <= 2:
            xn = 0.0
        if modsv2 <= 2:
            xm = 0.0

    for i in range(natmx):
        ai = aw[0][i]
        bi = bw[0][i]
        aj = aw[1][i]
        bj = bw[1][i]

        for j in range(natmx):
            val += (twopi * xn * xn / (3.0 * vnew)) * (ai * aw[0][j] * rc9 / 3.0 - bi * bw[0][j] * rc3)
            val += (twopi * xm * xm / (3.0 * vnew)) * (aj * aw[1][j] * rc9 / 3.0 - bj * bw[1][j] * rc3)
            val += (twopi * xn * xm / (3.0 * vnew)) * (ai * aw[1][j] * rc9 / 3.0 - bi * bw[1][j] * rc3) * 2.0

    return val


def isolut(iatom, nsolut, nlatm):
    for i in range(nsolut):
        if iatom <= nlatm[i]: break
    return i


def cappot(r2, capfk, capsq5, caprad):
    if r2 <= capsq5: return 0.0
    val = capfk * (pow(r2,0.5) - caprad - 5) ** 2
    return val


def sxpot(
    nm, movtyp, nstyp, nsvat, ac, nmov, nsatm, icut, ncutas, icutas, anew, ncent1, ncent2,
    asol, islab, edg2, edge, nsolut, icutat, iztyp, nrdfs, scutsq, sl2, su2, sul, nsatno, nofep,
    a, b, q, ityp, ityp1, ityp2,  anew1, anew2, aw, bw, qw, irfon, iewald, rffac, isolec,
    nrdfa1, nrdfa2, asol1, asol2, capfk, icapat, **kws
):
    if movtyp != 2:
        n1 = nstyp[nm]
        natoms = nsvat[n1]
        anm = [[0.0,0.0,0.0] for i in range(natoms)]
        for i in range(3):
            for j in range(natoms):
                anm[i][j] = ac[nm][j][i]
    else:
        n1 = nstyp[nmov]
        natoms = nsvat[n1]

    iyes = [1 for i in range(nsatm)]
    if icut < 2:
        wik = 1e20
        for j in range(ncutas[n1]):
            twik = 0.0
            nc = icutas[n1][j]
            for i in range(3):
                cnew = anew[nc][i]
                cnm = anm[i][nc]
                if movtyp == 1:
                    cnew = 0.5 * (anew[ncent1][i] + anew[ncent2][i])
                elif movtyp == 2:
                    cnm = 0.5 * (asol[ncent1][i] + asol[ncent2][i])
                x = abs(cnew - cnm)
                if not(islab == 1 and i == 2):
                    if x > edg2[i]:
                        x -= edge[i]
                twik += x**2
            wik = min(wik, twik)
    else:
        rmin = [1e20 for i in range(nsolut)]
        wik = 1e20
        for i in range(nsatm):
            bo = True
            if icut >= 3:
                if i in icutat:
                    bo = False
            else:
                if iztyp[i] == -1:
                    continue
            if bo: continue

            m = isolut(i)
            for l in range(ncutas[n1]):
                nc = icutas[n1][l]
                x = 0.0
                for j in range(3):
                    if movtyp != 2:
                        y = abs(anew[i][j] - anm[j][nc])
                    else:
                        y = abs(anew[nc][j] - asol[i][j])
                    if not(islab == 1 and j == 2):
                        if y > edg2[j]:
                            y -= edge[j]
                    x += y**2
                rmin[m] = min(rmin[m], x)
                wik = min(wik, x)

    rnew = [100.0 for i in range(nrdfs)]
    elj = [0.0, 0.0, 0.0]
    ecoul = [0.0, 0.0, 0.0]
    ec1 = ec2 = ec3 = 0.0
    ecsx = elsx = 0.0

    if wik > scutsq:
        return 0.0

    scal = 1.0
    if wik >= sl2:
        scal = sul * (su2 - wik)

    if icut != 0:
        j = 0
        for i in range(nsolut):
            if rmin[i] > scutsq:
                for k in range(j, j+nsatno[i]):
                    iyes[k] = 0
            j += nsatno[i]

    k0 = 0
    xedg, yedg, zedg = edge[0], edge[1], edge[2]
    if islab == 1:
        zedg = 0.0

    llim = 3
    if nofep == 1:
        llim = 1

    if movtyp != 2:
        for i in range(nsatm):
            if iztyp[i] == -1 or iyes[i] == 0:
                continue
            
            m = isolut(i)

            for l in range(1, llim + 1):
                xix = xiy = xiz = 0.0
                if l == 1:
                    #ii = ityp[i]
                    xx, yy, zz = anew[i][0], anew[i][1], anew[i][2]
                elif l == 2:
                    #ii = ityp1[i]
                    xx, yy, zz = anew1[i][0], anew1[i][1], anew1[i][2]
                elif l == 3:
                    #ii = ityp2[i]
                    xx, yy, zz = anew2[i][0], anew2[i][1], anew2[i][2]

                #aii, bii, qii = a[ii], b[ii], q[ii]
                aii, bii, qii = a[l-1][i], b[l-1][i], q[l-1][i]

                for j in range(natoms):
                    x = abs(xx - anm[0][j])
                    y = abs(yy - anm[1][j])
                    z = abs(zz - anm[2][j])

                    if j == 0:
                        if x > edg2[0]:
                            xix = xedg
                        if y > edg2[1]:
                            xiy = yedg
                        if z > edg2[2]:
                            xiz = zedg

                    rr = (x - xix)**2 + (y - xiy)**2 + (z - xiz)**2
                    r1 = pow(rr,0.5)
                    r6 = 1.0 / (rr * rr * rr)
                    el = scal * r6 * (aii * aw[n1][j] * r6 - bii * bw[n1][j])
                    elj[l-1] += el
                    ec = scal * qii * qw[n1][j]

                    if irfon == 0:
                        if iewald == 1:
                            raise NotImplementedError('Ewald not supported')
                        else:
                            ec /= r1
                    else:
                        ec *= (1.0 / r1 + rr * rffac)
                    ecoul[l-1] += ec

                    if l == 1 and m == isolec:
                        ecsx += ec
                        elsx += el

                    for kk in range(k0, nrdfs):
                        if i == nrdfa1[kk]:
                            if j == nrdfa2[kk]:
                                rnew[kk] = r1
                                k0 = kk + 1
                                break
    else:
        for i in range(nsatm):
            if iztyp[i] == -1 or iyes[i] == 0:
                continue

            m = isolut(i)

            for l in range(1, llim + 1):
                xix = xiy = xiz = 0.0
                if l == 1:
                    ii = ityp[i]
                    xx = asol[i][0]
                    yy = asol[i][1]
                    zz = asol[i][2]
                elif l == 2:
                    ii = ityp1[i]
                    xx = asol1[i][0]
                    yy = asol1[i][1]
                    zz = asol1[i][2]
                elif l == 3:
                    ii = ityp2[i]
                    xx = asol2[i][0]
                    yy = asol2[i][1]
                    zz = asol2[i][2]

                #aii, bii, qii = a[ii], b[ii], q[ii]
                aii, bii, qii = a[l-1][i], b[l-1][i], q[l-1][i]


                for j in range(natoms):
                    x = abs(xx - anew[j][0])
                    y = abs(yy - anew[j][1])
                    z = abs(zz - anew[j][2])

                    if j == 1:
                        if x > edg2[0]:
                            xix = xedg
                        if y > edg2[1]:
                            xiy = yedg
                        if z > edg2[2]:
                            xiz = zedg

                    rr = (x - xix) ** 2 + (y - xiy) ** 2 + (z - xiz) ** 2
                    r1 = pow(rr,0.5)
                    r6 = 1.0 / (rr * rr * rr)

                    el = scal * r6 * (aii * aw[n1][j] * r6 - bii * bw[n1][j])
                    elj[l] += el
                    ec = scal * qii * qw[n1][j]

                    if irfon == 0:
                        if iewald == 1:
                            raise NotImplementedError('Ewald not supported')
                        else:
                            ec /= r1
                    else:
                        ec *= (1.0 / r1 + rr * rffac)

                    ecoul[l] += ec

                    if l == 1:
                        if m == isolec:
                            ecsx += ec
                            elsx += el

                        for kk in range(k0, nrdfs + 1):
                            if i == nrdfa1[kk]:
                                if j == nrdfa2[kk]:
                                    rnew[kk] = r1
                                    k0 = kk + 1
                                    break

    if capfk == 0.0:
        return

    if icapat != 0:
        x = y = z = 0.0
        if movtyp == 1:
            for i in range(3):
                xx = ac[nm][0][i]
                x += (anew[icapat][i] - xx) ** 2
                y += (anew1[icapat][i] - xx) ** 2
                z += (anew2[icapat][i] - xx) ** 2
        else:
            for i in range(3):
                xx = anew[0][i]
                x += (asol[icapat][i] - xx) ** 2
                y += (asol1[icapat][i] - xx) ** 2
                z += (asol2[icapat][i] - xx) ** 2

        ec1 = cappot(x)

        if nofep == 0:
            ec2 = cappot(y)
            ec3 = cappot(z)

    # Final calculations
    sxpot = elj[0] + ecoul[0] + ec1

    if nofep == 0:
        es1 = elj[1] + ecoul[1] + ec2
        es2 = elj[2] + ecoul[2] + ec3
    else:
        es1 = sxpot
        es2 = sxpot
    return


def wxpot(
    movtyp, nstyp, iztyp, ac, nmov, nsatm, edg2, ncutat, ityp, ityp1, ityp2,
    anew, anew1, anew2, a, b, q, aw, bw, qw, asol, asol1, asol2, edge,
    nm, nsvat, modsv1, modsv2, icut, icutat, ncent1, ncent2, islab, isolec, irfon, iewald,
    nrdfs,  nrdfa1, nrdfa2, scutsq, nsatno, sl2, sul, su2, nofep, rffac, nsolut, capfk, icapat, **kws
):
    if movtyp != 2:     # done
        n1 = nstyp[nm]
        natoms = nsvat[n1]
        anm = [[0.0 for i in range(natoms)] for j in range(3)]
        for i in range(3):
            for j in range(natoms):
                anm[i][j] = ac[nm][j][i]
    else:
        n1 = nstyp[nmov]
        natoms = nsvat[n1]

    modsv = modsv1
    if n1 == 1:     # check whether is second solvent
        modsv = modsv2

    iyes = [1 for i in range(nsatm)]
    if icut < 2:
        wik = 0.0
        for i in range(3):
            cnew = anew[0][i]
            cnm = anm[i][0]
            if movtyp == 1:
                cnew = 0.5 * (anew[ncent1][i] + anew[ncent2][i])
            elif movtyp == 2:
                cnm = 0.5 * (asol[ncent1][i] + asol[ncent2][i])

            x = abs(cnew - cnm)
            if not (islab == 1 and i == 2):
                if x > edg2[i]:
                    x -= edge[i]
            wik += x**2
    else:   # done
        rmin = [100000000 for i in range(nsolut)]
        wik = 100000000
        for i in range(nsatm):
            if icut >= 3:
                bo = True
                for k in range(ncutat):
                    if i == icutat[k]:      # skip dummy atom
                        bo = False
                        break
                if bo: continue
            else:
                if iztyp[i] == -1:
                    continue

            x = 0.0
            for j in range(3):
                if movtyp != 2:
                    y = abs(anew[i][j] - anm[j][0])
                else:
                    y = abs(anew[0][j] - asol[i][j])

                if not (islab == 1 and j == 2):
                    if y > edg2[j]:
                        y -= edge[j]
                x += y**2

            #j = isolut(i)
            j = kws['isolute'][i]

            rmin[j] = min(rmin[j], x)
            wik = min(wik, x)

    rnew = [100.0 for i in range(nrdfs)]

    ec1 = 0.0
    ec2 = 0.0
    ec3 = 0.0
    ecsx = 0.0
    elsx = 0.0

    elj = [0.0, 0.0, 0.0]
    ecoul = [0.0, 0.0, 0.0]
    gkws = {
        'elj': elj,
        'ecoul': ecoul,
        'ec1': ec1,
        'ec2': ec2,
        'ec3': ec3
    }
    gkws['wik'] = wik

    if wik > scutsq:
        return gkws

    scal = 1.0
    if wik >= sl2:
        scal = sul * (su2 - wik)

    if icut != 0:
        j = 0   # done
        for i in range(nsolut):
            if rmin[i] > scutsq:
                l = j + nsatno[i]
                for k in range(j, l):
                    iyes[k] = 0
            j += nsatno[i]

    k0 = 0
    xedg = edge[0]
    yedg = edge[1]
    zedg = edge[2]
    if islab == 1:
        zedg = 0.0

    llim = 3
    if nofep == 1:
        llim = 1

    if movtyp != 2:     # done
        for i in range(nsatm):
            if iztyp[i] == -1 or iyes[i] == 0:
                continue

            #m = isolut(i)
            m = kws['isolute'][i]

            for l in range(llim):
                xix = xiy = xiz = 0.0
                if l == 0:
                    #ii = ityp[i]
                    xx, yy, zz = anew[i][0], anew[i][1], anew[i][2]
                elif l == 1:
                    #ii = ityp1[i]
                    xx, yy, zz = anew1[i][0], anew1[i][1], anew1[i][2]
                else:
                    #ii = ityp2[i]
                    xx, yy, zz = anew2[i][0], anew2[i][1], anew2[i][2]

                #aii, bii, qii = a[ii], b[ii], q[ii]
                aii, bii, qii = a[l][i], b[l][i], q[l][i]

                for j in range(natoms):
                    x = abs(xx - anm[0][j])
                    y = abs(yy - anm[1][j])
                    z = abs(zz - anm[2][j])

                    bo = True
                    if j == 0:
                        if x > edg2[0]:
                            xix = xedg
                        if y > edg2[1]:
                            xiy = yedg
                        if z > edg2[2]:
                            xiz = zedg
                        rr = (x - xix)**2 + (y - xiy)**2 + (z - xiz)**2
                        r1 = pow(rr,0.5)
                        r6 = 1.0 / (rr**3)
                        xlj = scal * r6 * (aii * aw[n1][0] * r6 - bii * bw[n1][0])
                        elj[l] += xlj
                        if l == 0 and m == isolec:
                            elsx += xlj
                        if modsv == 2:
                            #goto 140
                            bo = False
                    else:
                        rr = (x - xix)**2 + (y - xiy)**2 + (z - xiz)**2
                        r1 = pow(rr,0.5)

                    if bo:
                        coul = qii * qw[n1][j] * scal
                        if irfon == 0:
                            if iewald == 1:
                                raise NotImplementedError('`iewald==1`')
                            else:
                                coul /= r1
                        else:
                            coul = coul * (1.0/r1 + rr*rffac)

                        ecoul[l] += coul
                        if l == 0 and m == isolec:
                            ecsx += coul

                        if j > 1:
                            continue

                    if l == 0:  # done
                        for kk in range(k0, nrdfs):
                            #if i == nrdfa1[kk] and j == nrdfa2[kk]:
                            if i+1 == nrdfa1[kk] and j+1 == nrdfa2[kk]:     # input index starts at 1
                                rnew[kk] = r1
                                k0 = kk + 1
    else:
        for i in range(nsatm):
            if iztyp[i] == -1 or iyes[i] == 0:
                continue

            #m = isolut(i)
            m = kws['isolute'][i]

            for l in range(1, llim+1):
                xix = 0.0
                xiy = 0.0
                xiz = 0.0

                if l == 1:
                    #ii = ityp[i]
                    xx = asol[i][0]
                    yy = asol[i][1]
                    zz = asol[i][2]
                elif l == 2:
                    #ii = ityp1[i]
                    xx = asol1[i][0]
                    yy = asol1[i][1]
                    zz = asol1[i][2]
                elif l == 3:
                    #ii = ityp2[i]
                    xx = asol2[i][0]
                    yy = asol2[i][1]
                    zz = asol2[i][2]

                #aii, bii, qii = a[ii], b[ii], q[ii]
                aii, bii, qii = a[l-1][i], b[l-1][i], q[l-1][i]

                for j in range(natoms):
                    x = abs(xx - anew[j][0])
                    y = abs(yy - anew[j][1])
                    z = abs(zz - anew[j][2])

                    bo = True
                    if j == 1:
                        if x > edg2[0]:
                            xix = xedg
                        if y > edg2[1]:
                            xiy = yedg
                        if z > edg2[2]:
                            xiz = zedg

                        rr = (x - xix)**2 + (y - xiy)**2 + (z - xiz)**2
                        r1 = rr**0.5
                        r6 = 1.0 / (rr**3)
                        xlj = scal * r6 * (aii * aw[n1][0] * r6 - bii * bw[n1][0])
                        elj[l-1] += xlj

                        if l == 1 and m == isolec:
                            elsx += xlj

                        if modsv == 2:
                            #goto 200
                            bo = False
                    else:
                        rr = (x - xix)**2 + (y - xiy)**2 + (z - xiz)**2
                        r1 = rr**0.5

                    if bo:
                        coul = qii * qw[n1][j] * scal

                        if irfon == 0:
                            if iewald == 1:
                                raise NotImplementedError('`iewald`')
                            else:
                                coul /= r1
                        else:
                            coul *= (1.0 / r1 + rr * rffac)

                        ecoul[l-1] += coul

                        if l == 1 and m == isolec:
                            ecsx += coul

                        if j >= 2:
                            continue

                    if l == 1:
                        for kk in range(k0, nrdfs):
                            if i == nrdfa1[kk] and j == nrdfa2[kk]:
                                rnew[kk] = r1
                                k0 = kk + 1
                                break

    if capfk == 0.0:
        return gkws
    
    # not done
    if icapat != 0:
        x = y = z = 0.0
        if movtyp == 1:
            for i in range(3):
                xx = ac[nm][0][i]
                x += (anew[icapat][i] - xx) ** 2
                y += (anew1[icapat][i] - xx) ** 2
                z += (anew2[icapat][i] - xx) ** 2
        else:
            for i in range(3):
                xx = anew[0][i]
                x += (asol[icapat][i] - xx) ** 2
                y += (asol1[icapat][i] - xx) ** 2
                z += (asol2[icapat][i] - xx) ** 2

        ec1 = cappot(x)

        if nofep == 0:
            ec2 = cappot(y)
            ec3 = cappot(z)

    wxpot = elj[0] + ecoul[0] + ec1

    if nofep == 0:
        es1 = elj[1] + ecoul[1] + ec2
        es2 = elj[2] + ecoul[2] + ec3
    else:
        es1 = wxpot
        es2 = wxpot
    
    return gkws



def xxpot(
    edge, islab, nbc, nbcat1, nbcat2, rbc, anew, bck0, bck1, bck2, nofep, anew1, anew2, 
    ntdih, phinew, dihv0, rc0, rc1, rc2, isqm, itydi, **kws
):
    exxnew = exxne1 = exxne2 = exxnec = exxnel = ebcnew = ebcne1 = ebcne2 = 0.0
    edgg = [edge[0], edge[1], edge[2]]
    
    if islab == 1:
        edgg[2] = 0.0
    
    if nbc != 0:
        for i in range(nbc):
            j = nbcat1[i]
            k = nbcat2[i]
            r = rbc[i]

            if k == 9999:
                rr1 = r * ((anew[j][0] - bck0[i])**2 + (anew[j][1] - bck1[i])**2 + (anew[j][2] - bck2[i])**2)
                ebcnew += rr1

                if nofep == 0:
                    rr2 = r * ((anew1[j][0] - bck0[i])**2 + (anew1[j][1] - bck1[i])**2 + (anew1[j][2] - bck2[i])**2)
                    rr3 = r * ((anew2[j][0] - bck0[i])**2 + (anew2[j][1] - bck1[i])**2 + (anew2[j][2] - bck2[i])**2)
                    ebcne1 += rr2
                    ebcne2 += rr3
            
            elif j == 9999:
                j1 = int(bck0[i])
                j2 = int(bck1[i])
                j3 = int(bck2[i])
                
                if j1 == 0:
                    continue
                
                x = anew[j1][0]
                y = anew[j1][1]
                z = anew[j1][2]
                fac = 1.0
                
                if j2 != 0:
                    x += anew[j2][0]
                    y += anew[j2][1]
                    z += anew[j2][2]
                    fac = 2.0
                
                if j3 != 0:
                    x += anew[j3][0]
                    y += anew[j3][1]
                    z += anew[j3][2]
                    fac += 1.0
                
                x /= fac
                y /= fac
                z /= fac
                
                rr1 = r * ((anew[k][0] - x)**2 + (anew[k][1] - y)**2 + (anew[k][2] - z)**2)
                ebcnew += rr1
                ebcne1 += rr1
                ebcne2 += rr1
            
            elif r < 0.0:
                rr1 = ((anew[j][0] - anew[k][0])**2 + (anew[j][1] - anew[k][1])**2 + (anew[j][2] - anew[k][2])**2)
                rr2 = bck1[i]**2
                
                if rr1 > rr2:
                    rr3 = bck2[i] * (rr1**0.5 - bck1[i])**2
                    ebcnew += rr3
                else:
                    rr2 = bck0[i]**2
                    if rr1 > rr2:
                        continue
                    rr3 = bck2[i] * (rr1**0.5 - bck0[i])**2

                ebcnew += rr3
            else:
                rr1 = ((anew[j][0] - anew[k][0])**2 + (anew[j][1] - anew[k][1])**2 + (anew[j][2] - anew[k][2])**2)
                ebcnew += bck0[i] * (rr1**0.5 - r)**2

                if nofep == 0:
                    rr2 = ((anew1[j][0] - anew1[k][0])**2 + (anew1[j][1] - anew1[k][1])**2 + (anew1[j][2] - anew1[k][2])**2)
                    rr3 = ((anew2[j][0] - anew2[k][0])**2 + (anew2[j][1] - anew2[k][1])**2 + (anew2[j][2] - anew2[k][2])**2)
                    ebcne1 += bck1[i] * (rr2**0.5 - r)**2
                    ebcne2 += bck2[i] * (rr3**0.5 - r)**2
    
    if nofep == 1:
        ebcne1 = ebcnew
        ebcne2 = ebcnew

    if ntdih != 0:
        for i in range(ntdih):
            if itydi[i] == 500:
                p = phinew[i]
                x = dihv0[i][0]
                e = dihv0[i][3] * (p - x) ** 2
                ebcnew += e
                if rc0 != rc1:
                    ebcne1 += e
                if rc0 != rc2:
                    ebcne2 += e

    if isqm == 1:
        return 
    raise NotImplementedError('isqm!=1')


def ranu(irn):
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
    jrn = 124690

    icnt += 1
    if icnt == ichg:
        jrn = (jrn * jmul + jcon) % jmod
        rnj = float(jrn) / float(jmod)
        if rnj > 0.5:
            fac = 1.0 + 0.5 * rnj
        else:
            fac = 1.0 - 0.5 * rnj
        fac = float(icn0) * fac
        icon = int(fac)

        jrn = (jrn * jmul + jcon) % jmod
        rnj = float(jrn) / float(jmod)

        if rnj > 0.5:
            fac = 1.0 + 0.5 * rnj
        else:
            fac = 1.0 - 0.5 * rnj
        fac = float(ichg0) * fac
        ichg = int(fac)
        icnt = 0

    irn = (irn * imul + icon) % imod
    return float(irn) / float(imod)


import math
import random  # Used for generating random numbers

def rotate(
    angle,
    movtyp, nstyp, nat, nmov, modcus, ncents, iflxsl, isol, nrota1, nrota2, nfatm, nvdih, nvbnd, nvang,
    nlatm, anew, anew1, anew2, **kws
):
    iax = int(3 * ranu() + 1)
    if movtyp != 1:
        iat = 1
        ifat = 2
        ilat = nat
        if nstyp[nmov] == modcus:
            iat = ncents
            ifat = 1
            if iflxsl == 1:
                ilat = 3
    else:
        if isol == 1:
            iat = nrota1
        elif isol == 2:
            iat = nrota2
        elif isol > 2:
            iat = nfatm[isol]

        ifat = nfatm[isol]
        ilat = nlatm[isol]
        n = nvdih + nvbnd + nvang
        if n != 0:
            if (ilat - ifat) > 2:
                ilat = ifat + 2

    csa = math.cos(angle)
    sna = math.sin(angle)

    if iax == 1:
        jax, kax = 1, 2
    elif iax == 2:
        jax, kax = 2, 0
    elif iax == 3:
        jax, kax = 0, 1

    x = anew[iat][jax]
    y = anew[iat][kax]

    for i in range(ifat, ilat):
        if i != iat:
            xo = anew[i][jax] - x
            yo = anew[i][kax] - y
            x1 = xo * csa - yo * sna
            y1 = xo * sna + yo * csa
            anew[i][jax] = x + x1
            anew[i][kax] = y + y1

    if movtyp != 0:
        x = anew1[iat][jax]
        y = anew1[iat][kax]
        for i in range(ifat, ilat):
            if i != iat:
                xo = anew1[i][jax] - x
                yo = anew1[i][kax] - y
                x1 = xo * csa - yo * sna
                y1 = xo * sna + yo * csa
                anew1[i][jax] = x + x1
                anew1[i][kax] = y + y1

        x = anew2[iat][jax]
        y = anew2[iat][kax]
        for i in range(ifat, ilat):
            if i != iat:
                xo = anew2[i][jax] - x
                yo = anew2[i][kax] - y
                x1 = xo * csa - yo * sna
                y1 = xo * sna + yo * csa
                anew2[i][jax] = x + x1
                anew2[i][kax] = y + y1


def movsvn(
    nsvat, nstyp, nmov, modcus, ncents, rdel, rdl2, anew, ac, islab, edg2, edge, adel, adl2, iflxsl,
    modsv1, modsv2, noss, nmol, nbuse, nbor, nrdfs, knew, rinc, rnew, rdlmin, eij,
    esonol, ess, esonco, esonlo, ecsx, essc, essl, elsx, es1, es2, ess1, ess2,
    esol1, esol2, eold, npols, icussl, iewald, nopref, wnew, wi, wisum, wik, wkc, **kws
):
    natoms = nsvat[nstyp[nmov]]
    nat = natoms
    nc = 1
    if nstyp[nmov] == modcus:
        nc = ncents

    for i in range(3):
        rr = -rdel + ranu() * rdl2
        for j in range(natoms):
            anew[j][i] = ac[nmov][j][i] + rr

        if islab == 1 and i == 2:
            continue

        if anew[nc][i] <= edg2[i]:
            if anew[nc][i] < -edg2[i]:
                for j in range(natoms):
                    anew[j][i] += edge[i]
        else:
            for j in range(natoms):
                anew[j][i] -= edge[i]

    achg = -adel + ranu() * adl2

    rotate(achg)

    if iflxsl == 1:
        if modcus == nstyp[nmov - 1]:
            raise NotImplementedError('iflxsl == 1')

    e2 = 0.0
    emov = [1 for i in range(2000)]

    ntyp = modsv1
    if nstyp[nmov - 1] == 2:
        ntyp = modsv2

    if noss != 1:
        for i in range(nmol):
            if nbuse == 1:
                j = nbor[nmov][i]
                if j == 0:
                    break
            else:
                if i == nmov:
                    continue
                j = i

            if nsvat[1] == 0:
                if modsv1 > 2:
                    e = sspot(j)
                else:
                    e = wwpot(j)
            else:
                jtyp = modsv1
                if nstyp[j] == 1:
                    jtyp = modsv2

                if ntyp <= 2 and jtyp <= 2:
                    e = wwpot(j)
                else:
                    e = sspot(j)

            emov[j] = e
            e2 += e

    movtyp = 2

    if noss == 0:
        if ntyp > 2:
            esol = sxpot(0)
        else:
            esol = wxpot(0)
    else:
        raise NotImplementedError('`sxnoss`')

    movtyp = 0

    if nrdfs != 0:
        for i in range(nrdfs):
            knew[i] = int(rinc * (rnew[i] - rdlmin)) + 1

    e1 = 0.0

    for i in range(nmov-1):
        e1 += eij[nmol*i + nmov - i - (i*(i-1)) // 2]

    j = nmol*nmov - nmov - (nmov * (nmov-1)) // 2
    for i in range(nmov, nmol):
        e1 += eij[j + i + 1]

    emov[nmov] = 0.0
    esone = esonol + esol - ess[nmov]
    esonc = esonco + ecsx - essc[nmov]
    esonl = esonlo + elsx - essl[nmov]
    eson1 = esol1 + es1 - ess1[nmov]
    eson2 = esol2 + es2 - ess2[nmov]
    enew = eold + e2 - e1 + esol - ess[nmov]

    if npols != 0:
        raise NotImplementedError('`polslv`')

    if icussl == 1:
        if nstyp[nmov] == modcus:
            raise NotImplementedError('`esvint`')

    if iewald == 1:
        raise NotImplementedError('`iewald`')

    if nopref == 0:
        wnsum = 0.0
        wnmax = 0.0
        for i in range(nmol):
            wnew[i] = wi[i] * wisum

        wnew[nmov] = 1.0 / (wik + wkc)

        for i in range(nmol):
            wnsum += wnew[i]

        delta = 1.0 / wnsum

        for i in range(nmol):
            wnew[i] *= delta
            wnmax = max(wnmax, wnew[i])


def ranint(nmr, nran, nvar):
    k = 0
    xvar = float(nvar)
    while True:
        nmr[k] = int(xvar * ranu()) + 1
        bo = True
        for i in range(k):
            if nmr[k] == nmr[i]:
                bo = False
                break
        if bo:
            k += 1
            if k >= nran: break
    return sorted(nmr)


def movmol(
    movvol, nrdfs, nrdl, nsolut, nfatm, nlatm, icapat, anew, asol, asol1, asol2, anew1, anew2,
    ncent1, indsol, rdels, adels, ncent2, islab, edge, edg2, nvdih, nvbnd, nvang, maxovl, nsatm,
    cori, corf, phinew, phi, ivdsum, maxvar, nflip, idih4, flip, dihdl, dihdl2, twopi, nvddep,
    idd, phine1, phine2, bndnew, bnd, ivbsum, ibnd1, nvbdep, bndne1, bnr2, bnr1, bnrd0, bndne2,
    angnew, ang, ivasum, angdl, angdl2, bnddl, bnddl2, ibd, iang3, nvadep, iad, angne1, angne2,
    angt1, angt0, isqm, edihne, edihol, ebndne, ebndol, eangne, eangol, eold, tk25, enbne, enbol,
    igbsa, exxnew, exxold, eponw, epoold, iewald, natmx, nmol, modsv1, nstyp, modsv2, emov, ecsx,
    elsx, essc, essl,  ess, ess1, ess2, es1, es2, wik,  rnew, rinc, rdlemin, wnew, wkc, nopref,
    esinol, angt2, emovc, emovl, esmov, esmov2, esonol, beta, lhtsol, betlht, vdl2, vdel, ivxyz,
    vold, ac, nsvat, noss, eij, **kws
):
    if movvol == 1:
        # goto 300
        return

    idists = [ [0.0 for i in range(nrdl)] for j in range(nrdfs)]

    isol = 0
    if nsolut > 1:
        while True:
            isol = int(ranu() * nsolut) + 1
            if nfatm[isol] == nlatm[isol]:
                if icapat == nfatm[isol]:
                    continue
            break

    for j in range(3):
        for i in range(len(asol)):
            anew[i][j] = asol[i][j]
            anew1[i][j] = asol1[i][j]
            anew2[i][j] = asol2[i][j]

    # Determine displacement parameters
    icent = ncent1
    if isol <= 1 and indsol == 0:
        n = 1
        m = nlatm[1]
        if nsolut == 1:
            m = nlatm[0]
        rdl = rdels[0]
        adl = adels[0]
    else:
        rdl = rdels[isol]
        adl = adels[isol]
        if isol <= 1:
            icent = ncent2
        elif isol > 1:
            icent = nfatm[isol]
        n = nfatm[isol]
        m = nlatm[isol]

    for i in range(3):
        rr = -rdl + ranu() * (2 * rdl)
        for j in range(n, m):
            anew[j][i] = asol[j][i] + rr
            anew1[j][i] = asol1[j][i] + rr
            anew2[j][i] = asol2[j][i] + rr

        if islab == 1 and i == 2:
            continue

        if anew[icent][i] <= edg2[i]:
            if anew[icent][i] < -edg2[i]:
                delta = edge[i]
                for j in range(n, m):
                    anew[j][i] += delta
                    anew1[j][i] += delta
                    anew2[j][i] += delta
        else:
            delta = -edge[i]
            for j in range(n, m):
                anew[j][i] += delta
                anew1[j][i] += delta
                anew2[j][i] += delta

    achg = -adl + ranu() * (adl + adl)
    if (m-n) != 0:
        rotate(achg)
    
    if (nvdih + nvbnd + nvang) == 0:
        if maxovl == 1:
            qtrf(anew, anew1, cori, nsatm, 0)
            qtrf(anew, anew2, corf, nsatm, 0)
            for i in range(nsatm):
                anew1[i][0] = cori[i][0]
                anew1[i][1] = cori[i][1]
                anew1[i][2] = cori[i][2]
                anew2[i][0] = corf[i][0]
                anew2[i][1] = corf[i][1]
                anew2[i][2] = corf[i][2]
        # goto 260
        return
    
    iflp = 0
    if nvdih != 0:
        for i in range(nvdih):
            phinew[i] = phi[i]
        
        m = ivdsum[isol]
        if m == 0:
            # goto 140
            return
        
        nran = m
        if m > 3:
            nran = 2 + int(ranu() * (m-2))
            if nran > maxvar:
                nran = maxvar
            ranint(nmr, nran, m)
        else:
            nmr = [1, 2, 3]
        
        knt = 0
        j = 0
        nflip2 = 2 * nflip
        if nflip2 < 6:
            nflip2 = 6
        
        for i in range(nvdih):
            if isolut(idih4[i]) != isol:
                continue
            knt += 1
            if knt != nmr[j]:
                continue
            j += 1
            if flip[i] != 0:
                raise NotImplementedError()
            
            phinew[i] = phi[i]
            phinew[i] += -dihdl[i] + ranu() * dihdl2[i]
            
            if phinew[i] > twopi:
                phinew[i] -= twopi
            elif phinew[i] < 0:
                phinew[i] += twopi
        
        if iflp == 1:
            nflptr = -nflptr
        
        if nvddep != 0:
            for i in range(nvddep):
                ii = nvdih - nvddep + i
                iii = idd[i]
                phinew[ii] = phinew[iii]
        
        for i in range(nvdih):
            phine1[i] = phinew[i]
            phine2[i] = phinew[i]
    
    if nvbnd != 0:
        for i in range(nvbnd):
            bndnew[i] = bnd[i]
        
        m = ivbsum[isol]
        if m == 0:
            # goto 190
            return
        
        nran = m
        if m > 3:
            nran = 2 + int(ranu() * (m-2))
            if nran > maxvar:
                nran = maxvar
            ranint(nmr, nran, m)
        else:
            nmr = [1, 2, 3]
        
        scalinv = (0.80 + 0.2 * (1.0 / nran**3)) / math.sqrt(nran)
        iscale = 1
        if nvang > 0:
            iscale += 1
        if nvdih > 0:
            iscale += 1
        scalinv = (scalinv * (0.8 + 0.2 * (1.0 / iscale**3))) / math.sqrt(iscale)
        
        knt = 0
        j = 0
        for i in range(nvbnd):
            if isolut(ibnd1[i]) != isol:
                continue
            knt += 1
            if knt != nmr[j]:
                continue
            j += 1
            bndnew[i] = bnd[i] + scalinv * (ranu() * bnddl2[i] - bnddl[i])
        
        if nvbdep != 0:
            for i in range(nvbdep):
                ii = nvbnd - nvbdep + i
                iii = ibd[i]
                bndnew[ii] = bndnew[iii]
        
        for i in range(nvbnd):
            bndne1[i] = bndnew[i] + (bnr1[i] - bnrd0[i])
            bndne2[i] = bndnew[i] + (bnr2[i] - bnrd0[i])
    
    if nvang != 0:
        for i in range(nvang):
            angnew[i] = ang[i]
        
        m = ivasum[isol]
        if m == 0:
            # goto 240
            return
        
        nran = m
        if m > 3:
            nran = 2 + int(ranu() * (m-2))
            if nran > maxvar:
                nran = maxvar
            ranint(nmr, nran, m)
        else:
            nmr = [1, 2, 3]
        
        scalinv = (0.80 + 0.2 * (1.0 / nran**3)) / math.sqrt(nran)
        iscale = 1
        if nvbnd > 0:
            iscale += 1
        if nvdih > 0:
            iscale += 1
        scalinv = (scalinv * (0.8 + 0.2 * (1.0 / iscale**3))) / math.sqrt(iscale)
        
        knt = 0
        j = 0
        for i in range(nvang):
            if isolut(iang3[i]) != isol:
                continue
            knt += 1
            if knt != nmr[j]:
                continue
            j += 1
            angnew[i] = ang[i] + scalinv * (ranu() * angdl2[i] - angdl[i])
        
        if nvadep != 0:
            for i in range(nvadep):
                ii = nvang - nvadep + i
                iii = iad[i]
                angnew[ii] = angnew[iii]
        
        for i in range(nvang):
            angne1[i] = angnew[i] + (angt1[i] - angt0[i])
            angne2[i] = angnew[i] + (angt2[i] - angt0[i])
    
    makml2()
    if isqm != 1:
        eintra(x, 1)
        del_value = edihne - edihol + ebndne - ebndol + eangne - eangol
        if del_value > tk25:
            enew = eold + tk25
            return
    
    eintra(x, 2)
    del_value += enbne - enbol
    headfile = '~+`#@!'
    xztype = '~+%#@'
    headfile = 'head'    
    if headfile == 'head  ':
        with open(headfile, 'r') as f:
            xztemp, xztype = f.read().split()
    
    if igbsa == 1:
        raise NotImplementedError()
    
    if del_value > tk25:
        enew = eold + tk25
        return

    xxpot()
    del_value += exxnew - exxold

    polpot()
    del_value += eponw - epoold
    
    if del_value > tk25:
        enew = eold + tk25
        return
    
    if iewald == 1:
        raise NotImplementedError()
    
    esone = 0.0
    esonc = 0.0
    esonl = 0.0
    eson1 = 0.0
    eson2 = 0.0
    
    if natmx == 0:
        # goto 290
        return
    
    for i in range(nmol):
        modsv = modsv1
        ntyp = nstyp[i]
        
        if ntyp == 2:
            modsv = modsv2
        
        if modsv > 2:
            e = sxpot[i]
        else:
            e = wxpot[i]
        
        emov[i] = e
        emovc[i] = ecsx
        emovl[i] = elsx
        esone += e
        esonc += ecsx
        esonl += elsx
        esmov[i] = es1
        eson1 += es1
        esmov2[i] = es2
        eson2 += es2
        wnew[i] = wik
        
        if ntyp == 2:
            continue
        
        for j in range(nrdfs):
            k = int(rinc * (rnew[j] - rdlemin)) + 1
            if k <= nrdl:
                if k >= 1:
                    idists[k, j] += 1
    
    enew = eold + esone - esonol + del_value
    
    if qmname == 'g09u' and xztype == 'false':
        tk20 = 20.0 / beta
        del_value += esone - esonol
        if del_value < -tk20:
            qmname = 'g091'
        elif del_value <= tk20:
            x = ranu()
            xb = beta
            if isol == lhtsol:
                xb = betlht
            if lhtsol == 9999:
                xb = betlht
            emet = math.exp(-xb * del_value)
            if emet >= x:
                qmname = 'g092'
            else:
                qmname = 'g09d'
        else:
            qmname = 'g09l'
        
        if qmname == 'g091' or qmname == 'g092':
            # goto 255
            return
    
    if natmx == 0:
        # goto 470
        return
    
    delv = ranu() * vdl2 - vdel
    if ivxyz == 1:
        delv *= 0.75
    
    vnew = vold + delv
    nmovsv = nmov
    
    # Copy edge values
    oldedge = edge[:]
    
    if ivxyz != 1:
        fac = (vnew / vold) ** (1.0 / 3.0)
        slvfac = fac - 1.0
        
        for i in range(3):
            edge[i] = fac * edge[i]
            edg2[i] = 0.5 * edge[i]
        
        for j in range(3):
            for i in range(nmol):
                tmp = slvfac * ac[i][0][j]
                natoms = nsvat[nstyp[i]]
                
                for k in range(natoms):
                    ac[i][k][j] += tmp
                    
                soldel[j] = slvfac * (asol[ncent1][j] + asol[ncent2][j]) / 2.0
                tmp = soldel[j]
                
                for i in range(nsatm):
                    asol[i][j] += tmp
                    asol1[i][j] += tmp
                    asol2[i][j] += tmp
    else:
        if islab != 1:
            ivax = int(3.0 * ranu()) + 1
        else:
            ivax = int(2.0 * ranu()) + 1
        
        edge[ivax] = delv * edge[ivax] / vold
        edg2[ivax] = 0.5 * edge[ivax]
        
        slvfac = vnew / vold - 1.0
        
        for i in range(nmol):
            tmp = slvfac * ac[i][0][ivax]
            natoms = nsvat[nstyp[i]]
            for k in range(natoms):
                ac[i][k][ivax] += tmp
        
        soldel = [0.0, 0.0, 0.0]
        tmp = slvfac * (asol[ncent1][ivax] + asol[ncent2][ivax]) / 2.0
        soldel[ivax] += tmp
        for i in range(nsatm):
            asol[i][ivax] += tmp
            asol1[i][ivax] += tmp
            asol2[i][ivax] += tmp
    
    esone = esonc = esonl = eson1 = eson2 = enew = eone = 0.0
    for j in range(nrdfs):
        for i in range(nrdl):
            idists[i][j] = 0
    
    if noss == 1:
        # goto 425
        return
    
    knt = 1
    n = nmol - 1
    
    for i in range(n):
        k = i + 1
        nmov = i
        ntyp = modsv1
        
        if nstyp[i] == 1:
            ntyp = modsv2
        
        natoms = nsvat[nstyp[i]]
        for m in range(3):
            for kk in range(natoms):
                anew[kk][m] = ac[i][kk][m]
        
        for j in range(k, nmol):
            if nsvat[2] == 0:
                if modsv1 > 2:
                    e = sspot(j)
                else:
                    e = wwpot(j)
            else:
                jtyp = modsv1
                if nstyp[j] == 2:
                    jtyp = modsv2
                if ntyp <= 2 and jtyp <= 2:
                    e = wwpot(j)
                else:
                    e = sspot(j)
            
            eij[knt] = e
            enew += e
            knt += 1
            
            if i == nmovsv:
                emov[j] = e
            else:
                if j != nmovsv:
                    continue
                emov[i] = e
            eone += e
    
    movtyp = 1
    for j in range(3):
        for i in range(nsatm):
            anew[i][j] = asol[i][j]
            anew1[i][j] = asol1[i][j]
            anew2[i][j] = asol2[i][j]
    
    for i in range(nmol):
        modsv = modsv1
        ntyp = nstyp[i]
        
        if ntyp == 2:
            modsv = modsv2
        
        if modsv > 2:
            e = sxpot(i)
        else:
            e = wxpot(i)
        
        enew += e
        esone += e
        esonc += ecsx
        esonl += elsx
        eson1 += es1
        ess1[i] = es1
        eson2 += es2
        ess2[i] = es2
        wnew[i] = wik
        ess[i] = e
        essc[i] = ecsx
        essl[i] = elsx
        
        if ntyp == 2:
            continue
        
        for m in range(nrdfs):
            l = int(rinc * (rnew[m] - rdlemin)) + 1
            if l <= nrdl:
                if l >= 1:
                    idists[l][m] += 1
    
    enew += exxold + edihol + enbol + ebndol + eangol + esinol
    polpot()
    enew += eponw
    enew += ecut()
    
    if iewald == 1:
        raise NotImplementedError()
    
    if nopref == 0:
        wnsum = 0.0
        wnmax = 0.0
        
        for i in range(nmol):
            wnew[i] = 1.0 / (wnew[i] + wkc)
            wnsum += wnew[i]
        
        tmp = 1.0 / wnsum
        
        for i in range(nmol):
            wnew[i] *= tmp
            wnmax = max(wnmax, wnew[i])
    
    return


def enmtx(
    nmol, noss, icussl, igbsa,
    #slvbnd, nvbnd, nvang, nvdih, nabnd, naang, nadih,       # for `icussl==1`
    #modcus, nvsbnd, nvsang, nvsdih,
    **kws
):
    iewald = kws['iewald']
    ac = kws['ac']
    nstyp = kws['nstyp']
    nsvat = kws['nsvat']
    natmx = kws['natmx']
    asol = kws['asol']
    asol1 = kws['asol1']
    asol2 = kws['asol2']
    modsv1 = kws['modsv1']
    modsv2 = kws['modsv2']

    kws['eij'] = eij = [0.0 for i in range(nmol*(nmol-1))]

    enew = esbnne = esanhe = esdine = esnbne = esinne = 0.0

    if natmx != 0:  # done
        kws['movtyp'] = 0
        knt = 0
        for i in range(nmol):
            kws['nmov'] = i
            jtyp = nstyp[i]
            kws['natoms'] = natoms = nsvat[jtyp]        # special
            ntyp = modsv1 if jtyp != 2 else modsv2

            kws['anew'] = [[ac[i][l][m] for m in range(3)] for l in range(natoms)]

            for j in range(i+1, nmol):
                kws['nm'] = j
                if noss == 1:
                    pass
                elif nsvat[1] == 0: # not
                    if modsv1 > 2:
                        ckws = sspot(**kws)
                    else:
                        ckws = wwpot(**kws)
                else:   # done
                    jtyp = modsv1 if nstyp[j] != 1 else modsv2
                    if ntyp <= 2 and jtyp <= 2:
                        ckws = wwpot(**kws)
                    else:
                        ckws = sspot(**kws)

                e = 0.0
                eij[knt] = e
                enew += e
                knt += 1

        if icussl == 1:
            raise NotImplementedError('icussl==1')

    kws['anew']  = kws['solutesdata']['reference']['xyz']
    kws['anew1'] = kws['solutesdata']['first']['xyz']
    kws['anew2'] = kws['solutesdata']['second']['xyz']

    if igbsa == 1:
        raise NotImplementedError('igbsa==1')

    if natmx != 0:
        kws['movtyp'] = 1
        wkc = kws['wkc']
        essmin = kws['essmin']
        essinc = kws['essinc']
        wi = [0.0 for i in range(nmol)]
        esone = esonc = esonl = eson1 = eson2 = wisum = 0.0
        kws['ess'] = ess = [0.0 for i in range(nmol)]
        kws['ess1'] = ess1 = [0.0 for i in range(nmol)]
        kws['ess2'] = ess2 = [0.0 for i in range(nmol)]
        kws['essc'] = essc = [0.0 for i in range(nmol)]
        kws['essl'] = essl = [0.0 for i in range(nmol)]
        kws['ioess'] = ioess = [[0 for i in range(50)] for i in range(2)]
        for i in range(nmol):
            kws['nm'] = i
            modsv = modsv1 if nstyp[i] != 1 else modsv2
            if modsv > 2:
                ckws = sxpot(**kws)
            else:
                ckws = wxpot(**kws)
            wik = 1.0   # ckws['wik']
            wik = 1.0 / (wik + wkc)
            wi[i] = wik
            wisum += wik
            ess[i] = e
            essc[i] = 0.0 # ecsx
            essl[i] = 0.0 # elsx
            j = int((e - essmin) / essinc) + 1
            if 1 <= j <= 50:
                ioess[j][0] += 1
            eson1 += 0.0 # es1
            ess1[i] = 0.0 #es1
            eson2 += 0.0 # es2
            ess2[i] = 0.0 # es2
            esone += e
            esonc += 0.0 # ecsx
            esonl += 0.0 # elsx
        enew += esone
        wimax = 0.0
        for i in range(nmol):
            wi[i] /= wisum
            wimax = max(wimax, wi[i])

    movtyp = 0

    if nvbnd != 0:
        for i in range(nvbnd):
            x = bnd[i]
            bndnew[i] = x
            bndne1[i] = x + (bndr1[i] - bndr0[i])
            bndne2[i] = x + (bndr2[i] - bndr0[i])

    if nvang != 0:
        for i in range(nvang):
            x = ang[i]
            angnew[i] = x
            angne1[i] = x + (angt1[i] - angt0[i])
            angne2[i] = x + (angt2[i] - angt0[i])

    if nvdih != 0:
        for i in range(nvdih):
            x = phi[i]
            phinew[i] = phine1[i] = phine2[i] = x

    x = 0.0
    eintra(x, 0)
    enew += x

    if nabnd != 0:
        for i in range(nvbnd, ntbnd):
            bnd[i] = bndnew[i]

    if naang != 0:
        for i in range(nvang, ntang):
            ang[i] = angnew[i]

    if nadih != 0:
        for i in range(nvdih, ntdih):
            phi[i] = phinew[i]

    xxpot()
    enew += exxnew

    #polpot()
    #enew += eponew

    enew += ecut()

    if iewald == 1:
        raise NotImplementedError('iewald==1')



