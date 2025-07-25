
ESQ = 332.06
QT3 = [-0.834,  0.417,  0.417,   0.0]
QT4 = [   0.0,  0.52,   0.52,  -1.04]
SOLVENT_MODE = {
    1: {
        'solvent_name': 'TIP3P',
        'aoo': 582000.0,
        'coo': 595.0,
        'kq': 0,        # coulombic interaction index starts at 0
        'qq': [
            ESQ*QT3[0]*QT3[0], ESQ*QT3[0]*QT3[1], ESQ*QT3[0]*QT3[1],
            ESQ*QT3[0]*QT3[1], ESQ*QT3[1]*QT3[1], ESQ*QT3[1]*QT3[1],
            ESQ*QT3[0]*QT3[1], ESQ*QT3[1]*QT3[1], ESQ*QT3[1]*QT3[1],
        ]
    },
    2: {
        'solvent_name': 'TIP4P',
        'aoo': 600000.0,
        'coo': 610.0,
        'kq': 1,        # coulombic interaction index starts at 1, notice charge on OW is 0.0
        'qq': [
            ESQ*QT4[1]*QT4[1], ESQ*QT4[1]*QT4[1], ESQ*QT4[1]*QT4[3],    # 1 -> 1,2,3    index
            ESQ*QT4[1]*QT4[1], ESQ*QT4[1]*QT4[1], ESQ*QT4[1]*QT4[3],    # 2 -> 1,2,3
            ESQ*QT4[1]*QT4[3], ESQ*QT4[1]*QT4[3], ESQ*QT4[3]*QT4[3],    # 3 -> 1,2,3
        ]
    },
}


#    key:     [Name,    modenum,    atom_names,                  atom_indexes_in_parfile]
SOLVENT_PARS = {
    'tip3p' : ['TIP3P',    1,    ('OW','HW','HW'),               (111,112,112,100)      ],
    'tip4p' : ['TIP4P',    2,    ('OW','HW','HW','M'),           (113,114,114,115)      ],
    'ch3oh' : ['CH3OH',    3,    ('O','H','CH3'),                (78,79,80,100)         ],
    'ch3cn' : ['CH3CN',    4,    ('C','N','CH3'),                (95,94,96,100)         ],
    'meome' : ['MeOMe',    5,    ('O','CH3','CH3'),              (108,109,109,100)      ],
    'c3h8'  : ['C3H8',     6,    ('CH2','CH3','CH3'),            (71,68,68,100)         ],
    'chcl3' : ['CHCl3',    7,    ('CH','CL','CL','CL'),          (120,121,121,121)      ],
    'mecl2' : ['MeCl2',    8,    ('CH2','CL','CL'),              (118,119,119,100)      ],
    'thf'   : ['THF',      9,    ('O','CH2','CH2','CH2','CH2'),  (108,110,71,110,71)    ],
    'argon' : ['Argon',   10,    ('AR',),                        (103,)                 ],
    'ccl4'  : ['CCl4',    11,    ('C','CL','CL','CL','CL'),      (122,123,123,123,123)  ],
    'dmso'  : ['DMSO',    12,    ('S','O','CH3','CH3'),          (124,125,126,126)      ],
    'tip5p' : ['TIP5P',   13,    ('O','HW','HW','LP','LP'),      (97,98,98,99,99)       ],
    # special cases
    ''      : ['None',     0,    tuple(),                        tuple()                ],
    'none'  : ['None',     0,    tuple(),                        tuple()                ],
    'other' : ['OTHER',   20,    tuple(),                        tuple()                ],
    'gbsa'  : ['GB/SA',   99,    tuple(),                        tuple()                ],
    'gb/sa' : ['GB/SA',   99,    tuple(),                        tuple()                ],
}


PAR_SETTINGS = {
    # key       value   explanation
    'svmod1' :  ('',    'primary solvent name'),
    'xstren' :  (0.0,   'XSTREN'),
    'qmname' :  ('',    'quantum name'),
    'cmname' :  ('',    'calculation name'),
    'qmscale':  (0.0,   'scale factor'),
    'ich0'   :  (0,     'net charge for Reference solute'),
    'ich1'   :  (0,     'net charge for First solute'),
    'ich2'   :  (0,     'net charge for Second solute'),
    'slfmt'  :  ('',    'solute format'),
    'solor'  :  ('',    'source to get solvent coordinate'),
    'slab'   :  ('',    'solvent origin line to remove Z-periodicity'),
    'vxyz'   :  ('',    'whether move volume independently in X, Y and Z directions'),
    'icapat' :  (0,     'solute atom index (from 1) to define center of cap'),
    'caprad' :  (0.0,   "cap's radius"),
    'capfk'  :  (0.0,   "cap's force constant, unit in kcal/mol/A^2"),
    'optim'  :  ('',    'solute optimization'),
    'ftol'   :  (0.0,   'optimization convergence criterion, unit in kcal/mol'),
    'newzm'  :  ('',    'where to save optimization results'),
    'nsaevl' :  (0,     'maximum number of potential energy function evaluations in Simulated Annealing'),
    'nsatr'  :  (0,     'number of iterations before temperature reduction in Simulated Annealing'),
    'satemp' :  (0.0,   'temperature for Simulated Annealing'),
    'sart'   :  (0.0,   'temperature reduction factor in simulated annealing'),
    'icsopt' :  ((0,0,0,0,0,0,0,0), 'option settings for conformational searching (8I)'),
    'csetol' :  (0.0,   'energy criterion for a unique conformer in conformational searching'),
    'csrms'  :  (0.0,   'RMS tolerance in superposition of conformers in conformational searching'),
    'csrtol' :  (0.0,   'distance criterion for a unique conformer in conformational searching'),
    'csreup' :  (0.0,   'upper bound energy relative to the current lowest minimum in conformational searching'),
    'csaelo' :  (0.0,   'absolute lower bound energy in conformational searching'),
    'csaeup' :  (0.0,   'absolute upper bound energy in conformational searching'),
    'nmode'  :  (0,     'if it is set to 1 and ICALC is 2, 3, or 4, a normal coordinate analysis is performed'),
    'nsym'   :  (0,     'symmetry number in normal mode analysis'),
    'nconf'  :  (0,     'number of equivalent conformers in normal mode analysis'),
    'fscale' :  (0.0,   'frequency scale factor in normal mode analysis'),
    'frqcut' :  (0.0,   'cutoff frequency in normal mode analysis'),
    'nmol'   :  (0,     'number of primary solvent molecules'),
    'ibox'   :  (0,     'number of solvent molecules got from pre-stored box'),
    'boxcut' :  (0.0,   'extension number of the farthest solute atom in the x, y, z directions'),
    'svmod2' :  ('',    'second solvent name'),
    'nmol2'  :  (0,     'number of second solvent molecules'),
    'ncent1' :  (0,     'atom index used to define the centers of solutes 1'),
    'ncent2' :  (0,     'cutoff or atom index based on `icut`'),
    'ncents' :  (0,     'atom index used in a custom solvent molecule for cutoffs'),
    'icutas' :  ((0,0,0,0,0),   'solvent atom index used for solvent-solvent, solute-solvent cutoff (5I)'),
    'nrota1' :  (0,     'atom index for solutes 1 rotated about (rigid-body)'),
    'nrota2' :  (0,     'atom index for solutes 2 rotated about (rigid-body)'),
    'icut'   :  (0,     'cutoff mode'),
    'icutat' :  ((0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),   'solute atom indexes used for `icut=4` (15I)'),
    'irecnt' :  (0,     'whether recenter solute in middle of box'),
    'indsol' :  (0,     'whether let move solute move independently'),
    'izlong' :  (0,     'whether oriented solute in the longest axis along the Z-direction'),
    'maxovl' :  (0,     'for FEP, whether let perturbed solutes maximally overlaid on the reference solute'),
    'noxx'   :  (0,     'whether not evaluate solute-solute interactions'),
    'noss'   :  (0,     'whether not evaluate solvent-solvent interactions'),
    'nobndv' :  (0,     'whether not have bond variations'),
    'noangv' :  (0,     'whether not have angle variations'),
    'nonebn' :  (0,     'whether not determine covalent neighbors by using interatomic distance, use inputs instead'),
    'nrdf'   :  (0,     'number for solvent-solvent RDFs, unique atom RDFs is only allowed upto this number'),
    'nrdfs'  :  (0,     'number for solute-solvent RDFs, unique atom RDFs is only allowed upto this number'),
    'nrdfs1' :  ((0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),   'for solvent-solvent RDFs (15I)'),
    'nrdfs2' :  ((0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),   'for solvent-solvent RDFs (15I)'),
    'nrdfa1' :  ((0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),   'for solute-solvent RDFs (15I)'),
    'nrdfa2' :  ((0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),   'for solute-solvent RDFs (15I)'),
    'rdlmin' :  (0.0,   'minimum for RDF'),
    'rdlinc' :  (0.0,   'increment for RDF'),
    'eprmin' :  (0.0,   'minimum for solvent-solvent pair distribution'),
    'eprinc' :  (0.0,   'increment for solvent-solvent pair distribution'),
    'edmin'  :  (0.0,   'minimum for solvent-solvent energy distribution'),
    'edinc'  :  (0.0,   'increment for solvent-solvent energy distribution'),
    'essmin' :  (0.0,   'minimum for solute-solvent energy distribution'),
    'essinc' :  (0.0,   'increment for solute-solvent energy distribution'),
    'ebsmin' :  (0.0,   'minimum for solute-solvent pair distribution'),
    'ebsinc' :  (0.0,   'increment for solute-solvent pair distribution'),
    'nvchg'  :  (0,     'frequency for volume move'),
    'nschg'  :  (0,     'frequency for solute move'),
    'maxvar' :  (0,     'maximum number for bond/angle/dihedral variations'),
    'nconsv' :  (0,     'frequency for saving full system coordinates'),
    'nbuse'  :  (0,     'whether to use neighbor lists for evaluating solvent-solvent interactions'),
    'nocut'  :  (0,     'whether not to add cutoff correction for solvent-solvent LJ interactions'),
    'nosmth' :  (0,     'whether not to smooth solute-solvent and solvent-solvent interactions near cutoff'),
    'vdel'   :  (0.0,   'maximum volume change in angstrom^3'),
    'wkc'    :  (0.0,   'constant for preferential sampling using `1.0/(r^2+wkc)`'),
    'rdel'   :  (0.0,   'solvent translation in angstrom'),
    'adel'   :  (0.0,   'solvent rotation in degree'),
    'rsolv'  :  (0.0,   'solvent radius for SASA calculation'),
    'rdels1' :  (0.0,   'translation range in angstrom for solute 1'),
    'adels1' :  (0.0,   'rotation range in degree for solute 1'),
    'rdels2' :  (0.0,   'translation range in angstrom for solute 2'),
    'adels2' :  (0.0,   'rotation range in degree for solute 1'),
    'lhtsol' :  (0,     'solute index to be heated at temperature `tlht`'),
    'tlht'   :  (0.0,   'temperature'),
    'rcut'   :  (0.0,   'cutoff for solvent-solvent interactions'),
    'scut'   :  (0.0,   'cutoff for solute-solvent interactions'),
    'cutnb'  :  (0.0,   'cutoff for intrasolute nonbonded pairs'),
    'iewald' :  (0,     'whether use Ewald summation'),
    't'      :  (0.0,   'external temperature in Celcius'),
    'p'      :  (0.0,   'external pressure in atm'),
    'torcut' :  (0.0,   'cutoff for individual torsional energy term'),
    'scllj'  :  (0.0,   'factor to scale solute LJ sigma and epsilon values'),
    'dielec' :  (0.0,   'dielectric constant'),
    'scl14c' :  (0.0,   'factor to scale Coulombic 1-4 intramolecular nonbonded interactions'),
    'scl14l' :  (0.0,   'factor to scale LJ 1-4 intramolecular nonbonded interactions'),
    'dielrf' :  (0.0,   'cutoff for reaction field'),
    'pltfmt' :  ('',    'format to output for solute and solvent coordinates'),
    'isolec' :  (0.0,   'solute index where to decompose solute-solvent and solute-solute Coulomb and LJ energies'),
    'polscx' :  (0.0,   'include solute-solute polarization'),
    'polscs' :  (0.0,   'include custom solvent polarization'),
    # missing pars
    'nrdl'   :  (100,   'number of bins for each RDF histogram'),
}


AMBER_SURFACE_AREA = {
    'CT':1, 'C':2,  'CA':3, 'CM':3, 'CC':3, 'CV':3, 'CW':3, 'CR':3, 'CB':3, 'C*':3,
    'CN':3, 'CK':3, 'CQ':3, 'N':2,  'NA':2, 'NB':2, 'NC':2, 'N*':2, 'N2':2, 'N3':2,
    'OW':2, 'OH':2, 'OS':2, 'O':2,  'O2':2, 'S':4,  'SH':4, 'P':4,  'H':2,  'HW':2,
    'HO':2, 'HS':2, 'HA':3, 'HC':1, 'H1':1, 'H2':1, 'H3':1, 'H4':3, 'H5':3, 'HP':1,
    'OA':4, 'C=':3, 'CO':1, 'NT':2, 'CY':1, 'SY':4, 'OY':2, 'SA':4, 'CZ':1, 'NZ':2,
    'C2':1, 'C3':1, 'F':4,  'Cl':4, 'Br':4, 'I':4,  'LP':2, 'CS':3, 'CU':3, 'CH':1,
    'SZ':4, 'CX':3, 'C4':1, 'C7':1, 'C8':1, 'C9':1, 'Be':2, 'B':2,  'Al':2, 'Si':1,
    'Li':2, 'Na':2, 'K':2,  'Rb':2, 'Cs':2, 'Mg':2, 'Ca':2, 'Sr':2, 'Ba':2, 'DM':1,
    'Cu':2, 'Zn':2, 'He':1, 'Ne':1, 'Ar':1, 'Kr':1, 'Xe':1, 'NY':2, 'NO':2, 'ON':2,
    'C$':2, 'N$':2, 'C!':3, 'O$':2, 'S=':4, 'C+':2, 'P+':2, 'N=':2, 'NM':2, 'C:':3,
    'C#':3, 'XC':4, 'XB':4, 'XI':4, 'NS':2, 'NX':2, 'CP':3, 'CF':1,
    # special cases
    '':1,   ' ':1,  '  ':1,
}


# atomic number 0 is the special case, halogen XP site
ATOMIC_RADIUS = {
    0:0.44 ,  1:0.32,  2:0.93,  3:1.23,  4:0.9 ,  5:0.82,  6:0.77,  7:0.73,  8:0.7 ,   9:0.6,
    10:0.71, 11:1.54, 12:1.36, 13:1.18, 14:1.11, 15:1.06, 16:1.02, 17:1.0 , 18:0.98, 19:2.03,
    20:1.74, 21:1.44, 22:1.32, 23:1.22, 24:1.18, 25:1.17, 26:1.17, 27:1.16, 28:1.15, 29:1.17,
    30:1.25, 31:1.26, 32:1.22, 33:1.2 , 34:1.16, 35:1.16, 36:1.12, 37:2.16, 38:1.91, 39:1.62,
    40:1.45, 41:1.34, 42:1.3 , 43:1.27, 44:1.25, 45:1.25, 46:1.28, 47:1.34, 48:1.48, 49:1.44,
    50:1.41, 51:1.4 , 52:1.36, 53:1.33, 54:1.31, 55:2.35, 56:1.98, 57:1.69, 58:1.65, 59:1.65,
    60:1.64, 61:1.63, 62:1.62, 63:1.85, 64:1.61, 65:1.59, 66:1.59, 67:1.58, 68:1.57, 69:1.56,
    70:1.56, 71:1.56, 72:1.44, 73:1.34, 74:1.3 , 75:1.28, 76:1.26, 77:1.27, 78:1.3 , 79:1.34,
    80:1.49, 81:1.48, 82:1.47, 83:1.46, 84:1.46, 85:1.45, 86:1.45, 87:1.5 , 88:1.5 ,  89:1.6,
    90:1.65, 91:1.7 ,  92:1.7,  93:1.7,  94:1.7,  95:1.7,  96:1.7, 97:1.7 , 98:1.7 ,  99:0.0,
}


# sequence on periodic table, starting from Hydrogen, index 0
ATOMIC_WEIGHTS = [
    1.00790,   4.00260,   6.94100,   9.01218,   10.81000,
    12.01100,  14.00670,  15.99940,  18.998403, 20.17900,  22.98977,  24.30500,
    26.98154,  28.08550,  30.97376,  32.06000,  35.45300,  39.94800,  39.09830,
    40.08000,  44.95590,  47.90000,  50.94150,  51.99600,  54.93800,  55.84700,
    58.93320,  58.70000,  63.54600,  65.38000,  69.73500,  72.59000,  74.92160,
    78.96000,  79.90400,  83.80000,  85.46780,  87.62000,  88.90590,  91.22000,
    92.90640,  95.94000,  98.90620,  101.1700,  102.9055,  106.4000,  107.8680,
    112.4100,  114.8200,  118.6900,  121.7500,  127.6000,  126.9045,  131.3000,
    132.9054,  137.3300,  138.9055,  140.1200,  140.9077,  144.2400,  145.0000,
    150.4000,  151.9600,  157.2500,  158.9254,  162.5000,  164.9034,  167.2600,
    168.9342,  173.0400,  174.9670,  178.4900,  180.9479,  183.8500,  186.2070,
    190.2000,  192.2200,  195.0900,  196.9665,  200.5900,  204.3700,  207.2000,
    208.9804,  209.0000,  210.0000,  222.0000,  223.0000,  226.0254,  227.0000,
    232.0381,  231.0359,  238.0290,  237.0482,  244.0000,  243.0000,  247.0000,
    247.0000,  251.0000,  254.0000,  257.0000,  258.0000,  259.0000,  260.0000,
]


AMBER_TO_OPLSAA_SYNONYM_ATOM_TYPE = {
    'CT':'CT',  'C ':'C ',  'CA':'CA',  'CM':'CA',  'CC':'CA',  'CV':'CA',  'CW':'CA',
    'CR':'CA',  'CB':'CA',  'C*':'CA',  'CN':'CA',  'CK':'CA',  'CQ':'CA',  'N ':'N ',
    'NA':'NA',  'NB':'NC',  'NC':'NB',  'N*':'NA',  'N2':'NA',  'N3':'CT',  'OW':'OS',
    'OH':'OH',  'OS':'OS',  'O ':'O ',  'O2':'O ',  'S ':'S ',  'SH':'SH',  'P ':'S ',
    'H ':'H ',  'HW':'H' ,  'HO':'HO',  'HS':'HS',  'HA':'HA',  'HC':'HC',  'H1':'HC',
    'H2':'HC',  'H3':'HC',  'H4':'HA',  'H5':'HA',  'HP':'HC',  'C+':'CA',  'C=':'CM',
    'CO':'CT',  'NT':'CT',  'CY':'CT',  'SY':'S ',  'OY':'O ',  'P+':'P ',  'CZ':'C:',
    'NZ':'NC',  'C2':'CT',  'C3':'CT',  'F ':'F ',  'Cl':'Cl',  'Br':'Br',  'I ':'I ',
    'LP':'LP',  'CS':'CA',  'CU':'CV',  'CH':'CT',  'SZ':'SY',  'CX':'CA',  'C4':'CT',
    'C7':'CM',  'C8':'CM',  'C9':'CM',  'Be':'CZ',  'B ':'CA',  'Al':'S ',  'Si':'S ',
    'NY':'N2',  'NO':'N ',  'ON':'O ',  'C$':'C' ,  'N$':'N ',  'C!':'CA',  'O$':'OS',
    'S=':'SY',  'N=':'NC',  'NM':'N ',  'C:':'CZ',  'C#':'CM',  'XC':'XC',  'XB':'XB',
    'XI':'XI',  'NS':'NA',  'NX':'NA',  'CP':'CW',  'CF':'CT',  'OA':'OS',  'SA':'S ',
}


BONDED_ATOMIC_RADIUS = {
    'CT':0.765,  'C ':0.67 ,  'CA':0.7  ,  'CM':0.67,  'CC':0.7  ,  'CV':0.7  ,  'CW':0.7  ,
    'CR':0.7  ,  'CB':0.7  ,  'C*':0.7  ,  'CN':0.7 ,  'CK':0.7  ,  'CQ':0.7  ,  'N ':0.665,
    'NA':0.68 ,  'NB':0.69 ,  'NC':0.64 ,  'N*':0.71,  'N2':0.64 ,  'N3':0.71 ,  'OW':0.61 ,
    'OH':0.64 ,  'OS':0.64 ,  'O ':0.56 ,  'O2':0.58,  'S ':1.045,  'SH':1.045,  'P ':0.97 ,
    'H ':0.34 ,  'HW':0.34 ,  'HO':0.31 ,  'HS':0.3 ,  'HA':0.38 ,  'HC':0.33 ,  'H1':0.33 ,
    'H2':0.33 ,  'H3':0.33 ,  'H4':0.38 ,  'H5':0.38,  'HP':0.33 ,  'C+':0.7  ,  'C=':0.74 ,
    'CO':0.76 ,  'NT':0.68 ,  'CY':0.75 ,  'F ':0.58,  'Cl':1.02 ,  'Br':1.18 ,  'I ':1.38 ,
    'LP':0.0  ,  'SY':1.02 ,  'OY':0.42 ,  'CZ':0.61,  'NZ':0.55 ,  'P+':1.055,  'CH':0.765,
    'C2':0.765,  'C3':0.765,  'CS':0.74 ,  'CU':0.74,  'SZ':1.02 ,  'CX':0.7  ,  'C4':0.765,
    'C7':0.67 ,  'C8':0.67 ,  'C9':0.67 ,  'Be':0.9 ,  'B ':0.83 ,  'Al':1.18 ,  'Si':1.11 ,
    'NY':0.64 ,  'NO':0.72 ,  'ON':0.5  ,  'C$':0.67,  'N$':0.665,  'C!':0.72 ,  'O$':0.7  ,
    'S=':0.96 ,  'N=':0.64 ,  'NM':0.665,  'C:':0.63,  'C#':0.67 ,  'XC':0.58 ,  'XB':0.42 ,
    'XI':0.42 ,  'NS':0.68 ,  'NX':0.68 ,  'CP':0.7 ,  'CF':0.765,  'OA':0.66 ,  'SA':1.02 ,
}


ATOMIC_ELECTRONEGATIVITY = [
    2.2,  0.0,  1.0,  1.5,  2.0,  2.6, 3.05, 3.5,  3.9,  0.0,
    0.9,  1.2,  1.5,  1.9,  2.15, 2.6, 3.15, 0.0,  0.8,  1.0,
    1.3,  1.5,  1.6,  1.6,  1.5,  1.8, 1.8,  1.8,  1.9,  1.6,
    1.6,  1.9,  2.0,  2.45, 2.85, 0.0, 0.8,  1.0,  1.3,  1.6,
    1.6,  1.8,  1.9,  2.2,  2.2,  2.2, 1.9,  1.7,  1.7,  1.8,
    2.05, 2.30, 2.65, 0.0
]


PARS_CMD = {
    # key         default           explanation
    'rc0'       : (0.0,             'reference lambda'),
    'rc1'       : (0.5,             'perturbation double-window first lambda'),
    'rc2'       : (1.0,             'perturbation double-window second lambda'),
    'ncons'     : (0,               'number of total configurations'),
    'icalc'     : (1,               'type of calculation'),
    'inew'      : (1,               'indicator, whether a new run'),
    'iprint'    : (1,               'level of print setting'),
    'infile'    : ('fepin',         'file to save intermediate results for continuous run'),
    'upfile'    : ('fepup',         'file to save intermediate results for continuous run'),
    'save'      : ('fepsv',         'file to save final results'),
    'average'   : ('fepav',         'file to save averaged results'),
    'slvzmat'   : ('slvzmat',       'file to custom solvent zmatrix'),
    'bangpar'   : ('oplsaa.sb',     'file for bond/angle parameters'),
    'waterbox'  : ('watbox',        'file for water boxes'),
    'org1box'   : ('org1box',       'file for non aqueous water box 1'),
    'org2box'   : ('org2box',       'file for non aqueous water box 2'),
    'summary'   : ('fepsum',        'file for summary reults'),
    'output'    : ('fepota',        'file for output file'),
    'pltfile'   : ('fepplta',       'file for final coordiante file'),
    'iflink'    : ('gausslink.E01', 'file for intermediate scripts'),
    'zmatrix'   : ('pmfzmat',       'file for solute zmatrix'),
    'pars'      : ('pmfpar',        'file for parameter settings'),
    'parameter' : ('tmppar',        'file for compatible option, concatenation of @pars & @bangpar'),
}


# collection of energies that will be finally calculated
PAR_ENERGIES = {
    # key         default   explanation
    'esone'     : (0.0,     'solute-solvent LJ and coulombic interactions'),
    'esonco'    : (0.0,     'solute-solvent coulombic interactions'),
    'esonlo'    : (0.0,     'solute-solvent LJ interactions'),
    'esonol'    : (0.0,     'solute-solvent LJ and coulombic interactions -- old (esone)'),
    'esx'       : (0.0,     'solute-solvent LJ and coulombic interactions -- update (esonol)'),
    'esx1'      : (0.0,     'solute-solvent LJ and coulombic interactions'),
    'esx2'      : (0.0,     'solute-solvent LJ and coulombic interactions'),
    'esxco'     : (0.0,     'solute-solvent coulombic interactions'),
    'esxlj'     : (0.0,     'solute-solvent LJ interactions'),

    'exx'       : (0.0,     'solute-solute interactions'),
    'exx1'      : (0.0,     'solute-solute interactions'),
    'exx2'      : (0.0,     'solute-solute interactions'),
    'exxco'     : (0.0,     'solute-solute coulombic interactions'),
    'exxlj'     : (0.0,     'solute-solute LJ interactions'),

    'exxnew'    : (0.0,     'solute LJ and coulombic interactions'),
    'exxne1'    : (0.0,     'solute LJ and coulombic interactions'),
    'exxne2'    : (0.0,     'solute LJ and coulombic interactions'),
    'exxnel'    : (0.0,     'solute LJ interactions'),
    'exxnec'    : (0.0,     'solute coulombic interactions'),
    'exxold'    : (0.0,     'solute LJ and coulombic interactions -- old (exxnew)'),
    'exxol1'    : (0.0,     'solute LJ and coulombic interactions -- old (exxne1)'),
    'exxol2'    : (0.0,     'solute LJ and coulombic interactions -- old (exxne2)'),
    'exxoll'    : (0.0,     'solute LJ interactions -- old (exxnel)'),
    'exxolc'    : (0.0,     'solute coulombic interactions -- old (exxnec)'),

    'epol'      : (0.0,     'solvent polarizable energy'),
    'epol1'     : (0.0,     'solvent polarizable energy'),
    'epol2'     : (0.0,     'solvent polarizable energy'),
    'etpol'     : (None,    'list for solvent polarizable energy'),
    'etpol1'    : (None,    'list for solvent polarizable energy'),
    'etpol2'    : (None,    'list for solvent polarizable energy'),

    'ebndne'    : (0.0,     'solute bond stretching energy'),
    'ebnne1'    : (0.0,     'solute bond stretching energy'),
    'ebnne2'    : (0.0,     'solute bond stretching energy'),
    'ebndol'    : (0.0,     'solute bond stretching energy -- old (ebndne)'),
    'ebnol1'    : (0.0,     'solute bond stretching energy -- old (ebnne1)'),
    'ebnol2'    : (0.0,     'solute bond stretching energy -- old (ebnne2)'),
    'ebnd'      : (0.0,     'solute bond stretching energy -- update (ebndol)'),
    'ebnd1'     : (0.0,     'solute bond stretching energy -- update (ebnol1)'),
    'ebnd2'     : (0.0,     'solute bond stretching energy -- update (ebnol2)'),

    'esbnne'    : (0.0,     'solvent bond stretching energy'),
    'esanne'    : (0.0,     'solvent angle bending energy'),
    'esdine'    : (0.0,     'solvent dihedral torsional energy'),
    'esnbne'    : (0.0,     'solvent LJ and coulombic energy'),
    'esbnol'    : (0.0,     'solvent bond stretching energy -- old (esbnne)'),
    'esanol'    : (0.0,     'solvent angle bending energy -- old (esanne)'),
    'esdiol'    : (0.0,     'solvent dihedral torsional energy -- old (esdine)'),
    'esnbol'    : (0.0,     'solvent LJ and coulombic energy -- old (esnbne)'),
    'esinne'    : (0.0,     'sum = esbnne + esanne + esdine + esnbne'),

    'ebcnew'    : (0.0,     'solvent restrained energy'),
    'ebcne1'    : (0.0,     'solvent restrained energy'),
    'ebcne2'    : (0.0,     'solvent restrained energy'),
    'ebcold'    : (0.0,     'solvent restrained energy -- old (ebcnew)'),
    'ebcol1'    : (0.0,     'solvent restrained energy -- old (ebcne1)'),
    'ebcol2'    : (0.0,     'solvent restrained energy -- old (ebcne2)'),

    'eangne'    : (0.0,     'solute angle bending energy'),
    'eanne1'    : (0.0,     'solute angle bending energy'),
    'eanne2'    : (0.0,     'solute angle bending energy'),
    'eangol'    : (0.0,     'solute angle bending energy -- old (eangne)'),
    'eanol1'    : (0.0,     'solute angle bending energy -- old (eanne1)'),
    'eanol2'    : (0.0,     'solute angle bending energy -- old (eanne2)'),
    'eang'      : (0.0,     'solute angle bending energy -- update (eangol)'),
    'eang1'     : (0.0,     'solute angle bending energy -- update (eanol1)'),
    'eang2'     : (0.0,     'solute angle bending energy -- update (eanol2)'),

    'edihne'    : (0.0,     'solute dihedral torsional energy'),
    'edine1'    : (0.0,     'solute dihedral torsional energy'),
    'edine2'    : (0.0,     'solute dihedral torsional energy'),
    'edihol'    : (0.0,     'solute dihedral torsional energy -- old, (edihne)'),
    'ediol1'    : (0.0,     'solute dihedral torsional energy -- old, (edine1)'),
    'ediol2'    : (0.0,     'solute dihedral torsional energy -- old, (edine2)'),
    'edih'      : (0.0,     'solute dihedral energy -- update (edihol)'),
    'edih1'     : (0.0,     'solute dihedral energy -- update (ediol1)'),
    'edih2'     : (0.0,     'solute dihedral energy -- update (ediol2)'),

    'eslj'      : (0.0,     'solvent-solvent LJ energy'),
    'esco'      : (0.0,     'solvent-solvent coulombic energy'),
    'esljol'    : (0.0,     'solvent-solvent LJ energy -- old (eslj)'),
    'escool'    : (0.0,     'solvent-solvent coulombic energy -- old (escool)'),

    'enbne'     : (0.0,     'solute-solute QM energy'),
    'enbne1'    : (0.0,     'solute-solute QM energy'),
    'enbne2'    : (0.0,     'solute-solute QM energy'),
    'enbol'     : (0.0,     'solute-solute QM energy -- old (enbne)'),
    'enbol1'    : (0.0,     'solute-solute QM energy -- old (enbne1)'),
    'enbol2'    : (0.0,     'solute-solute QM energy -- old (enbne2)'),
}


# a list containts where to find default setting files, sequence low to high priority
PYBOSS_DATA_FILE_PATH = [
    '~/.config/PyBOSS/data',
    '~/PyBOSS/data',
]



