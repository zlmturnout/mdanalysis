# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
#

"""
MDAnalysis topology tables
==========================

The module contains static lookup tables for atom typing etc. The
tables are dictionaries that are indexed by the element.

.. autodata:: atomelements
.. autodata:: masses

The original raw data are stored as multi-line strings that are
translated into dictionaries with :func:`kv2dict`. In the future,
these tables might be moved into external data files; see
:func:`kv2dict` for explanation of the file format.

.. autofunction:: kv2dict

The raw tables are stored in the strings

.. autodata:: TABLE_ATOMELEMENTS
.. autodata:: TABLE_MASSES
"""

def kv2dict(s, convertor=str):
    """Primitive ad-hoc parser of a key-value record list.

    * The string *s* should contain each key-value pair on a separate
      line (separated by newline). The first white space after the key
      separates key and value.

    * Empty lines are allowed.

    * Comment lines (starting with #) are allowed.

    * Leading whitespace is ignored.

    The *convertor* is a function that converts its single argument to
    a valid Python type. The default is :func:`str` but other
    possibilities are :func:`int` (for integers) or :func:`float` for
    floating point numbers.
    """
    d = {}
    lines = s.splitlines()
    for line in lines:
        line = line.lstrip()
        values = line.split(None,1)
        if len(values) == 0 or line.startswith("#"):
            continue
        d[values[0]] = convertor(values[1])
    return d

#: Table with hard-coded special atom names, used for guessing atom types
#: with :func:`MDAnalysis.topology.core.guess_atom_element`.
TABLE_ATOMELEMENTS = """
# translation of atomnames to types/element
# based on CHARMM and AMBER usage with a little bit of GROMOS
# NOTE: CL might be ambiguous and is interpreted as chloride!

# --------- ------------------
# atomname   element
# --------- ------------------

# Calcium
CAL          CA
C0           CA
CA2+         CA

# Cesium
CES          CS

# Chloride
CLA          CL
CLAL         CL
CL           CL
CL-          CL

# Iron
FE           FE

# Lithium
LIT          LI
LI           LI
LI+          LI
QL           LI

# Magnesium
MG           MG
MG2+         MG

# Noble gases
## XXX collides with NE, HE in Arg  XXX
## XXX so we remove the noble gases XXX
##HE           HE
##NE           NE

# Potassium
POT          K
K+           K
QK           K

# Sodium
SOD          NA
NA           NA
NA+          NA
QN           NA

# Zink
ZN           ZN

# Copper
CU           CU

# Cesium
CS           CS
CS+          CS
CES          CS

# Cerium??
QC           CE

# Rubidium
RB           RB
QR           RB

# special carbons (Amber?)
BC           C
AC           C

AU            AU

# other types are guessed from the name; see
# topology.core.guess_atom_elements()
"""

#: Dictionary with hard-coded special atom names, used for guessing atom types
#: with :func:`MDAnalysis.topology.core.guess_atom_type`.
atomelements = kv2dict(TABLE_ATOMELEMENTS)

#: Plain-text table with atomic masses in u.
TABLE_MASSES = """
# masses for elements in atomic units (u)
# (taken from CHARMM and Gromacs atommass.dat)
#------------ -----------
# atomtype    mass
#------------ -----------
AC    227.028
AG    107.8682
AL    26.981539
AM    243
AR    39.948
AS    74.92159
AT    210
AU    196.96654
B     10.811
BA    137.327
BE    9.012182
BH    262
BI    208.98037
BK    247
BR    79.90400
C     12.01100
CA    40.08000
CD    112.411
CE    140.11600
CF    251
CL    35.45000
CM    247
CO    58.9332
CR    51.9961
CS    132.90000
CU    63.54600
DB    262
DY    162.5
ER    167.26
ES    252
EU    151.965
F     18.99800
FE    55.84700
FM    257
FR    223
GA    69.723
GD    157.25
GE    72.61
H     1.00800
HE    4.00260
HF    178.49
HG    200.59
HO    164.93032
HS    265
I     126.90450
IN    114.82
IR    192.22
K     39.10200
KR    83.8
LA    138.9055
LI    6.941
LR    262
LU    174.967
MD    258
MG    24.30500
MN    54.93805
MO    95.94
MT    266
N     14.00700
NA    22.989768
NB    92.90638
ND    144.24
NE    20.17970
NI    58.6934
NO    259
NP    237.048
O     15.99900
OS    190.2
P     30.97400
PA    231.0359
PB    207.2
PD    106.42
PM    145
PO    209
PR    140.90765
PT    195.08
PU    244
RA    226.025
RB    85.46780
RE    186.207
RF    261
RH    102.9055
RN    222
RU    101.07
S     32.06000
SB    121.757
SC    44.95591
SE    78.96
SG    263
SI    28.0855
SM    150.36
SN    118.71
SR    87.62
TA    180.9479
TB    158.92534
TC    98
TE    127.6
TH    232.0381
TI    47.88
TL    204.3833
TM    168.93421
U     238.0289
V     50.9415
W     183.85
XE    131.29
Y     88.90585
YB    173.04
ZN    65.37000
ZR    91.224
"""

#: Dictionary table with atomic masses in u, indexed by the element from
#: :data:`atomelements`.
masses = kv2dict(TABLE_MASSES, convertor=float)
