# -*- coding: utf-8 -*-
'''
Created on 2013-11-19

@author: Xin Zhang
'''
import math

from matplotlib import pyplot
from numpy import linalg
import numpy
from numpy.lib.function_base import trapz
from numpy.ma.core import exp


PI = 3.1415926
H_BAR = 1.0545887e-34
SPEED_OF_LIGHT = 2.99792458e8
ELECTRN_CHARGE = 1.6021892e-19
ELECTRN_MASS = 9.109534e-31
VACUUM_PERMITTIVITY = 8.85418e-12
KB = 1.3806488e-23
T = 300

def q2xy(Q):
    for y in range(0, 100):
        yy = float(y) / 100.0
        x = 1 - 0.47 * yy
        if Eg(x, yy) <= (H_BAR * 2 * PI * SPEED_OF_LIGHT / (Q * 1e-6)):
            break
    return [x, yy]

def playWithXY(a, b, c, d, x, y):
    return a * (1 - x) * y + b * x * y + c * (1 - x) * (1 - y) + d * x * (1 - y)

def Eg(x, y):
    E = 1.35 - 1.17 * y + 0.668 * (1 - x) - 0.069 * y * (1 - x) + 0.18 * y ** 2 + 0.03 * (1 - x) * y ** 2 + 0.758 * (1 - x) ** 2 - 0.322 * y * (1 - x) ** 2
    return ELECTRN_CHARGE * E


'''
注意：S1和S2的误差是问题的主要来源！不同的文献有不同的公式，可能最后还是一个经验型的东西。
'''
def s_1(x, y):
    av = -1.0 / 3.0 * (c11(x, y) + 2 * c12(x, y)) * degdp(x, y)
    return -2 * av * (1 - c12(x, y) / c11(x, y)) * strain(x, y)

def s_2(x, y):
    return -b(x, y) * (1 + 2 * c12(x, y) / c11(x, y)) * strain(x, y)

def c11(x, y):
    return 1e11 * playWithXY(11.8, 8.329, 14.12, 10.22, x, y) * 0.1

def c12(x, y):
    return 1e11 * playWithXY(5.38, 4.526, 6.253, 5.76, x, y) * 0.1

def degdp(x, y):
    return 1e-6 * ELECTRN_CHARGE / 1e5 * playWithXY(11.5, 10.0, 11.0, 8.5, x, y)

def b(x, y):
    return ELECTRN_CHARGE * playWithXY(-1.7, -1.8 , -1.5, -2.0, x, y)

def strain(x, y):
    return -(latticeConstant(x, y) - latticeConstant(1, 0)) / latticeConstant(1, 0)

def latticeConstant(x, y):
    return 1e-10 * playWithXY(5.6533, 6.0584 , 5.4512, 5.8688, x, y)

def delta(x, y):
    return ELECTRN_CHARGE * playWithXY(0.34 , 0.43, 0.10, 0.10, x, y)


h = 5.5e-9
H = 10e-9
wx = 0.8
wy = 0.8
[bx, by] = q2xy(1.25)
LdIII = 1e-12
k = 1
length = 0.05e-9
offset = 0.6
z = numpy.arange(-(H + h) / 2 , (H + h) / 2, length)
steps = z.size

'''
composition
'''
def compositionIII(Ld):
    r = z.copy()
    for i in range(0, z.size):
        r[i] = wx + (bx - wx) * (2 - math.erf((h / 2 - z[i]) / Ld / 2) - math.erf((h / 2 + z[i]) / Ld / 2)) / 2
    return r

def compositionV(Ld):
    r = z.copy()
    for i in range(0, z.size):
        r[i] = wy + (by - wy) * (2 - math.erf((h / 2 - z[i]) / Ld / 2) - math.erf((h / 2 + z[i]) / Ld / 2)) / 2
    return r

CIn = compositionIII(LdIII)
CAs = compositionV(LdIII * k)

'''
potential energy
'''
def vc(x, y):
    E = Eg(x, y)
    return offset * (E - E[E.argmin()]) - (1 - offset) * s_1(x, y) + E[E.argmin()]

def vhh(x, y):
    E = Eg(x, y)
    return (1 - offset) * (E - E[E.argmin()] - s_1(x, y)) + s_2(x, y)

def vlh(x, y):
    E = Eg(x, y)
    return (1 - offset) * (E - E[E.argmin()] - s_1(x, y)) + 0.5 * (s_2(x, y) + delta(x, y)) - 0.5 * numpy.sqrt(delta(x, y) ** 2 - 2 * s_2(x, y) * delta(x, y) + 9 * s_2(x, y) ** 2)

Vc = vc(CIn, CAs)
Vhh = vhh(CIn, CAs)
Vlh = vlh(CIn, CAs)

'''
effective mass
'''
def mc(x, y):
    return ELECTRN_MASS * playWithXY(0.0632 , 0.0213, 0.17, 0.077, x, y)

def mhhPerp(x, y):
    return ELECTRN_MASS * playWithXY(0.5 , 0.41, 0.54, 0.12, x, y)

def mlhPerp(x, y):
    return ELECTRN_MASS * playWithXY(0.088 , 0.024, 0.16, 0.12, x, y)
Mc = mc(CIn, CAs)
Mhh = mhhPerp(CIn, CAs)
Mlh = mlhPerp(CIn, CAs)


'''
energy level
'''
def levelFDM(V, m):
    A = numpy.zeros(steps)
    B = numpy.zeros(steps)
    C = numpy.zeros(steps)  # @UnusedVariable
    W = numpy.zeros(steps * steps).reshape(steps, steps)
    A[steps - 1] = H_BAR ** 2 / (2 * length ** 2 * m[steps - 1])
    B[0] = H_BAR ** 2 / (2 * length ** 2 * m[0])
    for i in range(0, steps - 1):
        A[i] = H_BAR ** 2 / (length ** 2 * (m[i] + m[i + 1]))
    for i in range(1, steps):
        B[i] = H_BAR ** 2 / (length ** 2 * (m[i] + m[i - 1]))
    C = A + B + V
    W[0, 0] = C[0]
    W[0, 1] = -A[0]
    W[steps - 1, steps - 1] = C[steps - 1]
    W[steps - 1, steps - 2] = -B[steps - 1]
    for i in range(1, steps - 1):
        W[i, i] = C[i]
        W[i, i + 1] = -A[i]
        W[i, i - 1] = -B[i]
    return numpy.sort(linalg.eigvalsh(W))

cLevels = levelFDM(Vc, Mc)
hhLevels = levelFDM(Vhh, Mhh)
lhLevels = levelFDM(Vlh, Mlh)


'''
Density Of State
'''
def roucBulk(E):
    Ec = (H_BAR ** 2 * PI ** 2) / (2 * mc(wx, wy) * h ** 2)
    Etemp = (E - Ec + numpy.absolute(E - Ec)) / 2
    return (2 * mc(wx, wy)) ** 1.5 / (2 * PI ** 2 * H_BAR ** 3) * numpy.sqrt(Etemp)

def rouhhBulk(E):
    Ehh = (H_BAR ** 2 * PI ** 2) / (2 * mhhPerp(wx, wy) * h ** 2)
    Etemp = (E - Ehh + numpy.absolute(E - Ehh)) / 2
    return (2 * mhhPerp(wx, wy)) ** 1.5 / (2 * PI ** 2 * H_BAR ** 3) * numpy.sqrt(Etemp)

def roulhBulk(E):
    Elh = (H_BAR ** 2 * PI ** 2) / (2 * mlhPerp(wx, wy) * h ** 2)
    Etemp = (E - Elh + numpy.absolute(E - Elh)) / 2
    return (2 * mlhPerp(wx, wy)) ** 1.5 / (2 * PI ** 2 * H_BAR ** 3) * numpy.sqrt(Etemp)

def roucQW(E):
    r = roucBulk(E)
    t = 0
    r0 = r[0]
    for i in range(0, r.size):
        if E[i] < cLevels[t]:
            r[i] = 0 
        elif E[i] >= cLevels[t] and E[i] < cLevels[t + 1]:
            r[i] = r0
        elif E[i] >= cLevels[t] and E[i] >= cLevels[t + 1]:
            r0 = r[i]
            if t < cLevels.size - 2:
                t = t + 1
    return r

def rouhhQW(E):
    r = rouhhBulk(E)
    t = 0
    r0 = r[0]
    for i in range(0, r.size):
        if E[i] < hhLevels[t]:
            r[i] = 0 
        elif E[i] >= hhLevels[t] and E[i] < hhLevels[t + 1]:
            r[i] = r0
        elif E[i] >= hhLevels[t] and E[i] >= hhLevels[t + 1]:
            r0 = r[i]
            if t < hhLevels.size - 2:
                t = t + 1
    return r

def roulhQW(E):
    r = roulhBulk(E)
    t = 0
    r0 = r[0]
    for i in range(0, r.size):
        if E[i] < lhLevels[t]:
            r[i] = 0 
        elif E[i] >= lhLevels[t] and E[i] < lhLevels[t + 1]:
            r[i] = r0
        elif E[i] >= lhLevels[t] and E[i] >= lhLevels[t + 1]:
            r0 = r[i]
            if t < lhLevels.size - 2:
                t = t + 1
    return r

Etmp = numpy.arange(0 , Eg(wx, wy) * 10, Eg(wx, wy) / 1000)
pyplot.figure()
pyplot.plot(Etmp, roucBulk(Etmp), ':b', label='roucBulk')
pyplot.plot(Etmp, roucQW(Etmp), '-b', label='roucQW')
pyplot.plot(Etmp, rouhhBulk(Etmp), ':g', label='rouhh')
pyplot.plot(Etmp, rouhhQW(Etmp), '-g', label='rouhhQW')
pyplot.plot(Etmp, roulhBulk(Etmp), ':r', label='roulh')
pyplot.plot(Etmp, roulhQW(Etmp), '-r', label='roulhQW')
pyplot.xlabel('$E$')
pyplot.ylabel('$rou$')
pyplot.legend()


'''
Quasi-Fermi Levels and fermi inversion factor
'''
delta_n = 2e17
EfcminusEc = -Eg / 2
while 1:
    nc = trapz(roucQW(Etmp) / (1 + exp((Etmp - EfcminusEc) / KB / T)));
    if nc > delta_n:
        break
    EfcminusEc = EfcminusEc + Eg / 1000;
EhhminusEfv = Eg / 2;
while 1:
    nhh = trapz(rouhhQW(Etmp) / (1 + exp((Etmp - EhhminusEfv) / KB / T)));
    if nhh < delta_n:
        break
    EhhminusEfv = EhhminusEfv - Eg / 1000;
ElhminusEfv = Eg / 2;
while 1:
    nlh = trapz(roulhQW(Etmp) / (1 + exp((Etmp - ElhminusEfv) / KB / T)));
    if nlh < delta_n:
        break;
    ElhminusEfv = ElhminusEfv - Eg / 1000;

step = Eg(wx, wy) / 100
E = numpy.arange(Eg(wx, wy), 1.1 * max([(EfcminusEc + EhhminusEfv + Eg) (EfcminusEc + ElhminusEfv + Eg)]), step)

E2minusEc = (mrhh / mc(wx, wy)) * (E - Eg)
fc = 1 / (1 + exp((E2minusEc - EfcminusEc) / KB / T))
EhhminusE1 = (mrhh / mhhPerp(wx, wy)) * (E - Eg)
fv = 1 / (1 + exp((EhhminusEfv - EhhminusE1) / KB / T))
fghh = fc - fv

E2minusEc = (mrlh / mc(wx, wy)) * (E - Eg)
fc = 1 / (1 + exp((E2minusEc - EfcminusEc) / KB / T))
ElhminusE1 = (mrlh / mlhPerp(wx, wy)) * (E - Eg)
fv = 1 / (1 + exp((ElhminusEfv - ElhminusE1) / KB / T))
fglh = fc - fv;

'''
optical matrix element
'''
OME = (ELECTRN_MASS / mc(wx, wy) - 1) * ELECTRN_MASS * Eg / 2



pyplot.show()
