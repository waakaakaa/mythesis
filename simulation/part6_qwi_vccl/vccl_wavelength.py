# -*- coding: utf-8 -*-
'''
Created on 2013-11-19

@author: Xin Zhang
'''
import cmath
import math

from matplotlib import pyplot
import numpy
from numpy.lib.polynomial import roots
from scipy.optimize.minpack import fsolve


C11 = 0.755
C12 = C11 - 1
C22 = C11
C21 = C12

n1 = 3.215
n2 = 3.215

r1 = (n1 - 1) / (n1 + 1)
r2 = (n1 - 1) / (n1 + 1)

c = 2.99792458e8

delta_f1 = 100e9
L1 = 1e6 * c / (2 * n1 * delta_f1)
f0 = 193400
m1 = 1934
lambda0 = c / f0 / 1000
gain0 = -10000 * math.log(r1 * r2) / L1
m2 = 2150
delta_f2 = f0 * 1e9 / m2
L2 = 1e6 * c / (2 * n2 * delta_f2)

e = 1.6e-19  # 电子的电量为1e-19 C
A = 1 / (2.776 * 1e-9)  # nonradiative coefficient
B = 1 * 2e-10 * 1e-6  # bimolecular recombination coefficient
C = 1 * 8e-29 * 1e-12  # Auger recombination coefficient
dn_dN = -5.97e-27  # index derivative wrt carrier density(m3)
dalpha_dN = 2.56e-21  # abosorption derivative wrt carrier density(m2)
tau_p = 0.08  # confinement factor
L3 = L2 / 3 * 2  # length of dbr section(172.7708μm)
w2 = 3e-6  # 波导层宽度（3μm）
d2 = 114e-9  # 波导层厚度（114nm）
V3 = L3 * w2 * d2 * 1e-6  # 波长切换区体积

def VCCL_threshold_dynamic(x):
    k = 2 * math.pi / x[0]
    g1 = x[1] / 20000
    g2 = g1 * L1 / L2
    y1 = r1 * r2 * cmath.exp(2 * g1 * L1 + 2 * 1j * k * n1 * L1)
    y2 = r1 * r2 * cmath.exp(2 * g2 * L2 + 2 * 1j * (k * delta_n3 * L3 + k * n2 * L2))
    yc = C11 * y1 + C22 * y2 - (C11 * C22 - C12 * C21) * y1 * y2 - 1
    yr = yc.real
    yi = yc.imag
    return [yr, yi]

I3 = numpy.arange(0, 0.8, 0.001)  # 波长转换区电流(A)
lambda1 = I3.copy()
for ii in range(0, I3.size):
    a = C
    b = B
    c1 = A
    d = -I3[ii] / e / V3
    D = [a, b, c1, d]
    NN3 = roots(D)
    for r in [0, 1, 2]:
        if (NN3[r].imag == 0) and (NN3[r] >= 0):
            N3 = NN3[r]
    delta_n3 = tau_p * dn_dN * N3  # the carrier-induced index change
 
    # Calculate threshold gains of various modes
    wl = numpy.arange(1, 10, 1)
    threshold_gain = numpy.arange(1, 10, 1)
    for m in range(m1 - 4, m1 + 5):
        lmd = 1e6 * c / (m * delta_f1)
        X1 = fsolve(VCCL_threshold_dynamic, [lmd, 0.95 * gain0])
        wl[m - m1 + 4] = X1[0] * 1000
        threshold_gain[m - m1 + 4] = X1[1]
    lambda1[ii] = wl[threshold_gain.argmin()]
pyplot.figure()
pyplot.plot(I3[0:I3.size-2], lambda1[0:I3.size-2], ':k+')
pyplot.xlabel('$Wavelength switching section current (A)$')
pyplot.ylabel('$laser wavelength (nm)$')
pyplot.show()
