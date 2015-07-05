# -*- coding: utf-8 -*-
'''
Created on 2013-11-18

@author: Xin Zhang
'''

import cmath
import math

from matplotlib import pyplot
import numpy
from scipy.constants.constants import pi
from scipy.optimize.minpack import fsolve


C11 = cmath.sqrt(0.904)
C12 = cmath.sqrt(1 - C11 * C11)
C22 = cmath.sqrt(0.904)
C21 = cmath.sqrt(1 - C22 * C22)

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

def VCCL_threshold_static(x):
    k = 2 * pi / x[0]
    g1 = x[1] / 20000
    g2 = g1 * L1 / L2
    y1 = r1 * r2 * cmath.exp(2 * g1 * L1 + 2 * 1j * k * n1 * L1)
    y2 = r1 * r2 * cmath.exp(2 * g2 * L2 + 2 * 1j * k * n2 * L2)
    yc = C11 * y1 + C22 * y2 - (C11 * C22 - C12 * C21) * y1 * y2 - 1
    yr = yc.real
    yi = yc.imag
    return [yr, yi]

X = fsolve(VCCL_threshold_static, [lambda0, gain0])

# Calculate the effective "reflectivity" of the second cavity as seen by the first cavity, and vice versa
wavelength = numpy.arange(1530, 1570, 0.0005)
k = 2000 * pi / wavelength
g1 = X[1] / 20000
g2 = g1 * L1 / L2
y2 = numpy.exp(2 * g2 * L2) * numpy.exp(2 * 1j * k * n2 * L2)
eta2 = C11 + C12 * C21 * r1 * r2 * y2 / (1 - C22 * r1 * r2 * y2);
ETA2 = (numpy.absolute(eta2)) ** 2;
y1 = numpy.exp(2 * g1 * L1) * numpy.exp(2 * 1j * k * n1 * L1);
eta1 = C22 + C12 * C21 * r1 * r2 * y1 / (1 - C11 * r1 * r2 * y1);
ETA1 = (numpy.absolute(eta1)) ** 2;

pyplot.figure()
pyplot.plot(wavelength, ETA1, ':k', label='ETA1')
pyplot.plot(wavelength, ETA2, '-k', label='ETA2')
pyplot.xlabel('$Wavelength (nm)$');
pyplot.ylabel('$Threshold Gain (cm^{-1})$');
pyplot.legend()

# Calculate threshold gains of various modes
wl = numpy.zeros(102).reshape(2, 51)
threshold_gain = numpy.zeros(102).reshape(2, 51)
for m in range(m1 - 25, m1 + 25):
    lmd = 1e6 * c / (m * delta_f1)
    X1 = fsolve(VCCL_threshold_static, [lmd, 0.8 * gain0])
    wl[0, m - m1 + 25] = X1[0] * 1000
    threshold_gain[0, m - m1 + 25] = X1[1]
for m in range(m2 - 25, m2 + 25):
    lmd = 1e6 * c / (m * delta_f2)
    X1 = fsolve(VCCL_threshold_static, [lmd, 0.8 * gain0])
    wl[1, m - m2 + 25] = X1[0] * 1000
    threshold_gain[1, m - m2 + 25] = X1[1]
pyplot.figure()
pyplot.plot(wl[0, 0:50], threshold_gain[0, 0:50], 'ko')
pyplot.plot(wl[1, 0:50], threshold_gain[1, 0:50], ':ko')
pyplot.xlabel('$Wavelength (nm)$');
pyplot.ylabel('$|\eta|^2$');

# Calculate the threshold difference of adjacent modes as a function of coupling coefficients
lmd = 1e6 * c / ((m1 + 1) * delta_f1)
coupling = numpy.zeros(1001).reshape(1, 1001)
threshold_diff = numpy.zeros(3003).reshape(3, 1001)
wl_diff = numpy.zeros(3003).reshape(3, 1001)
for p in range(0, 1000):
    coupling[0, p] = p * 0.001
    C12 = numpy.sqrt(coupling[0, p])
    C11 = numpy.sqrt(1 - coupling[0, p])
    C21 = numpy.sqrt(coupling[0, p])
    C22 = numpy.sqrt(1 - coupling[0, p])
    X = fsolve(VCCL_threshold_static, [lambda0, 0.75 * gain0])
    X1 = fsolve(VCCL_threshold_static, [lmd, 0.75 * gain0])
    threshold_diff[0, p] = X[1]
    threshold_diff[1, p] = X1[1]
    threshold_diff[2, p] = X1[1] - X[1]
    wl_diff[0, p] = X[0] * 1000
    wl_diff[1, p] = X1[0] * 1000
    wl_diff[2, p] = (X1[0] - X[0]) * 1000
fig = pyplot.figure()
ax = fig.add_subplot(211)
ax.plot(coupling[0, 0:1000], threshold_diff[0, 0:1000], '-k')
ax.plot(coupling[0, 0:1000], threshold_diff[1, 0:1000], ':k')
pyplot.xlabel('$Cross-Coupling Coefficient$')
pyplot.ylabel('$Threshold Gain (cm^{-1})$')
ax = fig.add_subplot(212)
ax.plot(coupling[0, 0:1000], threshold_diff[2, 0:1000], '-k')
pyplot.xlabel('$Cross-Coupling Coefficient$')
pyplot.ylabel('$Threshold Difference (cm^{-1})$')
pyplot.show()
