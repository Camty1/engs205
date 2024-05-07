#!/bin/python3
import numpy as np
import scipy as sp
import cmath

k = 13.84
ds = 0.130526
x = np.array([k * ds])
int_J0 = np.zeros(x.shape)
int_H0 = np.zeros(x.shape)
struve_0 = np.zeros(x.shape)
struve_1 = np.zeros(x.shape)
max_k = 20


int_H0 = x * sp.special.hankel1(0, x) + np.pi / 2 * x * (
    sp.special.struve(0, x) * sp.special.hankel1(1, x)
    - sp.special.struve(1, x) * sp.special.hankel1(0, x)
)

g1 = 1j / (4 * k) * (int_H0 - sp.special.hankel1(1, x))
g2 = 1j / (4 * k) * sp.special.hankel1(1, x)

print(g1, g2)
