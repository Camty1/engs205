#!/bin/python3
import numpy as np
import scipy as sp
import cmath

x = np.arange(0, 2.05, 0.1)
int_H0 = np.zeros(x.shape)

int_H0 = x * sp.special.hankel1(0, x) + np.pi / 2 * x * (
    sp.special.struve(0, x) * sp.special.hankel1(1, x)
    - sp.special.struve(1, x) * sp.special.hankel1(0, x)
)

print(int_H0)
