#!/bin/python3
import numpy as np
import scipy.special as sp

num_quad = 10

analytical = -0.123368 - 0.212261 * 1j
quadrature = np.loadtxt("quadrature/deg_" + str(num_quad) + ".gqd", delimiter=",")
quad_points = quadrature[:, 1]
quad_weights = quadrature[:, 2]

f = lambda x: -1j / 4 * sp.hankel1(0.0, 1 + x) * (1 - x) / 2
nu = lambda x, a, b, c, d: a * x**3 + b * x**2 + c * x + d
dnu = lambda x, a, b, c: 3 * a * x**2 + 2 * b * x + c

nu_tilde = -1
LHS = np.array(
    [
        [6 * nu_tilde, 2, 0, 0],
        [3 * nu_tilde**2, 2 * nu_tilde, 1, 0],
        [1, 1, 1, 1],
        [-1, 1, -1, 1],
    ],
    dtype=np.float64,
)
RHS = np.array([[0], [0], [1], [-1]], dtype=np.float64)
coeffs = np.linalg.solve(LHS, RHS)

a = coeffs[0].item()
b = coeffs[1].item()
c = coeffs[2].item()
d = coeffs[3].item()

standard_quad_int = 0
telles_quad_int = 0
for k in range(num_quad):
    z = quad_points[k]
    w = quad_weights[k]
    standard_quad_int += f(z) * w
    telles_quad_int += f(nu(z, a, b, c, d)) * dnu(z, a, b, c) * w

print(
    np.abs((standard_quad_int - analytical) / analytical),
    np.abs((telles_quad_int - analytical) / analytical),
)
