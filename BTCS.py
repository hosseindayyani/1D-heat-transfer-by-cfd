"""
solving the problem for bar with two boundary condition
first we should build our matrix

Tb  =   base temp(c)

k   =   conduction cof(W/m.k)

alfa=   thermal diffusivity(m^2/s)

Ac  =   cross-section area(cm^2)

L   =   len of fin(cm)

Tinf=   ambient temp(c)

alfa=   thermal conductivity

h   =   convection cof(W/m^2K)

T0  =   initial fin temp(c)

R   =   diameter



"""
# ---------------#
# import package#
# ---------------#

import numpy as np

import scipy.sparse.linalg as spla

import matplotlib.pyplot as plt

# ========================
# constant value and cof
# ========================
N = 50

cp = 0.91 * (10 ** 3)

rho = 2700

Tb = 120 + 273


alfa = 97.1 * (10 ** -6)

h = 35

T0 = 30 + 273

D = 1 * 10 ** -2

Ac = np.pi * ((D ** 2) / 4)

p = np.pi * D

L = 0.1

dx = L / (N-1)

teta = 1

A = np.zeros((N-2, N-2), dtype=float)

B = np.zeros((N-2, 1), dtype=float)

T = np.full((N, 1), T0, dtype=float)

T[0] = Tb

Tinf = 15 + 273

dt = 0.25

Tini = np.full((3, 1), T0, dtype=float)

k = np.zeros((N, 1), dtype=float)

cof8 = np.zeros((N, 1), dtype=float)

def  coeff():

    alfaperim = 1.829 * (10**-4)

    beta = 0.0245

    k = T/(((1.829 * (10**-4)) * (T**2)) + beta)

    k1 = (((alfaperim * (T**2 + beta))-(2 * alfaperim * (T**2)))/((
            (alfaperim * (T**2)) + beta)**2))

    cof1 = (k * dt * teta) / ((dx ** 2) * (rho * cp))

    cof2 = k / h

    cof3 = (h * p * teta * dt) / (rho * cp * Ac)

    cof4 = (k * dt * (1 - teta)) / (rho * cp * (dx ** 2))

    cof5 = k / (2 * dx * h * p)

    cof6 = (h * p * (1-teta) * dt) / (rho * cp * Ac)

    cof7 = (h * p * dt) / (rho * cp * Ac)

    cof8[1:N-1] = ((k1[1:N-1] * (T[2:N] - T[0:N-2])) * dt) / (4 * (
            dx**2) * rho * cp)

    cof8[0] = ((k1[0] * (T[1] - T[0])) * dt)/(2 * (dx**2) * rho * cp)

    cof8[N-1] = ((k1[N-1] * (T[N-1] - T[N-2])) * dt)/(2 * (
            dx**2) * rho * cp)

    return cof1, cof2, cof3, cof4, cof5, cof6, cof7, cof8


# ===========================
# build cof matrix
# ===========================

def cof():
    cf = coeff()
    cof1 = cf[0]
    cof3 = cf[2]
    cof5 = cf[4]
    cof8 = cf[7]
    for i in range(0, N-2):
        if i == 0:
            A[i, i] = 1 + (2 * cof1[i+1][0]) + cof3

            A[i, i+1] = -1 * cof1[i+1] - cof8[i+1]
        elif i == N-3:

            A[i, i-1] = ((-1 * cof1[i+1]-(1 * cof8[i+1])) * (
                    -1 * cof5[i+1] / (1 + 3 * cof5[i+1]))) - cof1[i+1]

            A[i, i] = ((-1 * cof1[i+1]-(1 * cof8[i+1])) * ((+4 * cof5[i+1]) / (
                1 + 3 * cof5[i+1]))) + (1 + (2 * cof1[i+1]) + cof3)

        else:
            A[i, i-1] = -1 * cof1[i+1] + cof8[i+1]

            A[i, i] = 1 + (2 * cof1[i+1]) + cof3

            A[i, i+1] = -1 * cof1[i+1] - cof8[i+1]

    A[2, 0] = 0

    return A


# =========================
# build RHS matrix
# =========================

def rhs():
    cf = coeff()
    cof4 = cf[3]
    cof6 = cf[5]
    cof1 = cf[0]
    cof7 = cf[6]
    for i in range(0, N-2):
        if i == 0:
            B[i, 0] = (T[i] * cof4[i+1]) + (T[i+2] * cof4[i+1]) + (
                (T[i+1] * (1 - (2 * cof4[i+1]) - cof6)) + (
                    cof7 * Tinf)) + (cof1[i+1] * T[i]) - (
                        cof8[i+1] * T[i])

        elif i == N-3:
            B[i, 0] = (T[i] * cof4[i+1]) + (T[i+2] * cof4[i+1]) + (T[i + 1] * (
                    1 - (2 * cof4[i+1]) - cof6)) + (cof7 * Tinf)
        else:
            B[i, 0] = (T[i] * cof4[i+1]) + (T[i+1] * cof4[i+1]) + (
                    T[i+1] * (1 - (2 * cof4[i+1]) - cof6)) + (
                              cof7 * Tinf)

    return B


# ================================
# update step  for T4
# ================================


def update():
    for i in range(0, N-2):
        T[i+1] = ss[0][i]

    cf = coeff()
    cof5 = cf[4]

    T[N-1] = (Tinf + (4 * cof5[N-1] * (ss[0][N-3])) - ((
        ss[0][N-4]) * cof5[N-1])) / (
            1 + 3 * cof5[N-1])
    return T


# =============
# Main loop
# =============

for i in range(0, 50000):

    if i == 0:

        teta = 1

        A = cof()

        B = rhs()

        ss = spla.gmres(A, B, tol=1e-5)

        T = update()

    else:
        teta = 1

        A = cof()

        B = rhs()

        ss = spla.gmres(A, B, tol=1e-5)

        T = update()

    print(np.transpose(T-273))


# ============
#  plot
# ============
X = np.linspace(0, 10,N)
plt.plot(X, T-273, 'b', label="TEMP")
plt.grid(True)

plt.show()

