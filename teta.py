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
N = 20

cp = 0.91 * (10 ** 3)

rho = 2700

Tb = 120 + 273

k = 237

alfa = 97.1 * (10 ** -6)

h = 35

T0 = 30 + 273

D = 1 * 10 ** -2

Ac = np.pi * ((D ** 2) / 4)

p = np.pi * D

L = 0.1

dx = L / (N-1)

teta = 0.9

m = (h * p) / (k * Ac)

A = np.zeros((N-2, N-2), dtype=float)

B = np.zeros((N-2, 1), dtype=float)

T = np.full((N, 1), T0, dtype=float)

T[0] = Tb

Tinf = 15 + 273

dt = 0.01

Tini = np.full((3, 1), T0, dtype=float)

cof1 = (k * dt * teta) / ((dx ** 2) * (rho * cp))

cof2 = k / h

cof3 = (h * p * teta * dt) / (rho * cp * Ac)

cof4 = (k * dt * (1 - teta)) / (rho * cp * (dx ** 2))

cof5 = k / (2 * dx * h * p)

cof6 = (h * p * (1 - teta) * dt) / (rho * cp * Ac)

cof7 = (h * p * dt) / (rho * cp * Ac)


# ===========================
# build cof matrix
# ===========================

def cof():
    for i in range(0, N-2):
        if i == 0:
            A[i, i] = 1 + (2 * cof1) + cof3

            A[i, i+1] = -1 * cof1
        elif i == N-3:

            A[i, i-1] = (-1 * cof1 * (
                        -1 * cof5 / (1 + 3 * cof5))) - cof1

            A[i, i] = (-1 * cof1 * ((+4 * cof5) / (1 + 3 * cof5))) + \
                      (1 + (2 * cof1) + cof3)
        else:
            A[i, i-1] = -1 * cof1

            A[i, i] = 1 + (2 * cof1) + cof3

            A[i, i+1] = -1 * cof1

    A[2, 0] = 0

    return A


# =========================
# build RHS matrix
# =========================

def rhs():
    for i in range(0, N-2):
        if i == 0:
            B[i, 0] = (T[i] * cof4) + (T[i+2] * cof4) + (
                (T[i+1] * (1 - (2 * cof4) - cof6)) + (
                    cof7 * Tinf)) + (cof1 * T[i])
        elif i == N-3:
            B[i, 0] = (T[i] * cof4) + (T[i+2] * cof4) + (T[i+1] * (
                    1 - (2 * cof4) - cof6)) + (cof7 * Tinf)
        else:
            B[i, 0] = (T[i] * cof4) + (T[i+2] * cof4) + (
                    T[i+1] * (1 - (2 * cof4) - cof6)) + (
                              cof7 * Tinf)



    return B


# ================================
# update step  for T4
# ================================


def update():
    for i in range(0, N-2):
        T[i+1] = ss[0][i]


    T[N-1] = (Tinf + (4 * cof5 * (ss[0][N-3])) - ((ss[0][N-4]) * cof5)) / (
            1 + 3 * cof5)
    return T


# =============
# Main loop
# =============

for i in range(0, 10000):

    if i == 0:

        teta = 0.9

        A = cof()

        B = rhs()

        ss = spla.gmres(A, B, tol=1e-5)

        T = update()

    else:
        teta = 0.9

        A = cof()

        B = rhs()

        ss = spla.gmres(A, B, tol=1e-5)

        T = update()



    print(np.transpose(T - 273))

# ============
#  plot
# ============
X = np.linspace(0, 10,N)
p1 = plt.plot(X, T-273,label="teta = 0.9" )

plt.legend(loc=' left', bbox_to_anchor=(1, 1))
plt.grid(True)

plt.show()

