
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
# ===========================================================================================================================
# import package
# ======================================================================================================================

import numpy as np
import matplotlib.pyplot as plt

# ==============================================================================
# constant value and cof
# ==============================================================================

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

N = 10

L = 0.1

dx = L / N

m = (h * p) / (k * Ac)

resd = np.full((N), 0, dtype=float)

mm = np.full((N), 0, dtype=float)

s = np.full((N), 0, dtype=float)

T = np.full((N), T0, dtype=float)

X = np.linspace(0, 10, 10)

T[0] = Tb

Tinf = 15 + 273

dt = 0.25

Tini = np.full((N), T0, dtype=float)

cof1 = (k * dt) / (rho * cp * (dx ** 2))

cof2 = (h * p * dt) / (rho * cp * Ac)

cof3 = k / (2 * dx * h * p)


# =========================
# open file for saving data
# ==========================



# =========================
# build RHS matrix
# =========================


def res():
    for j in range(1, N - 1):
        resd[j] = cof1 * (
            (T[j + 1] - 2 * T[j] + T[j - 1])) - cof2 * (
                          T[j] - Tinf)

    return resd


# ================================
# update step  for T4
# ================================


def update():
    for kk in range(0, N - 1):
        ff = res()
        mm[kk] = T[kk] + ff[kk]

    mm[N - 1] = (Tinf + (4 * cof3 * mm[N - 2]) - (
            mm[N - 3] * cof3)) / (
                        1 + 3 * cof3)

    return mm


# =============
# Main loop
# =============
TT = np.zeros((2400, 10), dtype=float)
for i in range(0, 2400):
    T = update()
    TT[i, :] = T
    with open("FTCS.txt", "w") as FTCS:
        FT = ""

        FT += str(np.transpose(T - 273))
        FT += "\n"
        FTCS.write(FT)
        FTCS.close()

   # if T[1]-273==115.46901658 or T[1]-273>115.46901658  :
      # print(i)
    print(np.transpose(T - 273))

# ============
#  plot
# ============

for i in range(0, 2400):
    if i%120==0:
      plt.plot(X, TT[i, :]-273,  label="TEMP")
      plt.grid(True)

plt.show()

