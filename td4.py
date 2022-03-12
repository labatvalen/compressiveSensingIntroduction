import random as rand
import numpy as np
import math
import matplotlib.pyplot as plt
import cs_TD1


# EXERCICE 2

def phi_1(M, N):
    phi = np.zeros((M, N))

    for m in range(M):
        for n in range(N):
            phi[m, n] = rand.random()

    return phi


def phi_2(M, N, p):
    phi = np.zeros((M, N))

    for m in range(M):
        for n in range(N):
            if rand.random() < p:
                phi[m, n] = -1
            else:
                phi[m, n] = 1

    return phi


def phi_3(M, N, p):
    phi = np.zeros((M, N))

    for m in range(M):
        for n in range(N):
            if rand.random() < p:
                phi[m, n] = 0
            else:
                phi[m, n] = 1

    return phi


def phi_4(M, N):
    phi = np.zeros((M, N))

    for m in range(M):
        for n in range(N):
            phi[m, n] = np.random.normal(0, 1/math.sqrt(M), 1)

    return phi


def phi_5(M, N, p):
    phi = np.zeros((M, N))

    for m in range(M):
        for n in range(N):
            if rand.random() < p:
                phi[m, n] = 0
            else:
                phi[m, n] = rand.random()

    return phi


def mutual_consistency(phi, D):
    (M, N) = phi.shape
    mu = np.zeros((M, N))

    for i in range(M):
        for j in range(N):
            phi_i = phi[i, :]
            D_j = D[:, j]
            mu[i, j] = math.sqrt(N) * abs(np.vdot(phi_i, D_j)) / (np.linalg.norm(phi_i, 2) * np.linalg.norm(D_j, 2))

    return np.amax(mu)


M = [25, 30, 45, 50, 100, 150, 200, 250]
N = 500

psy = cs_TD1.transfoCosDis(N)
mu = np.zeros((5, len(M)))

for i in range(len(M)):
    mu[0, i] = mutual_consistency(phi_1(M[i], N), psy)
    mu[1, i] = mutual_consistency(phi_2(M[i], N, 0.5), psy)
    mu[2, i] = mutual_consistency(phi_3(M[i], N, 0.5), psy)
    mu[3, i] = mutual_consistency(phi_4(M[i], N), psy)
    mu[4, i] = mutual_consistency(phi_5(M[i], N, 0.5), psy)

print(mu)


N = 500
fe = 400
t = [i/fe for i in range(1,N+1)]
f_0 = 50
sinus = [math.sin(2 * math.pi * f_0 * T) for T in t]
f_p = 100
tporteuse = [i/fe for i in range(1, len(sinus) + 1)]
porteuse = [math.cos(2 * math.pi * f_p * Tporteuse) for Tporteuse in tporteuse]
x = [sinus[i] * porteuse[i] for i in range(len(sinus))]

y = np.zeros((5,5))
for u in range(len(mu)):
    y[u] = mu[u] * x

abscisses = [i for i in range(N)]
plt.plot(abscisses, y)
plt.show()

