## Importations
import numpy as np
import math
import cmath

## Constantes
D2 = np.array([[1,1,2,5,0,0,3,-2,1,2,2,2],[0,-1,-1,1,0,0,5,0,2,2,7,-1],[1,1,1,5,1,2,2,1,1,1,1,5],[1,5,2,2,5,0,-4,5,1,5,0,0],[0,2,2,1,1,0,0,0,0,4,-1,-2],[-1,2,2,2,-2,-3,-4,1,1,1,1,0]])
x1 = np.array([[4/3-math.sqrt(2)/2],[4/3+math.sqrt(2)/2],[2/3]])
x2=np.array([[-10],[-10],[1],[21],[0],[9]])
eps=1e-6
k_max=4000

## Fonctions

# Fonction colonne
# Renvoie la ieme colonne dune matrice
def colonne(D,i):
    i = i-1 # Les tableaux commencent a 0 en python
    taille = np.shape(D)
    M = taille[0] # Nombre de lignes de la matrice
    Di = np.zeros((M,1))
    for k in range (0,M):
        Di[k] = D[k,i]
    return Di

# Fonction choix_mk
# Renvoie la position de l'atome cherche
def choix_mk(D,Rk):
    taille = np.shape(D)
    M = taille[0] # Nombre de lignes de la matrice
    N = taille[1] # Nombre de colonnes de la matrice
    X = np.zeros((N,1))
    for j in range(0,N):
        col = colonne(D,j)
        colT = col.T
        X[j] = abs(np.dot(colT,Rk))/np.linalg.norm(col)
    print(X)
    for j in range(0,N):
        if (np.amax(X) == X[j]):
            mk = j
    return mk

# Fonction mp
# Algorithme du matching pursuit
def mp(D,kmax,eps,x0):
    taille = np.shape(D)
    M = taille[0] # Nombre de lignes de la matrice
    N = taille[1] # Nombre de colonnes de la matrice
    alpha = np.zeros((N,1))
    Rk = x0
    P=[]
    x = 0
    k = 0
    while (k < kmax) & (np.linalg.norm(Rk) > eps):
        mk = choix_mk(D,Rk)
        P = [P,mk]
        print('Les atomes selectionnes : ')
        print(P)
        a = np.dot(Rk.T,colonne(D,mk))
        b = a / (np.linalg.norm(colonne(D,mk))**2)
        c = b * colonne(D,mk)
        x = x + c
        print('x = ')
        print(x)
        alpha[mk] = alpha[mk] + b
        Rk = Rk - c
        print('Norme de Rk = ')
        print(np.linalg.norm(Rk))
        k = k+1
    niter = k
    residu_final = np.linalg.norm(Rk)
    print('x=')
    print(x)
    print("Nombre d'iterations :")
    print(niter)
    print('Residu final =  ')
    print(residu_final)
    print('La reprssentation parcimonieuse alpha =  ')
    print(alpha)

# Fonction ajout_colonne_fin
# Ajoute la ligne Di a la fin de la matrice D
def ajout_colonne_fin(D,Di):
    taille = np.shape(D)
    M = taille[0] # Nombre de lignes de la matrice
    N = taille[1] # Nombre de colonnes de la matrice
    DD = np.zeros((M,N+1))
    for i in range(0,M):
        for j in range (0,N):
            DD[i,j] = D[i,j]
    for k in range(0,M):
        DD[k,N] = Di[k]
    return DD


# Fonction omp
# Algorithme de l'othogonal matching pursuit
def omp(D,kmax,eps,x):
    taille = np.shape(D)
    M = taille[0] # Nombre de lignes de la matrice
    N = taille[1] # Nombre de colonnes de la matrice
    alpha=np.zeros((N,1))
    Rk = x
    P = []
    Dk = D
    k = 0
    while (k < kmax) & (np.linalg.norm(Rk) > eps):
        mk = choix_mk(D,Rk)
        P = [P,mk]
        Dk = ajout_colonne_fin(Dk,colonne(Dk,mk))
        zk = np.dot(np.linalg.pinv(Dk),x)
        print('zk=')
        print(zk)
        Rk = Rk - np.dot(Dk,zk)
        k = k+1
    niter = k
    residu_final = np.linalg.norm(Rk)
    alpha = zk
    print('x = ')
    print(x)
    print("Nombre d'iterations =  ")
    print(niter)
    print('residu final =  ')
    print(residu_final)
    print('La representation parcimonieuse alpha =  ')
    print(alpha)

## Affichage
print(omp(D2, k_max, eps, x2))