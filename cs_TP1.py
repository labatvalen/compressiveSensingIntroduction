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
    t = 2 # On fixe t = 2 ici
    taille = np.shape(D)
    M = taille[0] # Nombre de lignes de la matrice
    N = taille[1] # Nombre de colonnes de la matrice (nombre d'atomes)
    X = np.zeros((N,1))
    mk = []
    for j in range(0,N):
        col = colonne(D,j)
        colT = col.T
        X[j] = abs(np.dot(colT,Rk))/np.linalg.norm(col) # Contribution des atomes
    S = t * np.linalg.norm(Rk)/math.sqrt(N)
    #print(X)
    for j in range(0,N):
        if (X[j] > S):
            mk.append(j)
    return mk

def ajout_colonne_fin(D,Di,M):
    N = 0
    if (len(D) == 0):
        N = 0
    else:
        taille = np.shape(D)
        N = taille[1] # Nombre de colonnes de la matrice

    DD = np.zeros((M,N+1))
    for i in range(0,M):
        for j in range (0,N):
            DD[i,j] = D[i,j]
    for k in range(0,M):
        DD[k,N] = Di[k]
    return DD


def stomp(D,kmax,eps,x):
    taille = np.shape(D)
    M = taille[0] # Nombre de lignes de la matrice
    N = taille[1] # Nombre de colonnes de la matrice
    alpha=np.zeros((N,1))
    Rk = x
    P = []
    Dk = D
    k = 0
    mat = []
    while (k < kmax) & (np.linalg.norm(Rk) > eps):
        # Etape 2 : selection
        mk = choix_mk(D,Rk)
        taille_mk = len(mk)
        # Etape 3 : maj des indices P
        for z in range(0,taille_mk):
            P.append(mk[z])
        # Etape 4 : constriction de la matrice
        for i in range(0,taille_mk):
            indice = mk[i]
            col = colonne(D,indice)
            mat = ajout_colonne_fin(mat,col,M)
        # Etape 5 : probleme d'optimisation par MCO
        inverse = np.linalg.pinv(np.dot(mat.T,mat))
        p = np.dot(inverse,mat.T)
        alphatemp = np.dot(p,x)
        cpt = 0
        for y in range(0,taille_mk):
            indice = mk[y]
            alpha[indice] = alphatemp[cpt]
            cpt = cpt + 1
        Rk = x - np.dot(D,alpha)
        P = []
        for j in range(0,len(alpha)):
            if (alpha[j] != 0):
                P.append(j)
        k = k+1
    niter = k
    residu_final = np.linalg.norm(Rk)
    print('x = ')
    print(x)
    print("Nombre d'iterations =  ")
    print(niter)
    print('residu final =  ')
    print(residu_final)
    print('La representation parcimonieuse alpha =  ')
    print(alpha)

## Affichage
print(stomp(D2, k_max, eps, x2))


