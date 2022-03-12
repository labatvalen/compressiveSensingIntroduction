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

def ajout_colonne_fin(D,Di,M):
    N = 0
    if (D == []):
        print('non')
        N = 0
    else:
        print('oui')
        taille = np.shape(D)
        N = taille[1] # Nombre de colonnes de la matrice

    DD = np.zeros((M,N+1))
    for i in range(0,M):
        for j in range (0,N):
            DD[i,j] = D[i,j]
    for k in range(0,M):
        DD[k,N] = Di[k]
    return DD

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

def ksvd(X,K):
    

print ksvd(x1,D2)