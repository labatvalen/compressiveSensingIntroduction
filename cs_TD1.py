## Importations
import numpy as np
import math
import cmath

## Transformee en cosinus discrete

# Constantes
N = 5

# Fonction
def transfoCosDis(N):
    psy = np.zeros((N,N))
    for k in range(N):
        for l in range(N):
            if k == 0:
                psy[k,l] = math.sqrt(1./N)
            else:
                psy[k,l] = math.sqrt(2./N) * math.cos((math.pi*k*(1./2.*N))/N)
    return psy

def printTransfoFouDis():
    psy = transfoCosDis(N)
    print('Transformee en cosinus disrete :')
    print(psy)
    print(psy * np.transpose(psy))


## Transformee de Fourier discrete

# Constantes
N = 5

# Fonction
def transfoFouDis(N):
    c = []
    psy = np.zeros((N,N))
    for k in range(1,N+1):
        b = []
        for l in range(1,N+1):
            a = complex((1./math.sqrt(N)),0.) * cmath.exp(complex(0,1.) * complex( 2. * (k-1.) * (l-1.) / N, 0.) )
            b.append(a)
        c.append(b)
    return c

def printTransfoFouDis():
    psy = transfoFouDis(N)
    print('Transformee de Fourier disrete :')
    print(psy)