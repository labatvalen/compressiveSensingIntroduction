function Di = colonne(D,i)
    [M,N] = size(D)
    Di = zeros(M,1)
    for k = 1:M
        Di(k) = D(k,i)
    end
endfunction

function Di = ligne(D,i)
    [M,N] = size(D)
    Di = zeros(1,N)
    for k = 1:N
        Di(k) = D(i,k)
    end
endfunction

function DD = ajout_colonne_fin(D,Di)
    [M,N] = size(D)
    DD = zeros(M,N+1)
    for i = 1:M
        for j = 1:N
            DD(i,j) = D(i,j)
        end
    end
    for k = 1:M
        DD(k,N+1) = Di(k)
    end
endfunction

function [niter, residu_final,alpha] = orth_matching_pursuit(D,kmax,eps,x)
    [M,N] = size(D)
    alpha=zeros(N,1)
    
    m1 = choix_mk(D,x)
    P=[m1]
    
    z1 = (x'*colonne(D,m1))/norm(colonne(D,m1))^2
    
    Rk = x-z1*colonne(D,m1)
    phi = colonne(D,m1)
    
    k = 1
    while (k < kmax) & (norm(Rk) > eps)
        mk = choix_mk(D,Rk)
        P=[P,mk];
        phi = ajout_colonne_fin(phi,colonne(D,mk))
        zk = inv(phi'*phi)*phi'*x
        [Mp,Np] = size(P)
        for i = 1:Np
            alpha(P(i)) = zk(i)
        end
        Rk = x - phi*zk
        k = k+1
    end
    niter = k
    residu_final = norm(Rk);
endfunction

function D0=init_D0(I,J)
    D0 = zeros(I,J)
    for i = 1:J
        D0(1,i) = 0.25
        D0(2,i) = 0.25
        D0(3,i) = 0.25
        D0(4,i) = 0.25
    end
endfunction

function Ei = calcul_Ei(X,D,A,i,K)
    somme = 0
    for j = 1:K
        if ~(j == i) then
            somme = somme + colonne(D,j)*ligne(A,j)
        end
    end
    Ei = X - somme
endfunction


function bool = contains(X,valeur)
    [M,N] = size(X)
    bool = %F
    for i = 1:M
        for j = 1:N
            if X(i,j) == valeur
                bool = %T
            end
        end
    end
endfunction

function res = calcul_omega(wk)
    [l,c] = size(wk)
    res = zeros(max(wk),c)
    for i = 1:max(wk)
        for j = 1:c
            if contains(wk,i) then
                omega(i,j)=1
            end
        end
    end
endfunction

function Res = setColonne(X,i,C)
    
endfunction

function [A] = KSVD(X,K,L,D)
    [M,l] = size(X)
    eps=1e-5;
    k_max=4000;
    A = []
    for i = 1:l
        [niter, residu_final, alpha] = orth_matching_pursuit(D, k_max, eps, colonne(X,i))
        A = [A, alpha]
    end
    for i = 1:K
        wk = []
        for j = 1:l
            if ~(A(i,j) == 0) then
                wk = [wk, j]
            end
        end
        Ei = calcul_Ei(X,D,A,i,K)
        omega_i = calcul_omega(wk)
        Eri = Ei * omega_i
        [U Delta V] = svd(Eri)
        setColonne(D,i,colonne(U,1))
    end
endfunction




A = KSVD(X,3,100)




