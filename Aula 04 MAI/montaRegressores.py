import numpy as np


def montaRegressores(y, u, ny, nu, ne, L, N, qs, qe):
    """
    y  = vetor linha de dados de saída
    u  = vetor linha de dados de entrada
    ny = número de regressores de saída
    nu = número de regressores de entrada
    ne = número de regressores do ruído
    L  = grau de não-linearidade
    N  = número de observações
    qs = quantidade de entradas
    qe = quantidade de saídas
    """
    Py = np.zeros(ny)
    if nu != 0:
        Pu = np.zeros(nu)
    n = ny + nu

    # Quantidade de regressores por grau
    ni = np.zeros(L + 1, dtype=int)
    ni[0] = 1
    for i in range(1, L + 1):
        ni[i] = (ni[i - 1] * (qs * ny + qe * nu + ne + i)) // i

    # Vetores com condições iniciais nulas
    yP = np.concatenate([Py, y])
    if nu != 0:
        uP = np.concatenate([Pu, u])

    PL = np.zeros((N, n))
    P = np.zeros((N, ni[L] - 1))
    Y = np.zeros((N, 1))

    for k in range(N):
        # Geração matriz linear
        if nu != 0:
            PL[k, :] = np.concatenate([
                yP[ny + k - 1:k - 1:-1] if k > 0 else yP[ny + k - 1::-1][:ny],
                uP[nu + k - 1:k - 1:-1] if k > 0 else uP[nu + k - 1::-1][:nu]
            ])
        else:
            PL[k, :] = yP[ny + k - 1:k - 1:-1] if k > 0 else yP[ny + k - 1::-1][:ny]

        # Geração matriz não-linear
        indice1 = np.ones(n, dtype=int)
        indice2 = np.copy(indice1)
        P1 = PL[k, :].copy()
        P[k, :len(PL[k, :])] = PL[k, :]

        for l in range(2, L + 1):
            aux = 0
            for ii in range(len(indice2) - 1, -1, -1):
                aux = aux + indice1[ii]
                indice2[ii] = aux
            indice1 = indice2.copy()

            aux1 = np.outer(PL[k, :], P1)
            aux2 = len(P1)
            P1_new = []
            for iii in range(n):
                for j in range(aux2 - indice2[iii], aux2):
                    P1_new.append(aux1[iii, j])
            P1 = np.array(P1_new)
            P[k, ni[l - 1] - 1:ni[l] - 1] = P1

        Y[k, 0] = y[k]

    # P = np.column_stack([np.ones((N, 1)), P])

    return P, Y, ni
