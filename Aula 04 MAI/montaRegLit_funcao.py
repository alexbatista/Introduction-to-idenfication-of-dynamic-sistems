import sympy as sp


def montaRegLit_funcao(ny, nu, L):
    k = sp.Symbol('k')

    y = [sp.Function('y')(k - i) for i in range(1, ny + 1)]

    if nu != 0:
        u = [sp.Function('u')(k - i) for i in range(1, nu + 1)]
    else:
        u = []

    tamYU = ny + nu

    indice1 = [1] * tamYU
    indice2 = list(indice1)

    PL1 = list(y + u)
    PL = list(PL1)
    P = list(PL1)

    for l in range(2, L + 1):
        aux = 0
        for i in range(len(indice2) - 1, -1, -1):
            aux = aux + indice1[i]
            indice2[i] = aux
        indice1 = list(indice2)

        aux1 = [[sp.conjugate(PL1[r]) * PL[c] for c in range(len(PL))] for r in range(tamYU)]
        aux2 = len(PL)
        PL = []
        for j in range(tamYU):
            for kk in range(aux2 - indice2[j], aux2):
                PL.append(aux1[j][kk])
        P = P + PL

    # P = [1] + P  # Caso não trabalhe com bias, retirar essa linha

    return P
