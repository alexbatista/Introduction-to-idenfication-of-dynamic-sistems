import sympy as sp

ny = 2  # Número de regressores de saída
nu = 2  # Número de regressores de entrada
L = 2   # Grau de não-linearidade

# Dados para Teste
k = sp.Symbol('k')
a = sp.Function('y')(k - 1)
b = sp.Function('y')(k - 2)
d = sp.Function('u')(k - 1)
e = sp.Function('u')(k - 2)

y = [a, b]
u = [d, e]

YX = y + u
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

    # aux1 = PL1^T * PL (outer product simbólico)
    aux1 = [[sp.conjugate(PL1[r]) * PL[c] for c in range(len(PL))] for r in range(tamYU)]
    aux2 = len(PL)
    PL = []
    for j in range(tamYU):
        for kk in range(aux2 - indice2[j], aux2):
            PL.append(aux1[j][kk])
    P = P + PL

P = [1] + P

for i, term in enumerate(P):
    print(f"P[{i}] = {term}")
