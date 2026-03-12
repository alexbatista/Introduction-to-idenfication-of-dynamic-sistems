# Relatório de Identificação de Sistemas Dinâmicos

Este relatório apresenta os resultados da identificação de sistemas dinâmicos por meio da análise de resposta ao degrau, aplicando seis métodos clássicos a cada um dos conjuntos de dados disponíveis.

Todos os métodos retornam um modelo **FOPDT** (*First-Order Plus Dead-Time*):

$$G(s) = \frac{K \cdot e^{-Ls}}{\tau s + 1}$$

onde:

- **K** – Ganho estático
- **τ (tau)** – Constante de tempo
- **L** – Tempo morto (atraso de transporte)

---

## Descrição dos Métodos

### a) Ziegler-Nichols

O método de **Ziegler-Nichols** (malha aberta) traça uma reta tangente no ponto de inflexão da curva de resposta ao degrau. O tempo morto *L* é determinado pela interseção da tangente com o valor inicial da saída, e a constante de tempo *τ* corresponde ao intervalo entre essa interseção e o ponto onde a tangente atinge o valor de regime permanente.

### b) Hägglund

O método de **Hägglund** utiliza dois pontos da resposta ao degrau: 28 % e 63 % de Δy. A constante de tempo é calculada como *τ = 1,5 × (t₆₃ − t₂₈)* e o tempo morto como *L = t₆₃ − τ*.

### c) Smith (1ª ordem)

O método de **Smith de 1ª ordem** é semelhante ao de Hägglund, porém emprega os limiares de 28,3 % e 63,2 % de Δy (correspondentes à constante de tempo de um sistema de 1ª ordem pura). As fórmulas são: *τ = (t₆₃,₂ − t₂₈,₃) / 0,572* e *L = t₆₃,₂ − τ*.

### c) Smith (2ª ordem)

O método de **Smith de 2ª ordem** calcula a razão *r = t₂₈,₃ / t₆₃,₂* e, por meio de uma tabela de interpolação publicada por Smith (1985), obtém o coeficiente *α* tal que *τ = α × (t₆₃,₂ − t₂₈,₃)* e *L = t₆₃,₂ − τ*.

### d) Sundaresan / Krishnaswamy

O método de **Sundaresan e Krishnaswamy** (1977) utiliza os pontos de 35,3 % e 85,3 % de Δy com as seguintes relações: *τ = 0,6669 × (t₈₅,₃ − t₃₅,₃)* e *L = 1,3 × t₃₅,₃ − 0,29 × t₈₅,₃*.

### e) Mollenkamp

O método de **Mollenkamp** emprega três pontos da resposta ao degrau: 20 %, 60 % e 90 % de Δy. Define-se *τ₁ = t₉₀ − t₂₀*, *τ₂ = t₆₀ − t₂₀* e *r = τ₂ / τ₁*, e então aplica-se uma aproximação polinomial para obter *τ* e *L* sem iteração.


---

## Critérios de Avaliação

Para avaliar o quão bem cada modelo identificado reproduz a resposta medida, foram calculados os seguintes índices de desempenho:

| Critério | Fórmula | Descrição |
|----------|---------|----------|
| **MSE** – Erro Médio Quadrático | $\text{MSE} = \frac{1}{N}\sum_{k=1}^{N} e_k^2$ | Média dos quadrados do erro. Penaliza erros grandes. |
| **IAE** – Integral do Módulo do Erro | $\text{IAE} = \int_0^T\|e(t)\|dt$ | Integral do valor absoluto do erro. |
| **ISE** – Integral dos Erros ao Quadrado | $\text{ISE} = \int_0^T e^2(t)\,dt$ | Integral do quadrado do erro; penaliza erros grandes. |
| **ITAE** – Integral do Módulo do Erro × Tempo | $\text{ITAE} = \int_0^T t \cdot \|e(t)\|\,dt$ | Penaliza erros persistentes no tempo. |

Quanto menor o valor de cada índice, melhor o ajuste do modelo.

---

## Resultados por Conjunto de Dados

### conjunto1.txt

#### Parâmetros Identificados

| Método | K | τ (tau) | L (tempo morto) |
|--------|---|---------|----------------|
| a) Ziegler-Nichols | 1.0030 | 0.1210 | 0.0000 |
| b) Hägglund | 1.0030 | 0.0962 | 0.0014 |
| c) Smith (1ª ordem) | 1.0030 | 0.1121 | 0.0000 |
| c) Smith (2ª ordem) | 1.0030 | 0.0906 | 0.0073 |
| d) Sundaresan / Krishnaswamy | 1.0030 | 0.1040 | 0.0000 |
| e) Mollenkamp | 1.0030 | -0.0003 | 0.0249 |

#### Métricas de Erro

| Método | MSE | IAE | ISE | ITAE |
|--------|-----|-----|-----|------|
| a) Ziegler-Nichols | 0.001340 | 0.024534 | 0.001210 | 0.006591 |
| b) Hägglund | 0.000115 | 0.008077 | 0.000104 | 0.003391 |
| c) Smith (1ª ordem) | 0.000561 | 0.016195 | 0.000506 | 0.004863 |
| c) Smith (2ª ordem) | 0.000247 | 0.010591 | 0.000223 | 0.003619 |
| d) Sundaresan / Krishnaswamy | 0.000177 | 0.009435 | 0.000160 | 0.003616 |
| e) Mollenkamp | N/A | N/A | N/A | N/A |

### conjunto2.txt

#### Parâmetros Identificados

| Método | K | τ (tau) | L (tempo morto) |
|--------|---|---------|----------------|
| a) Ziegler-Nichols | -2.1771 | 0.9024 | 2.8828 |
| b) Hägglund | -2.1771 | 0.4383 | 2.9842 |
| c) Smith (1ª ordem) | -2.1771 | 0.5099 | 2.9150 |
| c) Smith (2ª ordem) | -2.1771 | 0.1153 | 3.3095 |
| d) Sundaresan / Krishnaswamy | -2.1771 | 0.4760 | 3.0142 |
| e) Mollenkamp | -2.1771 | 0.0884 | 2.3951 |

#### Métricas de Erro

| Método | MSE | IAE | ISE | ITAE |
|--------|-----|-----|-----|------|
| a) Ziegler-Nichols | 0.028694 | 1.581661 | 0.344855 | 7.663118 |
| b) Hägglund | 0.015497 | 1.183143 | 0.186248 | 6.052895 |
| c) Smith (1ª ordem) | 0.015458 | 1.183904 | 0.185774 | 6.031664 |
| c) Smith (2ª ordem) | 0.038764 | 1.637671 | 0.465890 | 7.666934 |
| d) Sundaresan / Krishnaswamy | 0.015722 | 1.193206 | 0.188945 | 6.063945 |
| e) Mollenkamp | 0.262102 | 3.001313 | 3.150202 | 11.616685 |

### conjunto3.txt

#### Parâmetros Identificados

| Método | K | τ (tau) | L (tempo morto) |
|--------|---|---------|----------------|
| a) Ziegler-Nichols | 0.9973 | 1.1318 | 0.2123 |
| b) Hägglund | 0.9973 | 0.6444 | 0.3028 |
| c) Smith (1ª ordem) | 0.9973 | 0.7500 | 0.2002 |
| c) Smith (2ª ordem) | 0.9973 | 0.4231 | 0.5271 |
| d) Sundaresan / Krishnaswamy | 0.9973 | 0.4298 | 0.4102 |
| e) Mollenkamp | 0.9973 | 0.0380 | 0.3994 |

#### Métricas de Erro

| Método | MSE | IAE | ISE | ITAE |
|--------|-----|-----|-----|------|
| a) Ziegler-Nichols | 0.021292 | 0.668828 | 0.147081 | 1.419263 |
| b) Hägglund | 0.005392 | 0.332253 | 0.037249 | 0.752965 |
| c) Smith (1ª ordem) | 0.006810 | 0.377067 | 0.047043 | 0.829031 |
| c) Smith (2ª ordem) | 0.005120 | 0.322898 | 0.035367 | 0.647918 |
| d) Sundaresan / Krishnaswamy | 0.003501 | 0.287928 | 0.024181 | 0.618894 |
| e) Mollenkamp | 0.028276 | 0.572429 | 0.195326 | 0.798992 |

### conjunto4.txt

#### Parâmetros Identificados

| Método | K | τ (tau) | L (tempo morto) |
|--------|---|---------|----------------|
| a) Ziegler-Nichols | 1.9710 | 0.8875 | 2.1747 |
| b) Hägglund | 1.9710 | 0.7059 | 2.1561 |
| c) Smith (1ª ordem) | 1.9710 | 0.8229 | 2.0400 |
| c) Smith (2ª ordem) | 1.9710 | 0.2536 | 2.6093 |
| d) Sundaresan / Krishnaswamy | 1.9710 | 0.5362 | 2.2968 |
| e) Mollenkamp | 1.9710 | 0.0300 | 2.2979 |

#### Métricas de Erro

| Método | MSE | IAE | ISE | ITAE |
|--------|-----|-----|-----|------|
| a) Ziegler-Nichols | 0.018362 | 0.863920 | 0.146921 | 3.406013 |
| b) Hägglund | 0.010079 | 0.629579 | 0.080547 | 2.534373 |
| c) Smith (1ª ordem) | 0.013373 | 0.729594 | 0.106940 | 2.847478 |
| c) Smith (2ª ordem) | 0.035503 | 0.963930 | 0.284265 | 3.396815 |
| d) Sundaresan / Krishnaswamy | 0.010978 | 0.652460 | 0.087742 | 2.533189 |
| e) Mollenkamp | 0.112464 | 1.456684 | 0.900949 | 4.762649 |

### conjunto5.txt

#### Parâmetros Identificados

| Método | K | τ (tau) | L (tempo morto) |
|--------|---|---------|----------------|
| a) Ziegler-Nichols | 0.6672 | 0.9366 | 1.7128 |
| b) Hägglund | 0.6672 | 0.5013 | 1.8199 |
| c) Smith (1ª ordem) | 0.6672 | 0.5853 | 1.7385 |
| c) Smith (2ª ordem) | 0.6672 | 0.1689 | 2.1549 |
| d) Sundaresan / Krishnaswamy | 0.6672 | 0.3641 | 1.8990 |
| e) Mollenkamp | 0.6672 | 0.0292 | 1.8792 |

#### Métricas de Erro

| Método | MSE | IAE | ISE | ITAE |
|--------|-----|-----|-----|------|
| a) Ziegler-Nichols | 0.007849 | 0.491660 | 0.066866 | 1.586744 |
| b) Hägglund | 0.002784 | 0.297673 | 0.023718 | 0.961210 |
| c) Smith (1ª ordem) | 0.003166 | 0.319279 | 0.026973 | 1.030412 |
| c) Smith (2ª ordem) | 0.004111 | 0.330856 | 0.035022 | 0.976004 |
| d) Sundaresan / Krishnaswamy | 0.002381 | 0.280239 | 0.020285 | 0.885396 |
| e) Mollenkamp | 0.009774 | 0.443286 | 0.083270 | 1.226530 |

### conjunto6.txt

#### Parâmetros Identificados

| Método | K | τ (tau) | L (tempo morto) |
|--------|---|---------|----------------|
| a) Ziegler-Nichols | 0.0020 | 0.6155 | 4.1125 |
| b) Hägglund | 0.0020 | 0.2903 | 4.2094 |
| c) Smith (1ª ordem) | 0.0020 | 0.3375 | 4.1623 |
| c) Smith (2ª ordem) | 0.0020 | 0.0627 | 4.4372 |
| d) Sundaresan / Krishnaswamy | 0.0020 | 0.1412 | 4.3248 |
| e) Mollenkamp | 0.0020 | 0.0329 | 4.2651 |

#### Métricas de Erro

| Método | MSE | IAE | ISE | ITAE |
|--------|-----|-----|-----|------|
| a) Ziegler-Nichols | 0.000000 | 0.001896 | 0.000001 | 0.013354 |
| b) Hägglund | 0.000000 | 0.001537 | 0.000000 | 0.011549 |
| c) Smith (1ª ordem) | 0.000000 | 0.001582 | 0.000000 | 0.011755 |
| c) Smith (2ª ordem) | 0.000000 | 0.001534 | 0.000000 | 0.011479 |
| d) Sundaresan / Krishnaswamy | 0.000000 | 0.001474 | 0.000000 | 0.011220 |
| e) Mollenkamp | 0.000000 | 0.001663 | 0.000000 | 0.012052 |


---

## Resumo Comparativo — Melhor Método por Conjunto

A tabela abaixo indica, para cada conjunto de dados, qual método obteve o **menor valor** em cada critério.

| Conjunto | Melhor MSE | Melhor IAE | Melhor ISE | Melhor ITAE |
|----------|-----------|-----------|-----------|------------|
| conjunto1.txt | b) Hägglund | b) Hägglund | b) Hägglund | b) Hägglund |
| conjunto2.txt | c) Smith (1ª ordem) | b) Hägglund | c) Smith (1ª ordem) | c) Smith (1ª ordem) |
| conjunto3.txt | d) Sundaresan / Krishnaswamy | d) Sundaresan / Krishnaswamy | d) Sundaresan / Krishnaswamy | d) Sundaresan / Krishnaswamy |
| conjunto4.txt | b) Hägglund | b) Hägglund | b) Hägglund | d) Sundaresan / Krishnaswamy |
| conjunto5.txt | d) Sundaresan / Krishnaswamy | d) Sundaresan / Krishnaswamy | d) Sundaresan / Krishnaswamy | d) Sundaresan / Krishnaswamy |
| conjunto6.txt | d) Sundaresan / Krishnaswamy | d) Sundaresan / Krishnaswamy | d) Sundaresan / Krishnaswamy | d) Sundaresan / Krishnaswamy |


---

## Conclusão

Os resultados demonstram que diferentes métodos de identificação podem produzir parâmetros FOPDT distintos para o mesmo conjunto de dados, e que não existe uma única técnica universalmente superior. A avaliação pelos quatro critérios (MSE, IAE, ISE e ITAE) fornece uma visão multidimensional da qualidade de ajuste.

**Método com melhor desempenho geral:** `d) Sundaresan / Krishnaswamy`, que obteve o menor índice em **13** de 24 comparações realizadas (6 conjuntos × 4 critérios).

**Método com pior desempenho geral:** `a) Ziegler-Nichols`, que não obteve o primeiro lugar em nenhum critério para nenhum conjunto de dados avaliado.

**Desempenho por critério de avaliação:**

- **Erro Médio Quadrático (MSE):** melhor método em 3 conjunto(s): `d) Sundaresan / Krishnaswamy`.
- **Integral do Módulo do Erro (IAE):** melhor método em 3 conjunto(s): `b) Hägglund`.
- **Integral dos Erros ao Quadrado (ISE):** melhor método em 3 conjunto(s): `d) Sundaresan / Krishnaswamy`.
- **Integral do Módulo do Erro × Tempo (ITAE):** melhor método em 4 conjunto(s): `d) Sundaresan / Krishnaswamy`.

De modo geral, métodos de dois pontos bem calibrados (como Hägglund e Sundaresan/Krishnaswamy) mostraram-se robustos para sistemas FOPDT com diferentes relações L/τ. O método de Mollenkamp, embora baseado em três pontos, apresentou instabilidade numérica em alguns conjuntos (τ negativo), resultando em modelos inválidos. O método de Smith de 2ª ordem mostrou-se sensível a erros na estimativa dos pontos percentuais em dados ruidosos.
