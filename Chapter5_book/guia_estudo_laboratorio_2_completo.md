# Guia de estudo — Laboratório 2 completo (Capítulo 5)

Este documento acompanha o notebook `roteiro_laboratorio_2_completo.ipynb`, que implementa as **seis questões** do 2º Roteiro de Laboratório. A ideia é construir o entendimento em camadas: primeiro os conceitos, depois o que cada trecho de código faz e, por fim, como interpretar os gráficos e tabelas.

Se este for seu primeiro contato com identificação de sistemas, leia a Seção 1 (teoria) em ordem; se já passou pelo `guia_estudo_laboratorio_2.md`, você pode pular direto para as Seções referentes às Questões 3–6 (que são novas neste roteiro).

---

## 0. Visão geral da bancada

O laboratório aborda **identificação paramétrica** de sistemas dinâmicos em tempo discreto, percorrendo três modos de trabalho:

1. **Análise teórica** (Q1): saímos de funções de transferência contínuas $G(s)$, olhamos polos/zeros, discretizamos via ZOH e extraímos a equação às diferenças.
2. **Simulação controlada + EMQ** (Q2 e Q3): geramos dados sintéticos a partir de um modelo conhecido, ajustamos modelos ARX de várias ordens e validamos a escolha usando métricas objetivas.
3. **Identificação em dados reais / simulados com ruído** (Q4 a Q6): para seis conjuntos de dados (`dados_1.txt` … `dados_6.txt`) descobrimos a estrutura e os parâmetros de modelos **ARX** e **ARMAX**, recorrendo a mínimos quadrados estendidos, mínimos quadrados recursivos e critérios de informação (AIC/BIC).

Todas as ferramentas são implementadas à mão sobre NumPy, tanto para fixar os conceitos quanto para permitir leitura linha a linha.

---

## 1. Conceitos fundamentais (consolidados)

### 1.1 Sistemas lineares invariantes no tempo e função de transferência

Um **SLIT contínuo** relaciona uma entrada $u(t)$ e uma saída $y(t)$ através de equações diferenciais lineares com coeficientes constantes. No domínio de Laplace (condições iniciais nulas):

$$G(s) = \frac{Y(s)}{U(s)} = \frac{N(s)}{D(s)}.$$

- **Polos**: raízes de $D(s)$ — determinam os modos naturais.
- **Zeros**: raízes de $N(s)$ — moldam o transitório (sobressinal, antecipação/atraso).
- **Estabilidade BIBO**: todos os polos com parte real negativa.

### 1.2 Forma padrão de segunda ordem

$$G(s) = \frac{\omega_n^2}{s^2 + 2\zeta\omega_n s + \omega_n^2}.$$

Para $0 < \zeta < 1$: polos complexos $s = -\zeta\omega_n \pm j\omega_n\sqrt{1-\zeta^2}$, resposta subamortecida com oscilação.

### 1.3 Discretização ZOH e equação às diferenças

No controle digital a entrada é constante entre amostras. O **retentor de ordem zero (ZOH)** implementa exatamente isso. A discretização ZOH de $G(s)$ produz $G(z) = N(z)/D(z)$, que corresponde à **equação às diferenças causal**

$$\underbrace{\sum_{i=0}^{n_a} d_i\, y[k-i]}_{D(z)\,Y(z)} = \underbrace{\sum_{j=0}^{n_b-1} n_j\, u[k-1-j]}_{U(z)\,N(z)},$$

que é a forma implementada recursivamente no laboratório.

### 1.4 Identificação: do sinal ao modelo

Diferente da **simulação** (saímos do modelo e geramos $y$), a **identificação** parte de pares $(u[k], y[k])$ medidos e **estima** os parâmetros do modelo que melhor os explica.

**Estruturas paramétricas** comuns:

| Estrutura | Equação | Parâmetros |
|-----------|---------|-----------|
| ARX($n$)   | $y[k] = \sum a_i y[k-i] + \sum b_j u[k-j] + e[k]$ | $2n$ |
| ARMAX($n$) | $y[k] = \sum a_i y[k-i] + \sum b_j u[k-j] + \sum c_m e[k-m] + e[k]$ | $3n$ |
| OE         | $y[k] = \frac{B(z)}{A(z)} u[k] + e[k]$ | depende |
| BJ         | genérica | depende |

### 1.5 Estimador de Mínimos Quadrados (EMQ) — forma fechada

Para o modelo ARX linear nos parâmetros, escrevemos $Y = \Phi\theta + e$, com linhas de $\Phi$ formadas por $[y[k-1], \dots, y[k-n], u[k-1], \dots, u[k-n]]$ e minimizamos $\|Y - \Phi\theta\|^2$:

$$\hat\theta = (\Phi^T\Phi)^{-1}\,\Phi^T Y.$$

Na prática usamos `np.linalg.lstsq` — **numericamente mais estável** que inverter explicitamente.

### 1.6 Ruído: dinâmico vs sensor

| Tipo | Onde entra | Efeito no EMQ |
|------|-----------|---------------|
| **Dinâmico** ($y[k] = \text{eq limpa} + e[k]$, realimentado) | Dentro da recursão | EMQ é **consistente e não tendencioso** (hipótese-padrão do ARX). |
| **Sensor** ($Y[k] = y[k] + E[k]$, somado à saída pronta) | Fora da recursão | Regressor $\Phi$ carrega $E[k-i]$ → **errors-in-variables** → **viés severo**. |

Quando predomina ruído de sensor, usam-se **EMQE (estendido)**, **variáveis instrumentais** ou estruturas **Output-Error**.

### 1.7 Mínimos Quadrados Estendido (para ARMAX)

O regressor ARMAX contém $e[k-m]$, que é inobservável. A versão **iterativa** resolve isso:

1. Faz uma estimativa inicial do ruído (tipicamente $\hat e = 0$ ou resíduo do ARX).
2. Monta $\Phi$ com $\hat e$ no lugar de $e$.
3. Resolve $\hat\theta = (\Phi^T\Phi)^{-1}\,\Phi^T Y$.
4. Recalcula o resíduo: $\hat e = Y - \Phi\hat\theta$.
5. Repete até convergir.

### 1.8 EMQ Recursivo (RLS) e Estendido (RELS)

Permite atualizar $\hat\theta$ amostra a amostra, sem rearmar a matriz $\Phi^T\Phi$:

$$K[k] = \frac{P[k-1]\,\varphi[k]}{\lambda + \varphi[k]^T P[k-1]\,\varphi[k]}, \quad \hat\theta[k] = \hat\theta[k-1] + K[k]\,(y[k] - \varphi[k]^T\hat\theta[k-1]),$$

$$P[k] = \frac{1}{\lambda}\,(P[k-1] - K[k]\,\varphi[k]^T\,P[k-1]).$$

- $\lambda$ é o **fator de esquecimento** — $\lambda = 1$ dá peso igual a todos os dados; $\lambda < 1$ favorece dados recentes (útil para sistemas variantes no tempo).
- No **RELS**, $\varphi[k]$ inclui os resíduos passados $\hat e[k-m]$ — estende RLS para ARMAX.

### 1.9 Métricas de validação

- **SSE** (Sum of Squared Errors): $\sum (y[k] - \hat y[k])^2$ — valor absoluto do erro acumulado.
- **MSE**: $\mathrm{SSE}/N$ — erro médio (comparável entre datasets).
- **Coeficiente de correlação múltipla** $R^2 = 1 - \dfrac{\mathrm{SSE}}{\mathrm{SST}}$, com $\mathrm{SST} = \sum (y - \bar y)^2$. Está em $(-\infty, 1]$; quanto mais próximo de 1, melhor a fração de variância explicada.
- **Razão sinal-ruído**: $\mathrm{SNR} = 10\log_{10}\!\bigl(\mathrm{var}(y)/\mathrm{var}(y - \hat y)\bigr)$, em dB. Valores altos indicam resíduo pequeno comparado à variância do sinal.

### 1.10 Critérios de informação (AIC e BIC)

Quando se varre ordens, o ajuste (SSE) **sempre melhora** aumentando $n$. Os critérios penalizam a complexidade:

- **Akaike (AIC):** $\mathrm{AIC}(p) = N\ln(\mathrm{SSE}/N) + 2p$.
- **Bayes (BIC):** $\mathrm{BIC}(p) = N\ln(\mathrm{SSE}/N) + p\ln N$.

$p$ é o número de parâmetros ($2n$ para ARX, $3n$ para ARMAX). **BIC** é mais conservador em grandes $N$ (penalização $\propto \ln N$). Escolhe-se a estrutura que **minimiza** o critério.

---

## 2. Estrutura do notebook — questão por questão

### Questão 1 — Sistemas contínuos, polos, ZOH

**O que o código faz:**

1. Monta `ga` e `gb` com `ct.tf('s')`.
2. Imprime polos e ganho DC usando `ct.poles` e `ct.dcgain`.
3. Plota a resposta ao degrau unitário (`ct.step_response`).
4. Discretiza com `ct.c2d(sys, Ts, 'zoh')`, imprime $G(z)$ em forma fatorada e escreve a **equação às diferenças** em linguagem humana.
5. Sobrepõe degrau contínuo e degrau discreto ZOH para inspeção visual.

**O que observar nos resultados:**

- $G_a(s)$ tem polos em $s \in \{-1, -1 \pm j\}$, ganho DC = 1.
- $G_b(s)$ tem polos em $s = -0{,}5 \pm 1{,}5j$, ganho DC = 1 (para $2{,}5/(s^2+s+2{,}5)$ como no enunciado).
- Na sobreposição degrau contínuo/discreto, os marcadores caem **sobre** a curva contínua — ZOH com $T_s = 0{,}1$ s preserva a dinâmica.

**Conexão com a teoria:** polos com parte real menos negativa em $G_b$ $\Rightarrow$ envoltória cai mais devagar; $\zeta$ menor em $G_b$ $\Rightarrow$ maior sobressinal.

### Questão 2 — Simulação e EMQ

Três sub-etapas:

**2.1 Simulação (100 amostras, sem ruído)**: `simulate_diffeq` implementa a recursão $y[k] = (\sum num_j u[k-1-j] - \sum den_i y[k-i])/den_0$. É executada para `u_step` e `u_unif` em ambos os sistemas, gerando quatro saídas.

**2.2 EMQ ordem 1–5**: sobre $G_a$ com entrada uniforme (persistente), usamos `lse_fit_arx` para cada ordem e **plotamos o resíduo** em função de $k$, reportando média e RMS.
*Leitura esperada:* RMS cai a $\sim 10^{-16}$ na ordem 3 (ordem real) e estagna depois. O resíduo deixa de ter estrutura. No notebook gerado obtemos exatamente esse comportamento — confirmando que a implementação está correta.

**2.3 Monte Carlo (100 repetições, ordem 3)**: duas condições de ruído gaussiano ($\sigma = 0{,}05$):

- **Dinâmico** (dentro de `simulate_with_dynamic_noise`): $e[k]$ somado a cada passo, realimenta.
- **Sensor**: $E[k]$ somado de uma vez à saída limpa.

Geramos 100 realizações independentes, ajustamos ARX(3) em cada uma e analisamos média/desvio dos seis parâmetros $[a_1, a_2, a_3, b_1, b_2, b_3]$. O gráfico `errorbar` mostra a média $\pm$ desvio dos dois casos, com o valor verdadeiro como referência.

*Leitura esperada:* média dinâmica ≈ verdadeiro (viés pequeno); média sensor **afasta-se do verdadeiro** — ilustração concreta de errors-in-variables.

### Questão 3 — Validação por SSE, R², SNR

Pegamos os $\hat\theta$ da Questão 2 (em `thetas_by_order`) e os aplicamos a **novos dados** (`rng_val = 777`), calculando:

- SSE (absoluto), MSE (normalizado), $R^2$ (fração da variância explicada) e SNR em dB.

A função `predict_arx` implementa a predição um passo à frente: $\hat y[k] = \sum a_i y[k-i] + \sum b_j u[k-j]$ (usando $y$ verdadeiro no regressor, não $\hat y$ acumulado).

*Leitura esperada:* como os dados de Q2 foram gerados sem ruído, ordens $\ge 3$ produzem SSE numericamente nulo, $R^2 \to 1$ e SNR na faixa de centenas de dB. Ordens 1 e 2 aparecem com SNR modesto ($\sim 25$–63 dB). O gráfico sobrepõe $y$ verdadeiro e $\hat y$ para cada ordem — a separação visual some a partir da ordem 3.

### Questão 4 — ARX e ARMAX nos datasets 1 e 2

Para cada arquivo:

1. `load_dataset(path)` lê `dados_X.txt` como `(y, u)` — a **primeira coluna** é a saída, a **segunda** é a entrada.
2. `split_estim_val` divide 70 %/30 %.
3. `run_arx_armax` ajusta ARX e ARMAX para ordens 1–5, retornando MSE de estimação e validação, $R^2$ e SNR.
4. `imprime_tabela` gera a tabela exigida no enunciado (letra **b**).
5. `melhor_modelo` escolhe pelo menor MSE de validação.
6. Plota `MSE x ordem` (escala log) para visualizar parcimônia.

*Leitura esperada nos resultados:* no notebook,
- `dados_1.txt` → **ARMAX ordem 4** (MSE_val ≈ 9,0·10⁻³, $R^2 \approx 0{,}888$).
- `dados_2.txt` → **ARMAX ordem 3** (MSE_val ≈ 9,8·10⁻³, $R^2 \approx 0{,}933$).

O ARMAX tende a superar o ARX quando há autocorrelação no resíduo do ARX — sinal de que a estrutura `+ C(z) e[k]` captura dinâmica do ruído que o ARX atribuiria erradamente a $a_i$.

### Questão 5 — ARX/ARMAX nos datasets 3 e 4 + análise temporal

Mesma varredura da Questão 4 mais um passo adicional: **RLS de ordem 3** aplicado ao dataset completo, com `rls_arx(..., lam=1.0)`. Plotamos $\hat\theta[k]$ componente a componente.

**Como ler a trajetória:**

- Curvas que **estabilizam** após o transitório inicial → sistema **invariante** no tempo; o ajuste em lote é confiável.
- Curvas que continuam **derivando** ou **oscilando** → sistema **variante** (não-estacionário); convém usar $\lambda < 1$ ou segmentar os dados.

*Leitura esperada nos resultados:*
- `dados_3.txt` → melhor modelo ARX ordem 5; $R^2$ de validação moderado (~ 0,75).
- `dados_4.txt` → melhor modelo ARMAX ordem 3; $R^2$ alto (~ 0,98), indicando um sistema bem comportado.

### Questão 6 — AIC/BIC + RELS nos datasets 5 e 6

`varredura_aic_bic` ajusta ARX e ARMAX para cada ordem sobre o **dataset inteiro** e calcula SSE, AIC, BIC. A tabela permite escolher a estrutura que minimiza cada critério.

Depois roda `rels_armax` na ordem selecionada pelo **BIC** (mais parcimonioso) e plota a trajetória de todos os parâmetros $a_i, b_j, c_m$.

*Leitura esperada nos resultados:* no notebook, tanto `dados_5.txt` (1000 pontos) quanto `dados_6.txt` (500 pontos) são selecionados como **ARMAX ordem 4** por AIC e BIC — um raro cenário em que os dois concordam, sinal de que a estrutura é clara.

---

## 3. Mapa mental — conceito ↔ função ↔ gráfico

| Conceito | Função no notebook | Onde aparece | O que observar |
|----------|-------------------|--------------|----------------|
| Resposta a degrau contínuo | `ct.step_response` | Q1, Q5 (implícito) | Subida, oscilação, valor final = ganho DC |
| Polos e estabilidade | `ct.poles`, texto | Q1 | Todos com parte real negativa → estável |
| Discretização ZOH | `ct.c2d(sys, Ts, 'zoh')` | Q1 | Coef. ≈ Taylor do contínuo para $T_s$ pequeno |
| Equação às diferenças | `explicit_diffeq` (inline) | Q1 | Forma `D(z) y = N(z) u` escrita em prosa |
| Simulação recursiva | `simulate_diffeq` | Q2, Q3 | Saída determinística para degrau / uniforme |
| Matriz regressora ARX | `build_regressor_arx` | Q2-Q6 | Shape $(N-n, 2n)$ para ordem $n$ |
| EMQ padrão | `lse_fit_arx` | Q2, Q3, Q4-Q6 | $\theta$ via `lstsq`; resíduos ao lado |
| Ruído dinâmico | `simulate_with_dynamic_noise` | Q2 | `y[k] = acc + e[k]` realimenta |
| EMQE (ARMAX iterativo) | `armax_fit` | Q4-Q6 | Regressor com $\hat e[k-m]$, iterações até convergir |
| Predição 1 passo | `predict_arx`, `predict_armax` | Q3-Q6 | Usa $y$ verdadeiro no regressor |
| Métricas | `metrics` | Q3-Q5 | SSE, MSE, $R^2$, SNR de validação |
| AIC/BIC | `aic_bic` | Q6 | Penaliza complexidade |
| RLS | `rls_arx` | Q5 | $\hat\theta[k]$ evolui amostra a amostra |
| RELS | `rels_armax` | Q6 | RLS estendido para ARMAX |

---

## 4. Como interpretar cada tabela

### Tabela "MSE por ordem" (Q4, Q5)

```
ordem   ARX MSE est   ARX MSE val   ARX R² val   ARMAX MSE est   ARMAX MSE val   ARMAX R² val
```

Leituras típicas:

- **MSE_est diminui monotonicamente com $n$**: mais parâmetros sempre ajustam melhor o conjunto de estimação.
- **MSE_val pode aumentar em ordens altas**: é o **sobreajuste** (overfitting) se manifestando.
- **$R^2$ de validação próximo de 1**: modelo explica bem a variância em dados não vistos.
- **ARMAX melhora o ARX**: o resíduo do ARX tinha estrutura autocorrelacionada.

### Tabela "AIC/BIC por ordem" (Q6)

```
ordem   SSE ARX   AIC ARX   BIC ARX   SSE ARMAX   AIC ARMAX   BIC ARMAX
```

Escolhe-se a linha com **mínimo** de cada critério. Se AIC e BIC discordarem:

- **AIC** tende a pedir ordem mais alta (menos penalização para $N$ grande).
- **BIC** é consistente: se a ordem verdadeira está na família testada, converge a ela quase certamente com $N \to \infty$.

---

## 5. Decisões conceituais que você deve saber defender

Se alguém lhe perguntar, estas são respostas prontas:

1. **"Por que usar entrada uniforme $[-1,1]$ na identificação em vez de degrau?"**  
   Porque o degrau excita apenas o DC; para identificar $2n$ parâmetros é preciso um sinal **persistentemente excitante** — a uniforme distribui energia em toda a banda.

2. **"Por que o EMQ fica viesado com ruído de sensor?"**  
   Porque a matriz $\Phi$ passa a conter $Y[k-i] = y[k-i] + E[k-i]$. O ruído nas variáveis explicativas é o cenário **errors-in-variables**; o EMQ assume $\Phi$ determinística, logo viola a hipótese.

3. **"Como escolher entre ARX e ARMAX?"**  
   Se o resíduo do melhor ARX ainda tem autocorrelação significativa (padrões visíveis, RMS estagnando longe do ruído de medição), ARMAX tende a ajudar. Em termos práticos: compare MSE de validação; se ARMAX reduzir sensivelmente, vale o custo extra dos parâmetros $c_m$.

4. **"Quando usar RLS/RELS em vez do EMQ em lote?"**  
   Quando (a) os dados chegam em streaming, (b) o sistema é variante no tempo, ou (c) se quer ver **como** a identificação converge. Para sistemas estacionários com todos os dados em memória, o EMQ em lote é equivalente (com $\lambda = 1$) mas mais direto.

5. **"Por que AIC e BIC dão ordens diferentes às vezes?"**  
   Porque penalizam a complexidade de forma diferente: AIC favorece ajuste, BIC favorece parcimônia conforme $N$ cresce. Para aplicação prática, **prefira BIC** quando $N$ for grande (≥ 300 amostras).

---

## 6. Como replicar no futuro

1. **Ambiente**: Python 3.11+ com `numpy`, `matplotlib`, `control`, `nbconvert` (`pip install numpy matplotlib control nbconvert`).
2. **Arquivos de dados** em `dataset/dados_1.txt`..`dataset/dados_6.txt`, duas colunas sem cabeçalho (saída, entrada).
3. **Executar** o notebook na ordem — cada célula depende do estado criado anteriormente (especialmente `simulate_diffeq`, `lse_fit_arx`, etc. definidas na primeira célula de código).
4. **Regerar do zero**: se quiser reconstruir o `.ipynb` do zero, rode `python _build_notebook.py` — o script gerador é autoexplicativo.
5. **Sanity-check rápido**: depois de rodar Q1, confira `Polos G_a = [-1, -1±j]` e `Polos G_b = [-0.5 ± 1.5j]`. Se não bater, há erro de digitação na função de transferência.

---

## 7. Mini-exercícios de fixação

1. Recalcule, à mão, $G_a(0)$ e $G_b(0)$ e confirme o valor final do degrau.
2. Aumente `Ts` para 0,5 s na Questão 1 e verifique se os pontos discretos se afastam da curva contínua — por que isso acontece?
3. Na Questão 2, diminua `N` de 100 para 30 e observe se o EMQ ainda consegue recuperar a ordem real. Qual o limite para a ordem 3?
4. Na Questão 4, troque o `frac=0.7` para `0.5` e veja se a escolha de melhor modelo muda — discuta o impacto de ter menos dados para estimação.
5. Na Questão 6, modifique `lam` de `rels_armax` para 0,98 (esquecimento leve) e observe se a trajetória dos parâmetros muda de qualitativa — discuta.
6. Compute o coeficiente de correlação **do resíduo** (ACF) para verificar se o ARX escolhido deixou estrutura não capturada.

---

## 8. Referência rápida — funções do notebook

| Função | Propósito | Onde é usada |
|--------|-----------|--------------|
| `simulate_diffeq(num, den, u)` | Simula saída por equação a diferença | Q2, Q3 |
| `build_regressor_arx(y, u, n)` | Monta $\Phi$, $Y$ para EMQ ARX | Q2, Q3 |
| `lse_fit_arx(y, u, n)` | EMQ ARX em forma fechada | Q2-Q6 |
| `armax_fit(y, u, n)` | EMQE iterativo para ARMAX | Q4-Q6 |
| `predict_arx`, `predict_armax` | Predição um passo à frente | Q3-Q5 |
| `rls_arx(y, u, n, lam, P0)` | RLS com histórico de $\theta[k]$ | Q5 |
| `rels_armax(y, u, n, lam, P0)` | Mínimos Quadrados Recursivo Estendido | Q6 |
| `metrics(y_true, y_pred)` | SSE, MSE, $R^2$, SNR | Q3-Q5 |
| `aic_bic(sse, N, p)` | Critérios AIC e BIC | Q6 |
| `simulate_with_dynamic_noise` | Variante de `simulate_diffeq` com ruído realimentado | Q2 |
| `load_dataset(path)` | Lê `dados_X.txt` como $(y, u)$ | Q4-Q6 |
| `split_estim_val(y, u, frac)` | Particiona estimação/validação | Q4, Q5 |
| `run_arx_armax`, `imprime_tabela`, `melhor_modelo`, `varredura_aic_bic` | Orquestração de experimentos | Q4-Q6 |

---

*Documento gerado como guia de estudo alinhado ao `roteiro_laboratorio_2_completo.ipynb`. Ajuste caminhos de arquivo e parâmetros conforme sua cópia local do projeto.*
