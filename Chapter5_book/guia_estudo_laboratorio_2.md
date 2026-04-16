# Guia de estudo — Laboratório 2 (Capítulo 5)

Este documento acompanha o notebook `roteiro_laboratorio_2.ipynb`. A ideia é construir o entendimento **passo a passo**: primeiro os conceitos, depois o que cada trecho de código faz e como isso se conecta aos gráficos e impressões.

---

## 0. O que este laboratório faz (visão geral)

Você tem **dois sistemas lineares invariantes no tempo (SLIT)** descritos por **funções de transferência** \(G_a(s)\) e \(G_b(s)\) em **malha aberta** (sem realimentação). O roteiro:

1. Simula a **saída** de cada sistema quando a **entrada** vem de um arquivo medido (`dataset/dados_1.txt`).
2. Mostra a **resposta ao degrau unitário** no tempo contínuo.
3. Discute **polos, zeros, estabilidade e forma da resposta** (interpretação qualitativa).
4. **Discretiza** cada sistema com retentor de ordem zero (ZOH), obtém \(G(z)\), imprime a **equação às diferenças** e compara degrau **contínuo vs discreto**.
5. **Simula** a saída pela equação a diferença (100 amostras) para entrada **degrau** e **uniforme** $[-1,1]$, sem ruído.
6. Aplica o **Estimador de Mínimos Quadrados (EMQ)** para ajustar modelos ARX de ordem **1 a 5**, analisando o **resíduo** de cada ajuste.
7. Faz um estudo **Monte Carlo** (100 realizações) adicionando ruído gaussiano de duas formas — **dinâmico** (dentro da equação) e de **sensor** (somado à saída pronta) — para observar como cada tipo afeta **média e desvio** dos parâmetros estimados.

Use este guia na ordem das seções; cada seção prepara a seguinte.

---

## 1. Conceitos fundamentais (definições)

### 1.1 Sistema dinâmico, entrada e saída

Um **sistema dinâmico** (aqui, linear e invariante no tempo) relaciona uma **entrada** \(u(t)\) e uma **saída** \(y(t)\) ao longo do tempo. Em engenharia de controle, costuma-se trabalhar no **domínio de Laplace** para sistemas contínuos: \(U(s)\) e \(Y(s)\) são as transformadas de \(u\) e \(y\).

### 1.2 Função de transferência

A **função de transferência** \(G(s)\) é a relação


$$G(s) = \frac{Y(s)}{U(s)}
$$
supondo condições iniciais nulas. Ela é uma **razão de polinômios em \(s\)**:

- **Polos**: raízes do **denominador** (determinam modos naturais: exponenciais, oscilações amortecidas).
- **Zeros**: raízes do **numerador** (“embasam” ou cancelam efeitos de certas frequências; influenciam sobressinal e forma da resposta).

Variável simbólica \(s\) no Python: no notebook, `s = ct.tf('s')` cria o objeto que permite escrever polinômios em \(s\) e dividir para formar `ct.tf`.

### 1.3 Malha aberta

**Malha aberta** significa que a saída **não** é medida para corrigir a entrada (não há realimentação). Aqui você apenas aplica \(u\) e observa \(y = G(s)\,U(s)\) no modelo — sem controlador fechando a malha.

### 1.4 Estabilidade BIBO (entrada limitada → saída limitada)

Para um SLIT contínuo descrito por \(G(s)\) **racional e próprio**, o sistema é **BIBO estável** se **todos os polos** têm **parte real negativa** (ficam no **semiplano esquerdo** do plano \(s\)). Se algum polo tiver parte real positiva ou zero (casos patológicos à parte), a resposta pode divergir ou não ser limitada para qualquer entrada limitada.

### 1.5 Resposta forçada e resposta ao degrau

- **Resposta forçada**: saída quando a entrada é um sinal específico \(u(t)\) (no notebook, vindo do arquivo).
- **Degrau unitário** \(u(t)=1\) para \(t\ge 0\): é o teste clássico para ver **ganho em regime permanente**, **tempo de subida**, **sobressinal** e **tempo de acomodação**.

No `python-control`, `ct.forced_response(sys, T=..., U=...)` calcula a saída amostrada nos tempos `T` para a entrada `U`. Já `ct.step_response(sys)` calcula a resposta ao degrau (ele escolhe o vetor de tempo internamente, a menos que você passe `T=`).

### 1.6 Sistema de segunda ordem e parâmetros $\zeta$ e $\omega_n$

Forma padrão:


$$\frac{\omega_n^2}{s^2 + 2\zeta\omega_n s + \omega_n^2}
$$

- $\omega_n$: **frequência natural** (rad/s).
- $\zeta$: **fator de amortecimento** ($0<\zeta<1$ → resposta **subamortecida** com oscilação).

Quanto **menor** $\zeta$ (abaixo de ~0,5), em geral **maior** o sobressinal relativo. Polos complexos $s = -\zeta\omega_n \pm j\omega_n\sqrt{1-\zeta^2}$.

### 1.7 Discretização e ZOH

No controle digital, o computador atualiza a cada período \(T_s\). Muito comum supor um **retentor de ordem zero (ZOH)**: o sinal de controle **mantém-se constante** entre instantes $kT_s$ e $(k+1)T_s$.

A conversão **contínuo → discreto** com equivalência ZOH produz uma função de transferência em **\(z\)** (transformada \(Z\)), \(G(z)\), que relaciona sequências \(u[k]\) e \(y[k]\).

### 1.8 Equação às diferenças

$G(z) = \frac{N(z)}{D(z)}$ corresponde a uma **equação às diferenças** linear com coeficientes constantes, que o notebook imprime na forma causal (valores passados de \(y\) e \(u\)). É a versão “recursiva” do filtro/sistema discreto — útil para implementar ou simular no tempo.

### 1.9 Simulação recursiva da equação a diferença

Uma equação a diferença como
$$y[k] = a_1\,y[k-1] + a_2\,y[k-2] + \cdots + b_1\,u[k-1] + b_2\,u[k-2] + \cdots$$
é **recursiva**: cada novo valor de $y$ depende de **valores passados** de $y$ e $u$. Para simular numericamente, basta um **loop** com um acumulador:

1. Para cada instante $k$, some as contribuições das entradas atrasadas ($u[k-1], u[k-2], \ldots$).
2. Some as contribuições (com sinal trocado, se vierem do denominador) das saídas anteriores ($y[k-1], y[k-2], \ldots$).
3. Divida por `den[0]` (geralmente 1) e armazene em `y[k]`.

Nos primeiros passos, índices negativos são tratados como **zero** — isso é equivalente a dizer que o sistema parte do **repouso** (condições iniciais nulas).

### 1.10 Identificação de sistemas

**Identificação** é o processo oposto da simulação: dado pares de amostras $(u[k], y[k])$ medidas, estimar os parâmetros do modelo que melhor descrevem a dinâmica. Quando se assume a estrutura **ARX** (AutoRegressive with eXogenous input):
$$y[k] = \underbrace{a_1 y[k-1] + \cdots + a_n y[k-n]}_{\text{parte autorregressiva (AR)}} + \underbrace{b_1 u[k-1] + \cdots + b_n u[k-n]}_{\text{entrada exógena (X)}} + e[k]$$
a identificação se reduz a encontrar o vetor $\theta = [a_1, \ldots, a_n, b_1, \ldots, b_n]^T$ que minimize $\sum_k e[k]^2$.

### 1.11 Estimador de Mínimos Quadrados (EMQ)

Escrevendo o modelo ARX em forma matricial para todas as $N-n$ amostras disponíveis:
$$Y = \Phi\,\theta + e$$
onde
- $Y = [y[n+1], y[n+2], \ldots, y[N]]^T$ é o **vetor alvo** (saídas medidas);
- $\Phi$ é a **matriz regressora** — cada linha contém $[y[k-1], \ldots, y[k-n], u[k-1], \ldots, u[k-n]]$;
- $\theta$ é o **vetor de parâmetros** a estimar;
- $e$ é o **vetor de resíduos**.

A solução do EMQ minimiza $\|Y - \Phi\theta\|^2$:
$$\hat\theta = (\Phi^T \Phi)^{-1}\,\Phi^T Y$$
Na prática usa-se `np.linalg.lstsq(Phi, Y)` — **numericamente mais estável** que inverter $\Phi^T\Phi$ (que pode ser mal-condicionada). O resultado é o mesmo em teoria.

### 1.12 Resíduo — o que é e como interpretar

O **resíduo** é a diferença entre a saída medida e a predita pelo modelo:
$$r[k] = y[k] - \hat y[k]$$
- Modelo **correto** → resíduo ~ ruído branco (média ≈ 0, sem padrão temporal).
- Modelo com **ordem baixa demais** → resíduo mostra **estrutura** (oscilações, tendências) — há dinâmica não capturada.
- Modelo com **ordem alta demais** (sobreajuste) → resíduo pequeno, mas parâmetros extras ficam mal-determinados (alta variância em presença de ruído).

Indicadores úteis:
- **Média do resíduo**: próxima de zero se não há viés sistemático.
- **RMS do resíduo**: $\sqrt{\frac{1}{N}\sum r[k]^2}$ — grandeza global do erro, como um desvio padrão.

### 1.13 Entrada persistentemente excitante

Para identificar todos os parâmetros, a entrada $u[k]$ precisa **excitar** toda a banda do sistema. Um **degrau** excita apenas o componente DC e é **insuficiente** para estimar ordens altas: o EMQ produz $\Phi^T\Phi$ mal-condicionada (colunas quase linearmente dependentes). Uma entrada **uniforme em $[-1, 1]$** distribui energia em todas as frequências discretas — é dita **persistentemente excitante** e torna a identificação bem-posta.

### 1.14 Ruído dinâmico vs ruído de sensor

Duas formas comuns de ruído afetarem um sistema real:

| Tipo | Onde entra | Equação |
|------|-----------|---------|
| **Dinâmico** (equation error) | Dentro da dinâmica — em cada passo da recursão | $y[k] = (\text{equação limpa}) + e[k]$ |
| **Sensor** (output error) | No medidor, somado à saída já calculada | $Y[k] = y[k] + E[k]$ |

Diferença crítica: no ruído dinâmico, cada $e[k]$ passa pela recursão e afeta $y[k+1], y[k+2], \ldots$ via $y[k-1]$. No ruído de sensor, os $E[k]$ são independentes e **não se propagam** pela dinâmica.

### 1.15 Estudo Monte Carlo e viés do EMQ

Quando o ruído é aleatório, uma **única** identificação é apenas uma amostra: outro experimento com outra realização do ruído dará parâmetros diferentes. Em um estudo **Monte Carlo**, repete-se a identificação $N_{MC}$ vezes (aqui 100), cada uma com uma realização **independente** do ruído, e analisam-se:
- **Média dos parâmetros estimados** → indica se há **viés** (desvio sistemático em relação ao valor verdadeiro).
- **Desvio padrão** → indica a **variância** (incerteza) do estimador.

Resultados teóricos importantes:
- Com **ruído dinâmico branco** e modelo ARX correto, o EMQ é **consistente e não tendencioso** — a média converge para o verdadeiro.
- Com **ruído de sensor**, o regressor $\Phi$ contém $Y[k-i] = y[k-i] + E[k-i]$: as variáveis explicativas **também carregam ruído**. Isso é um problema de **errors-in-variables** e produz **viés severo** no EMQ. Para esse caso, usam-se variantes como **Variáveis Instrumentais** ou **Mínimos Quadrados Estendidos (EMQE)**.

---

## 2. Estrutura do notebook — célula por célula

### Célula 0 (markdown)

Enunciado: considerar os sistemas em **malha aberta** (apresentados no código seguinte).

### Célula 1 — Dados de entrada e resposta forçada

**Imports**

- `control as ct`: biblioteca de sistemas e funções de transferência.
- `matplotlib.pyplot`: gráficos.
- `numpy`: arrays, leitura de arquivo, interpolação.

**Função `load_and_resample_dataset`**

1. `np.loadtxt` lê duas colunas do arquivo (tempo e entrada). `unpack=True` separa em dois vetores.
2. `np.argsort(t_data)` ordena os índices pelo tempo (caso o arquivo não esteja ordenado).
3. `np.linspace` cria uma grade de tempos **uniformemente espaçada** com o mesmo número de pontos, do primeiro ao último tempo original.
4. `np.interp` **interpola** a entrada nos novos tempos — assim o passo no tempo fica regular, o que combina bem com simulação numérica.

**Definição dos sistemas**

- `ga` e `gb` são objetos `tf` (função de transferência) montados com o mesmo `s`.
- `systems` e `labels` permitem o loop sobre os dois sistemas com legendas em LaTeX.

**Simulação**

- Para cada `sys`, `ct.forced_response(sys, T=t_uniform, U=u_uniform)` devolve `time` e `response`: saída \(y(t)\) quando a entrada é `u_uniform` nos instantes `t_uniform`.
- O gráfico mostra **duas curvas**: saída de \(G_a\) e de \(G_b\) para a **mesma entrada medida**. Isso compara **dinâmicas diferentes** frente ao mesmo estímulo.

**Leitura conceitual do gráfico**

- Picos e vales refletem como cada \(G(s)\) **filtra e atrasa** o conteúdo de frequência da entrada.
- Diferenças entre as curvas vêm dos **polos/zeros** e do **ganho DC** de cada função.

### Célula 2 — Resposta ao degrau unitário

- `ct.step_response(sys)` sem `T` explícito: a biblioteca gera uma duração adequada.
- O laço plota \(G_a\) e \(G_b\) no **mesmo eixo** para comparar **subida**, **oscilação** e **valor final**.
- O `print` resume quantidade de pontos, intervalo de tempo e **valor aproximado no último instante** (`y_final`) — proxy do ganho em regime para o degrau (para sistemas estáveis com ganho DC finito).

**Conceito**: degrau = “ligar” uma referência constante; a resposta mostra o **comportamento transiente** típico de cada planta.

### Célula 3 (markdown) — Análise por polos (encaixe com os gráficos)

Esta célula **não executa código**; interpreta \(G_a\) e \(G_b\):

**\(G_a(s)\)**

- Denominador fatorado $(s+1)(s^2+2s+2)$: polos em \(s=-1\) e \(s=-1\pm j\). Todos com parte real negativa → **estável**.
- Numerador \(0.5(s+2)^2\): **dois zeros** em \(s=-2\).
- Ganho DC \(G_a(0)=1\).
- O par complexo dá modo oscilatório; o texto relaciona \(\omega_n\) e \(\zeta\) a partir da forma quadrática. Os **zeros** podem **reduzir sobressinal** e influenciar a forma em relação a um sistema só com aquele par de polos.

**\(G_b(s)\)**

- Segunda ordem, polos \(s=-0.5\pm j1.5\), sem zeros finitos no numerador (só constante).
- \(\zeta\) menor que o do modo oscilatório comparado no texto para \(G_a\) → espera-se **mais sobressinal** e oscilação mais evidente.
- Parte real dos polos **menos negativa** que \(-1\) → envoltória decai **mais devagar** que em \(G_a\) (em geral resposta “mais lenta”).
- \(G_b(0)=2.8/2.5=1.12\) → valor final do degrau **ligeiramente acima de 1**.

**Como isso “fecha” com a célula 2**

- Maior oscilação e estabelecimento mais lento em \(G_b\) batem com \(\zeta\) menor e polos mais próximos do eixo imaginário.
- \(G_a\) com ganho 1 no degrau vs \(G_b\) com ~1,12 explica **níveis finais diferentes**.

### Célula 4 — Discretização ZOH e equação às diferenças

- `Ts = 0.1` s: período de amostragem.
- `ct.c2d(ga, Ts, method='zoh')` calcula \(G_a(z)\) e \(G_b(z)\) assumindo ZOH na entrada contínua equivalente.

**Função `_poly_z_str`**

- Formata o numerador ou denominador como polinômio em \(z\) (potências decrescentes) para leitura humana.

**Função `explicit_Gz_and_diffeq`**

- Extrai coeficientes `num` e `den` do objeto discreto (`sys_d.num`, `sys_d.den`).
- Imprime \(G(z)=N/D\) explícito.
- Monta a **equação às diferenças** no formato causal usado pelo `python-control` para relação entrada-saída (combinações de \(y[k], y[k-1], \ldots\) e \(u[k-1], u[k-2], \ldots\)).

**Conceito**: passar de \(G(s)\) para \(G(z)\) é o passo entre **modelo contínuo** (PLC/microcontrolador trabalha em passos \(k\)) e **implementação/simulação discreta**.

### Célula 5 — Degrau: contínuo vs discreto

- `t_cont`: tempo contínuo denso (`linspace`).
- `t_disc`: instantes \(0, T_s, 2T_s, \ldots\) até `t_end`.
- Para cada par contínuo/discreto:
  - `ct.step_response(sys_c, T=t_cont)` — degrau no modelo contínuo.
  - `ct.step_response(sys_d, T=t_disc)` — degrau no modelo **ZOH** discreto nos mesmos instantes.

**Leitura do gráfico**

- A curva contínua é a referência física idealizada.
- Os marcadores discretos devem **acompanhar** a contínua nos instantes de amostragem se \(T_s\) for pequeno o suficiente; diferenças visíveis podem aparecer se \(T_s\) for grande em relação às constantes de tempo do sistema.

### Célula 6 (markdown) — Enunciado da Questão 2

Texto da segunda questão: gerar dados pela equação a diferença, identificar com EMQ e estudar o efeito de ruído. As três células seguintes implementam cada item.

### Célula 7 — Parte 1 da Questão 2: Simulação pela equação a diferença

**Objetivo:** produzir `y[k]` a partir dos coeficientes discretos, para duas entradas distintas, durante 100 passos, **sem ruído**.

**Imports e constantes**

- `N = 100`: número de amostras.
- `rng_in = np.random.default_rng(seed=42)`: gerador pseudoaleatório com semente fixa (reprodutibilidade).
- `num_a`, `den_a`, `num_b`, `den_b`: coeficientes já obtidos na Célula 4 (`ga_d`, `gb_d`). Duplicados aqui como constantes **para a célula ser independente**, mas o valor é idêntico ao que `ga_d.num/den` devolve.

**Função `simulate_diffeq(num, den, u)`**

Implementa o loop recursivo explicado em 1.9:

1. Inicializa `y = np.zeros(N)` (condições iniciais nulas).
2. Para cada `k`: soma contribuições de entradas atrasadas (`num[j] * u[k-1-j]`) e subtrai contribuições de saídas anteriores (`den[i] * y[k-i]`).
3. Divide por `den[0]` (= 1 neste caso).

Indícios negativos são ignorados pela checagem `if idx >= 0:` — equivale a assumir zero nesse passado.

**Entradas geradas**

- `u_step = np.ones(N)` → **degrau unitário**.
- `u_unif = rng_in.uniform(-1, 1, size=N)` → sinal uniforme com média ≈ 0, para identificação persistente (ver 1.13).

**Saídas calculadas**

Quatro combinações: `y_a_step`, `y_a_unif`, `y_b_step`, `y_b_unif`.

**Gráfico 2×2**

- Linha 1: $G_a$ (degrau à esquerda, uniforme à direita).
- Linha 2: $G_b$ nos mesmos cenários.

**Leitura do gráfico**

- **Degrau**: cada sistema tende ao seu ganho DC ($G_a \to 1$, $G_b \to 1.12$) com oscilação amortecida.
- **Uniforme**: saída "ruidosa" na aparência, mas determinística — é a entrada arbitrária que o sistema filtra.

### Célula 8 — Parte 2 da Questão 2: EMQ de ordem 1 a 5

**Objetivo:** ajustar modelos ARX de várias ordens à saída simulada e observar o comportamento dos resíduos.

**Função `build_regressor(y, u, n)`**

- Monta `Phi` de forma $(N-n) \times 2n$: cada linha contém `[y[k-1], ..., y[k-n], u[k-1], ..., u[k-n]]`.
- Monta `Y` com os valores `y[k]` correspondentes, para $k = n, \ldots, N-1$.
- O loop começa em `k = n` porque antes disso não há $n$ amostras passadas completas.

**Função `lse_fit(y, u, n)`**

- Chama `build_regressor` e resolve $\hat\theta = \arg\min \|Y - \Phi\theta\|^2$ via `np.linalg.lstsq` (ver 1.11).
- Calcula `y_hat = Phi @ theta` e `res = Y - y_hat`.

**Loop sobre ordens 1..5**

Para cada ordem:

1. Ajusta o modelo.
2. Computa `media = res.mean()` e `rms = sqrt(mean(res²))`.
3. Plota o resíduo vs. `k`, com linha horizontal vermelha na média.
4. Imprime uma tabela: `ordem | media | RMS | # parâmetros`.

**Leitura esperada do resultado**

Como $G_a$ é de 3ª ordem:

| Ordem | RMS | Interpretação |
|-------|-----|---------------|
| 1     | alto (~1e-3) | modelo simples demais — resíduo mostra estrutura oscilatória |
| 2     | médio (~1e-5) | melhora, mas ainda incompleto |
| **3** | **~1e-16** | **ordem correta** — resíduo vira apenas erro numérico |
| 4-5   | ~1e-16 | não melhora (sistema já foi capturado); parâmetros extras ≈ 0 |

**Mensagem-chave:** a ordem ideal é aquela a partir da qual o RMS **estagna**. Ir além só adiciona parâmetros mal-determinados.

### Célula 9 — Parte 3 da Questão 2: Monte Carlo com ruído

**Objetivo:** repetir 100 vezes a identificação de 3ª ordem em duas configurações de ruído e comparar **média e desvio** dos parâmetros estimados.

**Configuração**

- `sigma = 0.05`: desvio padrão do ruído gaussiano.
- `N_MC = 100`: número de realizações Monte Carlo.
- `n_model = 3`: ordem do modelo ajustado (igual à ordem real de $G_a$).

**Função `simulate_with_dynamic_noise(num, den, u, sigma, rng)`**

Cópia de `simulate_diffeq`, mas adiciona `e[k] ~ N(0, sigma²)` **dentro do loop**:

```python
y[k] = acc / den[0] + e[k]
```

Como `y[k]` volta a ser usado no próximo passo via `y[k-1]`, o ruído **passa pela dinâmica** (ver 1.14).

**Entrada fixa**

- `u_mc`: um único sinal uniforme — **igual em todas as 100 realizações**. Só o ruído varia, isolando o efeito dele.
- `y_limpo = simulate_diffeq(num_a, den_a, u_mc)`: pré-calculada para o caso sensor.

**Loop Monte Carlo (100 iterações)**

Em cada iteração `i`:

1. `rng_i = np.random.default_rng(seed=1000 + i)`: nova semente → realização independente.
2. **Caso A (dinâmico)**: simula com ruído injetado, ajusta ordem 3 → `thetas_dyn[i]`.
3. **Caso B (sensor)**: gera `E`, soma em `y_limpo`, ajusta ordem 3 → `thetas_sen[i]`.

**Análise estatística**

- `theta_true = np.concatenate([-den_a[1:], num_a])`: valores verdadeiros reorganizados no formato $[a_1, a_2, a_3, b_1, b_2, b_3]$ (note o sinal trocado em `den_a[1:]` — ver 1.11).
- Imprime tabela com média e desvio padrão de cada parâmetro nos dois casos.

**Gráfico `errorbar`**

- Bolinhas (ruído dinâmico) e quadrados (sensor) mostram **média** dos 100 ajustes, com **barras de erro = desvio padrão**.
- `×` preto marca o valor verdadeiro.

**Leitura esperada do resultado**

| Parâmetro | Verdadeiro | Média dyn ≈ verdadeiro? | Média sensor ≈ verdadeiro? |
|-----------|-----------:|:-----------------------:|:--------------------------:|
| $a_1, a_2, a_3$ | — | **sim** (viés pequeno) | **não** (viés severo) |
| $b_1, b_2, b_3$ | — | **sim** | aproximadamente, com maior viés |

**Conclusão** (ver 1.15): ruído dinâmico branco é **exatamente** a hipótese do EMQ, então ele permanece **não tendencioso**. Ruído de sensor contamina o regressor $\Phi$ (pois $Y[k-i]$ passa a ter ruído), caracterizando **errors-in-variables** — o EMQ clássico falha e precisa ser substituído por técnicas mais robustas (Variáveis Instrumentais, EMQE, Output Error).

---

## 3. Mapa mental: do conceito ao resultado na tela

| Conceito | Onde aparece no código | O que você observa |
|----------|-------------------------|---------------------|
| Entrada arbitrária \(u(t)\) | `load_and_resample_dataset`, `forced_response` | Curvas de saída para `dados_1.txt` |
| Degrau unitário | `step_response` (célula 2 e 5) | Subida, oscilação, valor final |
| Polos/zeros/estabilidade | Markdown célula 3 (+ teoria) | Explicação do formato das curvas |
| ZOH e \(G(z)\) | `c2d(..., 'zoh')` | Polinômios impressos em \(z\) |
| Implementação recursiva | `explicit_Gz_and_diffeq` | Equação às diferenças impressa |
| Fidelidade discreta | Célula 5 | Marcadores sobre a curva contínua |
| Simulação recursiva | `simulate_diffeq` (Célula 7) | Saídas `y_a_step`, `y_a_unif`, `y_b_step`, `y_b_unif` |
| Entrada persistente | `rng.uniform(-1, 1, 100)` | Gráficos 2×2 com sinais variados |
| Matriz regressora ARX | `build_regressor` (Célula 8) | `Phi` de forma $(N-n, 2n)$ |
| EMQ | `lse_fit` via `np.linalg.lstsq` | Vetor `theta` de parâmetros ajustados |
| Resíduo e ordem correta | Loop de ordens 1..5 (Célula 8) | RMS despenca na ordem real do sistema |
| Ruído dinâmico | `simulate_with_dynamic_noise` (Célula 9) | Viés pequeno no Monte Carlo |
| Ruído de sensor | `y_limpo + E` (Célula 9) | Viés grande no Monte Carlo |
| Monte Carlo | Loop `for i in range(N_MC)` | `errorbar` com média ± desvio vs. valor verdadeiro |

---

## 4. Como replicar o experimento no futuro

1. **Ambiente Python** com `numpy`, `matplotlib` e `control` (`pip install control` se necessário).
2. Coloque o arquivo de dados em `dataset/dados_1.txt` (duas colunas: tempo e entrada), relativo ao diretório de trabalho do notebook, ou ajuste o caminho em `load_and_resample_dataset(...)`.
3. Execute as células **na ordem** (as variáveis `systems`, `labels`, `ga`, `gb` são reutilizadas após a primeira célula de código).
4. Se mudar `Ts`, execute de novo a célula de `c2d` e a de comparação de degraus.
5. Para conferir à mão: verifique polos com álgebra (ou `ct.poles(ga)` / `ct.poles(gb)` em uma célula extra) e compare com o markdown.

---

## 5. Sugestão de mini-exercícios (fixação)

1. Usando apenas o denominador de \(G_a(s)\), confira as raízes (polos) e compare com o texto da célula 3.
2. Calcule \(G_a(0)\) e \(G_b(0)\) substituindo \(s=0\) nas funções e relacione com o valor final no gráfico do degrau.
3. Aumente `Ts` (por exemplo 0,5 s) e observe se os pontos discretos no último gráfico se afastam da curva contínua — interprete em termos de **perda de informação** entre amostras.
4. Escreva à mão a equação às diferenças impressa para um dos sistemas e simule um passo \(y[k]\) com \(u\) degrau e valores iniciais nulos.

---

## 6. Referência rápida — funções `python-control` usadas

| Função | Papel |
|--------|--------|
| `ct.tf('s')` | Variável simbólica \(s\) para montar razões de polinômios |
| `ct.forced_response(sys, T, U)` | Resposta a entrada `U` nos tempos `T` |
| `ct.step_response(sys, T=...)` | Resposta ao degrau; opcionalmente fixa o vetor de tempo |
| `ct.c2d(sys, Ts, method='zoh')` | Discretização com equivalência ZOH |

---

*Documento gerado como guia de estudo alinhado ao `roteiro_laboratorio_2.ipynb`. Ajuste caminhos de arquivo e parâmetros conforme sua cópia local do projeto.*
