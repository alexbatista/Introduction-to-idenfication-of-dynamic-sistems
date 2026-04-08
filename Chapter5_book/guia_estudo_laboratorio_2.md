# Guia de estudo — Laboratório 2 (Capítulo 5)

Este documento acompanha o notebook `roteiro_laboratorio_2.ipynb`. A ideia é construir o entendimento **passo a passo**: primeiro os conceitos, depois o que cada trecho de código faz e como isso se conecta aos gráficos e impressões.

---

## 0. O que este laboratório faz (visão geral)

Você tem **dois sistemas lineares invariantes no tempo (SLIT)** descritos por **funções de transferência** \(G_a(s)\) e \(G_b(s)\) em **malha aberta** (sem realimentação). O roteiro:

1. Simula a **saída** de cada sistema quando a **entrada** vem de um arquivo medido (`dataset/dados_1.txt`).
2. Mostra a **resposta ao degrau unitário** no tempo contínuo.
3. Discute **polos, zeros, estabilidade e forma da resposta** (interpretação qualitativa).
4. **Discretiza** cada sistema com retentor de ordem zero (ZOH), obtém \(G(z)\), imprime a **equação às diferenças** e compara degrau **contínuo vs discreto**.

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

### Célula 6

Vazia no notebook original (pode ser usada para suas próprias anotações ou experimentos).

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
