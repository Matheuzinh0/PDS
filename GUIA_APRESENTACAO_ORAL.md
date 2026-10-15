# GUIA DE APRESENTAÇÃO ORAL — DESAFIO 08 (PDS 2026.2)

**Tema:** Discretização de Filtros IIR por Invariância ao Impulso e Análise de Aliasing Espectral
**Formato:** 10 minutos por dupla = **7 a 8 min de exposição + 2 a 3 min de arguição**
**Arquivos:** `slides_desafio08.pdf` (14 slides) | figuras em `figuras/` | códigos em `src/`

---

## 0. ANTES DE APRESENTAR (checklist de 3 minutos)

```bash
cd desafio_08_invariancia_ao_impulso
python src/verificacao.py     # deve imprimir: TODAS AS VERIFICACOES PASSARAM (23/23)
python src/desafio08.py       # regera figuras/ e resultados/
python src/gerar_slides.py    # regera slides_desafio08.pdf
```

- [ ] Slides em PDF abertos **em tela cheia**, já no slide 1, fonte grande (o professor precisa ler os números).
- [ ] `resultados/saida_desafio08.txt` aberto numa aba (tabela de números para consulta rápida).
- [ ] Celular no silencioso e **um cronômetro de 7:30 rodando** (quem não apresenta controla o tempo e faz sinal em 5:00, 6:30 e 7:15).
- [ ] Decidido quem fala cada bloco (sugestão: A = itens 1–3; B = aliasing + reflexão).

---

## 1. DISTRIBUIÇÃO DO TEMPO (mapa de 7:30)

| # | Slide | Assunto | Tempo | Acumulado | Falante |
|---|-------|---------|-------|-----------|---------|
| 1 | 1 | Título / roteiro em 1 frase | 0:15 | 0:15 | A |
| 2 | 2 | Protótipo analógico e parâmetros | 0:40 | 0:55 | A |
| 3 | 3 | **Item 1** — ha(t) | 0:40 | 1:35 | A |
| 4 | 4 | ha(t) + amostras (fig01) | 0:20 | 1:55 | A |
| 5 | 5 | **Item 2** — coeficientes de H(z) | 0:50 | 2:45 | A |
| 6 | 6 | Polos e zeros (fig05) | 0:25 | 3:10 | A |
| 7 | 7 | **Item 3** — ganho de pico (fig03) | 0:40 | 3:50 | A |
| 8 | 8 | Mecanismo do aliasing (fig04) | 0:35 | 4:25 | B |
| 9 | 9 | Quando o aliasing importa (fig07) | 0:30 | 4:55 | B |
| 10 | 10 | **Reflexão** — passa-altas degenerado (fig06) | 0:45 | 5:40 | B |
| 11 | 11 | Reflexão — a matemática | 0:45 | 6:25 | B |
| 12 | 12 | Rigor extra + verificação | 0:25 | 6:50 | B |
| 13 | 13 | h[n] do filtro digital (fig02, **cortável**) | 0:15 | 7:05 | A |
| 14 | 14 | Conclusões (5 frases numeradas) | 0:35 | 7:40 | A |

**Regra de ouro:** se o cronômetro passou de **6:30** no slide 11, **pule o slide 13** e vá direto para as conclusões.
Se atrasar mais de 1 minuto, corte a explicação da identidade de Poisson (slide 12) e apenas diga a frase de uma linha.

---

## 2. NÚMEROS QUE PRECISAM SAIR DA BOCA (decore estes 10)

| Grandeza | Valor |
|---|---|
| ωn / fn | 12000π = 37 699,1118 rad/s = **6 kHz** |
| ζ  /  Q | 1/120 = 0,008333  /  **60** |
| σ = ζωn | 100π = 314,159 rad/s |
| ωd | 100π√14399 = 37 697,8028 rad/s (fd = 5999,79 Hz) |
| Ts / Fs | **25 µs / 40 kHz** |
| r (raio dos polos) | e^(−σTs) = **0,992176780** |
| θ (ângulo dos polos) | ωdTs = **0,29999π rad/amostra** |
| b0 , b1 | 2,5×10⁻⁵ ; −1,4747557×10⁻⁵ |
| a1 , a2 | −1,166426294 ; 0,9844147634 |
| Ganho de pico: analógico / digital | 1,5915494×10⁻³ (−55,9636 dB) / **1,6040863×10⁻³ (−55,8954 dB)** → **+0,0685 dB** |
| Zero de H(z) | z = 0,58990 (inventado pelo método) |
| Banda de −3 dB | 100,02 Hz (≈ f0/Q) |
| Contraprova de aliasing | 15 kHz → +0,19 dB | 12 kHz → **+6,1 dB** | 8 kHz → pico em 2 kHz |
| Passa-altas degenerado | \|H(1)\| = 1,0806 (+0,67 dB) → \|H(−1)\| = 0,0847 (−21,45 dB) |

---

## 3. ROTEIRO FALADO (o que dizer em cada slide)

### Slide 1 — Título (0:15)
> "Boa tarde. Nosso desafio foi projetar um filtro digital passa-faixa por **invariância ao impulso** a partir do protótipo analógico dado, analisar o **aliasing espectral** com Fs = 40 kHz e responder por que o método **falha** no passa-altas. Os códigos, as figuras e as verificações ficam disponíveis para consulta."

### Slide 2 — Protótipo e parâmetros (0:40)
- Escreva/aponte Ha(s) e diga: "o protótipo é de 2ª ordem com **um par de polos complexos conjugados e um zero na origem**".
- **Frase-chave:** "Comparando Ha(s) com a forma padrão, o termo de 1ª ordem dá 2ζωn = 200π, logo **ζ = 1/120 e Q = 60** — é um filtro **muito seletivo** (banda estreita)."
- Explique a definição do método: "**h[n] = Ts·ha(nTs)**; só depois disso calculamos H(z) = Σ h[n] z⁻ⁿ."

### Slide 3 — Item 1: ha(t) (0:40)
- "Decompomos em frações parciais: s/[(s+σ)²+ωd²] = (s+σ)/[(s+σ)²+ωd²] − σ/[(s+σ)²+ωd²]."
- Resultado: ha(t) = e^(−100πt)[cos(ωd t) − (1/√14399) sen(ωd t)] u(t).
- **Diga os dois pontos de atenção:** (i) ωd = 100π√14399 e não ωn (amortecimento desloca a oscilação); (ii) **ha(0⁺) = 1**, ou seja, existe um **salto em t = 0** — isso vai reaparecer no slide 12.
- "Conferimos simbolicamente no SymPy: o resultado exato sai em função de √119, e coincide com a forma numérica a 10⁻¹⁶."

### Slide 4 — fig01: ha(t) e as amostras (0:20)
- "Este gráfico é a **definição do método**: as hastes são h[n]/Ts = ha(nTs), e elas caem **exatamente sobre a curva** — o método não aproxima nada no tempo."
- "A envoltória é e^(−σt): constante de tempo 1/σ = **3,18 ms**, ou seja, o filtro tem memória longa."

### Slide 5 — Item 2: coeficientes de H(z) (0:50) ⭐ (rubrica: resolução numérica)
- Antes da fórmula, diga: "com **r = e^(−σTs)** e **θ = ωdTs**", e escreva H(z) = Ts[1 − r(cosθ + (σ/ωd)senθ)z⁻¹] / [1 − 2r cosθ z⁻¹ + r²z⁻²].
- Leia os coeficientes **devagar**: b0 = 2,5×10⁻⁵ (= Ts), b1 = −1,4747557×10⁻⁵, a1 = −1,166426294, a2 = 0,9844147634.
- "Polos em **0,99217678·e^(±j0,29999π)** → módulo menor que 1 ⇒ **estável**. O ângulo corresponde a 5999,79 Hz, coerente com os 6 kHz de projeto."
- "Observem um detalhe fino: o método **inventa** um zero em z = 0,58990 — os **zeros não são mapeados** por z = e^(sTs); só os polos são."
- Se perguntarem do ganho em DC: "não é zero, vale 1,2534×10⁻⁵ ≈ Ts/2 — voltaremos nisso no slide 12."

### Slide 6 — fig05: polos e zeros (0:25)
- "Os polos analógicos vão para o plano z por z = e^(sTs): raio e^(−σTs) < 1 ⇒ **estabilidade preservada**."
- "À direita comparo com Fs = 15 kHz: o polo gira para θ = 0,8π (mais perto de Nyquist) — é a mesma ressonância **comprimida** no eixo digital."
- Aponte que **não há zero em z = 1**: "o zero em s = 0 simplesmente desaparece — o que sobrou foi o zero em 0,59."

### Slide 7 — Item 3: ganho de pico (0:40) ⭐
- "ω0 = 2π·6000·Ts = **0,3π rad/amostra** — exatamente a ressonância, porque ωnTs = 0,3π."
- Analítico: "|Ha(jωn)| = 1/(2ζωn) = **1/(200π) = 1,59155×10⁻³** (−55,96 dB)."
- Digital: "|H(e^(j0,3π))| = **1,60409×10⁻³** (−55,90 dB)."
- **Resposta ao item:** "a diferença é **+0,0685 dB (+0,79 %)** em amplitude: **não houve alteração significativa**. O gráfico da esquerda traz um *zoom* justamente porque as curvas praticamente coincidem."
- "Motivo: as réplicas caem em 34 kHz e 46 kHz, onde |Ha| está ~50 dB abaixo do pico; Fs/2 = 20 kHz está muito acima de 6 kHz."

### Slide 8 — fig04: mecanismo do aliasing (0:35) ⭐ (rubrica: reflexão)
- "A TDFT do filtro digital é a **soma das réplicas** de Ha deslocadas de kFs."
- Painel (a), 40 kHz: "réplicas em 34 e 46 kHz, longe da banda-base — nada a temer."
- Painel (b), 15 kHz: "a imagem da componente de **−6 kHz cai em 9 kHz**; como o Nyquist é 7,5 kHz, ela **dobra de volta para 6 kHz**, somando-se ao pico: **+0,19 dB**."
- Painel (c), 8 kHz: "agora é a **própria ressonância** que fica acima de Nyquist: o pico digital aparece em **2 kHz** — o filtro virou outro filtro."

### Slide 9 — fig07: quando o aliasing importa (0:30)
- "Varremos Fs de 12 a 60 kHz: o erro do pico cai de **+6,1 dB** em 12 kHz para +0,19 dB em 15 kHz e +0,07 dB em 40 kHz."
- **Insight:** "o caso crítico é **Fs = 2f0 = 12 kHz**, onde a ressonância cai exatamente em Nyquist e a imagem se soma **em fase**: o pico **dobra** (2,026×, medido)."
- "Abaixo de 12 kHz a ressonância é subamostrada e reaparece em Fs − 6 kHz: em 8 kHz, em 2 kHz (painel b)."

### Slide 10 — Reflexão: passa-altas degenerado (0:45) ⭐⭐ (2,5 pontos da rubrica)
- "Agora o desafio de reflexão: o mesmo método aplicado a Ha(s) = s/(s+2000π)."
- "A transformada inversa é **ha(t) = δ(t) − 2000π e^(−2000πt)u(t)**. Amostrando: h[n] = Ts·δ[n] − Ts·2000π·αⁿu[n], com α = e^(−2000πTs) = 0,8546."
- "O H(z) resultante tem **ganho 1,0806 em DC (+0,67 dB)** e **0,0847 em Nyquist (−21,45 dB)**: o ganho **cai** com a frequência — o passa-altas virou um **passa-baixas**."
- Aponte o zero: "e note: o zero está em z ≈ −1,4×10⁻⁴, praticamente na origem — **não há zero em z = 1**, então a resposta não zera em DC."

### Slide 11 — Reflexão: a matemática (0:45) ⭐⭐
- "A causa é estrutural. A identidade exata do método é a soma de Poisson: H(e^jω) = Σk Ha(j(ω−2πk)/Ts) + Ts·ha(0⁺)/2."
- "No passa-baixas, |Ha| → 0 em alta frequência e a soma converge. No **passa-altas |Ha(jΩ)| → 1** quando Ω → ∞: **cada réplica vale 1**, a série **não converge** — energia de todo o eixo Ω é dobrada para dentro da banda digital."
- "Consequência prática: **não existe Fs que resolva**. É limitação do **método**, não da taxa de amostragem."
- "Segunda razão, algébrica: o que definia o passa-altas era o **zero em s = 0**. Na invariância ao impulso **zeros não mapeiam** (só polos, por z = e^(sTs)): somamos um termo Ts·δ[n] a um bloco passa-baixas de 1ª ordem e o resultado decresce."
- "Solução correta para passa-altas: **transformação bilinear com pré-deformação** — exatamente o método do Desafio 09."

### Slide 12 — Rigor extra e verificação (0:25)
- "Um detalhe que **só aparece na conta fina**: como ha(0⁺) = 1, a soma clássica de réplicas fica **exatamente Ts/2 = 12,5 µs abaixo** da resposta exata. Medimos o desvio: **1,2500×10⁻⁵** contra 1,2500×10⁻⁵ previsto — e isso explica por que o ganho em DC não é zero, mas ≈ Ts/2."
- "Rodamos **23 verificações independentes**: coeficientes por dois caminhos (fórmula fechada × resíduos dos polos), H(e^jω) por fórmula × TDFT × FFT, polos digitais = e^(paTs), pico dobrando em 12 kHz. Comando: `python src/verificacao.py` — sai **23/23**."

### Slide 13 — fig02: h[n] do filtro digital (0:15, cortável)
- "A resposta ao impulso do filtro digital: decai devagar porque os polos estão a 0,9922 do círculo unitário; a banda de −3 dB é de **100 Hz**, coerente com f0/Q = 6000/60."

### Slide 14 — Conclusões (0:35)
1. "Resolvemos os três itens de forma fechada e verificada: ha(t), os quatro coeficientes de H(z) e o ganho de pico."
2. "O filtro digital é **estável** (r = 0,9922 < 1) e de **2ª ordem**, com pico exatamente em 0,3π."
3. "Em Fs = 40 kHz o pico difere do analógico em **+0,07 dB**: **não houve alteração significativa**."
4. "O aliasing é governado por Fs/2f0: em 15 kHz já custa 0,19 dB, em **12 kHz dobra o pico** e abaixo disso a ressonância é subamostrada."
5. "E o ponto principal do desafio: a invariância ao impulso **degenera no passa-altas** porque as réplicas não decaem — a identidade de Poisson deixa de convergir."
- Fechamento: "Obrigado. Ficamos à disposição para as perguntas."

---

## 4. BANCO DE PERGUNTAS DA ARGUIÇÃO (2–3 min)

**P1. "Por que ζ = 1/120?"**
"O termo de 1ª ordem do denominador é 2ζωn. Como ωn = 12000π, temos 2ζωn = 200π ⇒ ζ = 200π/(2·12000π) = 1/120. Isso dá **Q = 60** — filtro estreito: banda de −3 dB teórica f0/Q = 100 Hz, e medimos **100,02 Hz**."

**P2. "Mostre que o pico analógico é 1/(200π)."**
"|Ha(jΩ)| = Ω/√((ωn²−Ω²)² + (2ζωnΩ)²). Em Ω = ωn o primeiro termo zera e sobra ωn/(2ζωn·ωn) = 1/(2ζωn) = 1/(200π) = 1,59155×10⁻³. O máximo exato é ligeiramente maior, em Ωr = ωn√(1−2ζ²) = 37 696,49 rad/s (5999,58 Hz): 1,5914942×10⁻³ — a diferença é de 3×10⁻⁵ dB, ou seja, desprezível."

**P3. "Por que o ganho do filtro digital em DC não é exatamente zero?"**
"Porque a soma Σh[n] é uma **soma de Riemann** de ∫ha(t)dt = Ha(0) = 0, que não é exata, **e** porque a soma de Poisson exige o termo de meia-descontinuidade Ts·ha(0⁺)/2 = 12,5 µs. Medimos |H(e^j0)| = 1,2534×10⁻⁵ ≈ Ts/2."

**P4. "Vocês diriam que houve alteração significativa do pico? Com que critério?"**
"Não. A diferença é **+0,0685 dB** (0,79 % em amplitude), contra um nível de −55,9 dB. Como referência, 0,07 dB é 40 vezes menor que o próprio ripple típico de 3 dB de uma banda passante. Se o critério for 'erro < 0,5 dB', passamos com folga em qualquer Fs acima de ~14 kHz."

**P5. "E se Fs fosse 20 kHz? E 12 kHz?"**
"Dobrando para 20 kHz o pico digital fica **1,016× o analógico** (−55,827 dB, ou seja, +0,14 dB de erro). Com Fs = 12 kHz a ressonância cai **exatamente em Nyquist** (6 kHz = Fs/2): a réplica de k = 1 cai em −6 kHz, ou seja, **em fase** com a componente original, e o pico **dobra**: medimos 3,225×10⁻³ = 2,026× o analógico (+6,1 dB)."

**P6. "Como vocês provam que o sistema digital é estável?"**
"Os polos são z = r·e^(±jθ) com r = e^(−σTs) = 0,9921768 < 1 — logo **ambos dentro do círculo unitário**. Equivalentemente, a2 = r² < 1 e |a1| < 1 + a2. Como o sistema é causal e estritamente próprio, estabilidade BIBO está garantida."

**P7. "Por que o zero de H(z) está em 0,5899 e não em z = 1 (mapeando s = 0)?"**
"Porque a invariância ao impulso **não mapeia zeros**. Só os polos obedecem z = e^(sTs). O numerador de H(z) vem da combinação r·(cosθ + (σ/ωd)senθ), que produz um zero real em −b1/b0 = 0,58990. É um artefato — e é justamente por isso que o método distorce a resposta em frequência de filtros que dependem de zeros (passa-altas, notch)."

**P8. "Qual a diferença entre invariância ao impulso e transformação bilinear?"**
"A invariância ao impulso **amostra a resposta ao impulso** e preserva a resposta no tempo; a TDFT fica como soma de réplicas ⇒ aliasing. A bilinear é um **mapeamento algébrico de frequências** (s = 2/Ts·(1−z⁻¹)/(1+z⁻¹)): não tem aliasing, mas **deforma** o eixo de frequência — por isso exige **pré-deformação** (pre-warping). Como o aliasing depende de |Ha| → 0, a invariância ao impulso só é segura para **passa-baixas/passa-faixa** limitados em banda."

**P9. "O que aconteceria se vocês usassem h[n] = ha(nTs) sem multiplicar por Ts?"**
"A resposta em frequência ficaria multiplicada por **1/Ts = 40 000**: o formato (posição do pico, banda, polos) continua idêntico, mas a escala de ganho muda. É uma convenção: a definição com Ts faz |H(e^jω)| ≈ |Ha(jΩ)| no intervalo de Nyquist, que é o que permite comparar os dois gráficos em dB diretamente."

**P10. "Como vocês garantiram que os números estão certos?"**
"Com `src/verificacao.py`, que checa por caminhos independentes: (i) ha(t) simbólica no SymPy contra a forma fechada; (ii) os quatro coeficientes pela fórmula e pelos **resíduos complexos** dos polos; (iii) H(e^jω) pela fórmula racional contra a **soma direta (TDFT)** e contra a **FFT** de 2 000 000 amostras; (iv) polos digitais = e^(paTs); (v) teste do dobramento em 12 kHz. São 23 verificações, todas passando."

**P11. "O que muda se eu pedir a banda de −3 dB do filtro digital?"**
"Medimos 100,02 Hz, contra 100 Hz previstos por f0/Q — coerente. No plano z, a banda está ligada ao raio dos polos: quanto mais perto de 1, mais estreita a banda."

**P12. "Por que a soma das réplicas não converge no passa-altas? Dê o argumento rigoroso."**
"Porque o termo geral da série é Ha(j(ω−2πk)/Ts) e |Ha(jΩ)| → 1 quando Ω → ∞ (o protótipo tem grau do numerador igual ao do denominador). A série **diverge** (não é absolutamente convergente), o que significa que energia de **todas** as frequências é dobrada para dentro de |ω| ≤ π. Já no passa-faixa do desafio, |Ha| ~ 1/Ω, a série converge condicionalmente — e aí a identidade vale, com o termo de meia-descontinuidade Ts·ha(0⁺)/2."

**P13. "Se o professor mostrar o espectro do sinal no domínio do tempo, o que ele deve ver?"**
"No caso de 15 kHz, um pico de 6 kHz **inflado em 0,19 dB** e a resposta começando a 'encostar' na borda de Nyquist (0,8π). Em 8 kHz, um pico em **2 kHz** que não existe no protótipo analógico."

---

## 5. PLANO B / FRASES DE SEGURANÇA

- **Se travar num número:** "esse valor está no nosso relatório, `resultados/relatorio_resultados.md`, e na tabela do slide" — e leia da tabela. Nunca invente.
- **Se o tempo estourar (sinal aos 6:30):** pule o slide 13, diga só: "a resposta ao impulso do filtro digital segue no material" e vá para as conclusões.
- **Se pedirem "mostre o código":** abra `src/desafio08.py` (seções numeradas 1..6 iguais aos slides) e `src/pds08_core.py` (funções `item1`, `discretiza_invariancia_impulso`, `soma_de_replicas`).
- **Se perguntarem de algo que não foi feito:** "não calculamos esse ponto; a metodologia está pronta para incluir — o passo seria ______" (diga o passo honestamente).
- **Se a dupla se atrapalhar:** sempre volte para as **três frases âncora**:
  1. "Somos um passa-faixa de 6 kHz com Q = 60, amostrado a 40 kHz."
  2. "O pico digital reproduz o analógico com erro de 0,07 dB."
  3. "No passa-altas o método degenera porque as réplicas não decaem."

---

## 6. RESUMO DE UMA LINHA POR ITEM (para conferência final)

| Item do enunciado | Resposta curta |
|---|---|
| 1. ha(t) | e^(−100πt)[cos(ωd t) − (1/√14399)sen(ωd t)]u(t), ωd = 100π√14399 |
| 2. H(z) | b = [2,5×10⁻⁵; −1,4747557×10⁻⁵], a = [1; −1,166426294; 0,9844147634] |
| 3. Pico em ω0 = 0,3π | 1,60409×10⁻³ (−55,895 dB) vs. 1,5915494×10⁻³ (−55,964 dB) → **+0,0685 dB: sem alteração significativa** |
| Reflexão | Passa-altas degenera: |H(1)| = 1,0806 > |H(−1)| = 0,0847; réplicas de ganho ≈ 1 ⇒ soma de Poisson diverge |
| Computacional | 7 figuras (fig01–fig07) + relatório + 23 verificações automáticas |
