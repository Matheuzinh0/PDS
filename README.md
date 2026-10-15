# Desafio 08 — PDS 2026.2 (UFPE / DES)

**Discretização de Filtros IIR por Invariância ao Impulso e Análise de Aliasing Espectral**

Protótipo analógico do enunciado:

$$H_a(s)=\frac{s}{s^2+200\pi s+(12000\pi)^2}=\frac{s}{s^2+2\zeta\omega_n s+\omega_n^2},
\qquad T_s=25\ \mu\text{s}\ (F_s=40\ \text{kHz})$$

---

## 1. Como executar

```bash
# ambiente (Python 3.13; a máquina é um Core 2 Duo, sem SSE4.2 → numpy 2.2.6)
python3 -m venv .venv_pds
.venv_pds/bin/pip install "numpy==2.2.6" "scipy==1.15.3" "matplotlib==3.10.0" "sympy==1.13.3"

# 1) checagens independentes (deve terminar com "TODAS AS VERIFICACOES PASSARAM")
.venv_pds/bin/python src/verificacao.py

# 2) resolução completa: imprime os resultados e gera figuras/ e resultados/
.venv_pds/bin/python src/desafio08.py

# 3) slides em PDF (slides_desafio08.pdf)
.venv_pds/bin/python src/gerar_slides.py
```

## 2. Arquivos

| Arquivo | Conteúdo |
|---|---|
| `src/pds08_core.py` | Núcleo de cálculo (funções puras, sem gráficos): ha(t), coeficientes de H(z), resposta em frequência, soma de Poisson, filtro passa-altas degenerado |
| `src/desafio08.py` | Script principal: responde os itens 1–3, o desafio de reflexão e gera as 7 figuras + relatório |
| `src/verificacao.py` | 23 verificações independentes (SymPy, resíduos, TDFT/FFT, dobramento em 12 kHz) |
| `src/gerar_slides.py` | Monta `slides_desafio08.pdf` (14 slides 16:9) com as figuras |
| `figuras/fig01..fig07` | Figuras em PNG (slides) e PDF (vetorial) |
| `resultados/relatorio_resultados.md` | Todos os números da resolução em tabelas |
| `resultados/saida_desafio08.txt` | Saída completa do script principal (log) |
| `GUIA_APRESENTACAO_ORAL.md` | **Roteiro falado com cronômetro de 7:30 + banco de perguntas da arguição** |
| `material/material_estudo_desafio08.tex` | **Material de estudo (LaTeX, 34 páginas): 16 capítulos de teoria + 48 exercícios resolvidos** |
| `material/material_estudo_desafio08.pdf` | PDF compilado do material de estudo |
| `material/gabarito_numerico.py` | Calcula **todos** os números que aparecem no material (fonte da verdade) |
| `material/gabarito_saida.txt` | Saída do gabarito numérico (para conferir o texto) |
| `slides_desafio08.pdf` | Entrega (via Python): 14 slides 16:9 gerados por `matplotlib` |
| `slides/apresentacao_desafio08.tex` | **Fonte LaTeX/Beamer da apresentação** (14 frames, notas do orador em `\note{}`) |
| `slides/apresentacao_desafio08.pdf` | PDF compilado do `.tex` (14 páginas) |
| `tools/acentuar_strings.py` | Utilitário interno (acentuação das strings do código) |

### 2.1 Slides em LaTeX (Beamer)

```bash
cd slides
pdflatex -interaction=nonstopmode apresentacao_desafio08.tex
pdflatex -interaction=nonstopmode apresentacao_desafio08.tex   # 2ª passada
```

```bash
# opção A: TeX Live do sistema
sudo apt install texlive-latex-recommended texlive-fonts-recommended \
                 texlive-lang-portuguese
# opção B: TinyTeX (sem sudo)
wget -qO- "https://yihui.org/tinytex/install-bin-unix.sh" | sh
export PATH="$HOME/.TinyTeX/bin/x86_64-linux:$PATH"
tlmgr install beamer pgf translator booktabs babel-portugues
```
### 2.2 Material de estudo (LaTeX)

```bash
cd material
pdflatex -interaction=nonstopmode material_estudo_desafio08.tex   # 3x (sumário/refs)
```

**34 páginas** de revisão, de Sinais e Sistemas até PDS, com **3 exercícios
resolvidos por assunto** (48 no total). Estrutura:

| Cap. | Assunto | Cap. | Assunto |
|---|---|---|---|
| 1 | Sinais e sistemas contínuos | 9 | Transformada $Z$: ROC e inversa |
| 2 | SLIT e convolução | 10 | Função de transferência digital |
| 3 | Laplace e $H(s)$ | 11 | TDFT/DFT/FFT, leakage, janelas |
| 4 | Filtros de 2ª ordem ($\zeta$, $Q$, pico) | 12 | Fase, atraso de grupo, fase mínima |
| 5 | Fourier: espectro e banda | 13 | FIR fase linear (tipos I–IV) e Kaiser |
| 6 | Amostragem, Nyquist e aliasing | 14 | Estruturas IIR e quantização |
| 7 | Reconstrução, ZOH e identidade de Poisson | 15 | **Invariância ao impulso** |
| 8 | Sistemas discretos, BIBO | 16 | **Bilinear e pré-deformação** |

Apêndices: A) formulário de bolso; B) o Desafio 08 resolvido em uma página;
C) como conferir cada número + tabela de ``pegadinhas''.

**Todos os números do material são calculados por**
`material/gabarito_numerico.py` (saída em `gabarito_saida.txt`), que inclui as
varreduras de bits da Seção 14 e a verificação da identidade de Poisson do
Capítulo 7. Se um valor do texto divergir do script, o texto está errado.

## 3. Resultados principais

**Item 1 — resposta ao impulso analógica**

$$h_a(t)=e^{-100\pi t}\left[\cos(\omega_d t)-\frac{1}{\sqrt{14399}}\sin(\omega_d t)\right]u(t),
\qquad \omega_d=100\pi\sqrt{14399}=37\,697{,}8028\ \text{rad/s}$$

* parâmetros: ωn = 12000π rad/s (6 kHz), **ζ = 1/120**, **Q = 60**, σ = 100π rad/s;
* forma exata (SymPy): (−√119·sen(1100√119·π·t)/1309 + cos(1100√119·π·t))·e^(−100πt);
* salto em t = 0: ha(0⁺) = 1 (importante na análise de aliasing).

**Item 2 — coeficientes de H(z)** (com r = e^(−σTs) = 0,992176780, θ = ωdTs = 0,29998958π)

$$H(z)=\frac{b_0+b_1z^{-1}}{1+a_1z^{-1}+a_2z^{-2}}$$

| b0 | b1 | a1 | a2 |
|---|---|---|---|
| 2,500000000×10⁻⁵ | −1,474755715×10⁻⁵ | −1,166426294 | 0,9844147634 |

Polos: z = 0,992176780·e^(±j0,29999π) (|z| < 1 → **estável**). Zero do método: z = 0,58990.

**Item 3 — ganho de pico em ω0 = 2π·6000·Ts = 0,3π rad/amostra**

| | ganho | dB |
|---|---|---|
| analógico, \|Ha(jωn)\| = 1/(2ζωn) = 1/(200π) | 1,5915494×10⁻³ | −55,9636 |
| digital, \|H(e^(j0,3π))\| | 1,6040863×10⁻³ | −55,8954 |

**Diferença: +0,0685 dB (+0,79 %) → não houve alteração significativa** (Fs/2 = 20 kHz ≫ 6 kHz;
réplicas em 34 kHz e 46 kHz, ~50 dB abaixo do pico).

**Aliasing (contraprova por varredura de Fs)**

| Fs | pico (dB) | freq. do pico | razão vs. analógico |
|---|---|---|---|
| 40 kHz | −55,895 | 5999,94 Hz | 1,008× |
| 20 kHz | −55,827 | 5999,65 Hz | 1,016× |
| 15 kHz | −55,777 | 5999,10 Hz (0,8π) | 1,022× |
| **12 kHz** | −49,829 | 6000,00 Hz (**Nyquist**) | **2,026×** |
| 8 kHz | −55,624 | 2000,20 Hz (dobra) | 1,040× |

**Desafio de reflexão — falha no passa-altas** Ha(s) = s/(s+2000π):

* ha(t) = δ(t) − 2000π·e^(−2000πt)u(t) → H(z) = Ts − Ts·2000π/(1 − αz⁻¹), α = 0,854636, zero em z ≈ −1,4×10⁻⁴;
* |H(1)| = 1,0806 (+0,67 dB) e |H(−1)| = 0,0847 (−21,45 dB): **o ganho decresce → virou passa-baixas**;
* causa: |Ha(jΩ)| → 1 quando Ω → ∞ ⇒ **as réplicas da identidade de Poisson não decaem** e a série diverge; nenhum Fs corrige. O zero em s = 0 **não** é transportado para z = 1 (zeros não mapeiam);
* solução correta para passa-altas: transformação bilinear com pré-deformação.

**Detalhe fino comprovado numericamente:** a identidade de Poisson exige o termo de
meia-descontinuidade, pois ha(0⁺) = 1:

$$H(e^{j\omega})=\sum_k H_a\!\left(j\frac{\omega-2\pi k}{T_s}\right)+\frac{T_s\,h_a(0^+)}{2}$$

O desvio medido das réplicas "puras" em relação à resposta exata foi 1,2500×10⁻⁵ contra
Ts/2 = 1,2500×10⁻⁵ previsto (erro relativo final 6,5×10⁻⁵). Esse mesmo termo explica por que
|H(e^j0)| ≈ Ts/2 = 1,2534×10⁻⁵ e não exatamente zero.

## 4. Figuras

| Figura | Conteúdo |
|---|---|
| `fig01_ha_t_e_amostras` | ha(t) analítica com as amostras h[n]/Ts = ha(nTs) sobrepostas |
| `fig02_resposta_impulso_digital` | h[n] do filtro digital (200 amostras) + envoltória |
| `fig03_magnitude_analogico_vs_digital` | magnitude em dB: analógico × digital (Fs = 40 kHz e 15 kHz) + zoom do pico |
| `fig04_aliasing_replicas` | mecanismo: réplicas kFs para Fs = 40/15/8 kHz, banda-base e dobramento |
| `fig05_polos_zeros` | plano s, plano z (Fs = 40 kHz) e comparação com Fs = 15 kHz |
| `fig06_falha_passa_altas` | o passa-altas que virou passa-baixas + réplicas que não decaem |
| `fig07_erro_aliasing_vs_Fs` | erro do pico (dB) e frequência do pico em função de Fs |

## 5. Verificações executadas (`src/verificacao.py`)

23/23 passaram: ha(t) simbólica × fechada; coeficientes por fórmula × resíduos dos polos;
H(e^jω) racional × TDFT × FFT(2·10⁶); identidade de Poisson (com e sem a correção de
meia-descontinuidade); pico analítico × digital; dobramento em 9 kHz → 6 kHz;
ressonância em Nyquist (Fs = 12 kHz) dobrando o pico; polos digitais = e^(paTs);
estabilidade (|z| < 1); Σh[n] ≈ 0.
