# Desafio 08 - PDS 2026.2 - Resultados numéricos

Gerado automaticamente por `src/desafio08.py`.

## Parâmetros do protótipo analógico

| Grandeza | Valor |
|---|---|
| $\omega_n$ | 12000π = 37699.111843 rad/s (f_n = 6000 Hz) |
| $\zeta$ | 1/120 = 0.0083333333 |
| $Q$ | 60 |
| $\sigma=\zeta\omega_n$ | 100π = 314.159265 rad/s |
| $\omega_d$ | 100π√14399 = 37697.802823 rad/s |
| $T_s$ | 25 µs (F_s = 40 kHz) |

## Item 1 - ha(t)

$$h_a(t)=e^{-100\pi t}\left[\cos(\omega_d t)-\frac{1}{\sqrt{14399}}\sin(\omega_d t)\right]u(t),\quad\omega_d=100\pi\sqrt{14399}$$

Forma exata (SymPy): $(-\sqrt{119}\sin(1100\sqrt{119}\pi t)/1309+\cos(1100\sqrt{119}\pi t))e^{-100\pi t}$

## Item 2 - Coeficientes de H(z)

$r=e^{-\sigma T_s}=0.992176780293$, $\theta=\omega_dT_s=0.942445070585$ rad ($0.2999895832\pi$ rad/amostra)

$$H(z)=\frac{b_0+b_1z^{-1}}{1+a_1z^{-1}+a_2z^{-2}}$$

| Coeficiente | Expressão | Valor |
|---|---|---|
| $b_0$ | $T_s$ | 2.500000000000e-05 |
| $b_1$ | $-T_s r(\cos\theta+\frac{\sigma}{\omega_d}\sin\theta)$ | -1.474755715063e-05 |
| $a_1$ | $-2r\cos\theta$ | -1.166426294332e+00 |
| $a_2$ | $r^2$ | 9.844147633517e-01 |

Polos: $z=0.992176780e^{\pm j0.942445071}$ (|z| < 1 : estável). Zero: $z=0.589902286$.

## Item 3 - Ganho de pico

| Filtro | Ganho de pico | dB |
|---|---|---|
| Analógico (em $\Omega=\omega_n$) | 1.5915494309e-03 | -55.9636 |
| Analógico (máximo exato) | 1.5914941640e-03 | -55.9639 |
| Digital em $\omega_0=0{,}3\pi$ | 1.6040862640e-03 | -55.8954 |

Diferença em relação ao pico analógico (-55.9639 dB): **+0.0685 dB** (+0.7912 %) ; em relação a $|H_a(j\omega_n)|$ (-55.9636 dB): +0.0682 dB -> aliasing desprezivel em F_s = 40 kHz.

Banda de -3 dB do filtro digital: 100.0200 Hz.

### Varredura da taxa de amostragem (aliasing)

| $F_s$ [kHz] | pico [dB] | freq. do pico [Hz] | razão vs. analógico |
|---|---|---|---|
| 40.0 | -55.895 | 5999.90 | 1.008 |
| 20.0 | -55.827 | 5999.65 | 1.016 |
| 15.0 | -55.777 | 5999.10 | 1.022 |
| 12.0 | -49.829 | 6000.00 | 2.026 |
| 8.0 | -55.624 | 2000.20 | 1.040 |

## Desafio de Reflexão - passa-altas

$H_{hp}(z)=\frac{(-1.570546327e-01)-2.136589998e-05z^{-1}}{1-0.854635999z^{-1}}$

| Frequência | Ganho | dB |
|---|---|---|
| DC (z=1) | 1.080570 | +0.6731 |
| Nyquist (z=-1) | 0.084671 | -21.4453 |

O ganho decresce com a frequência: o passa-altas degenerou em passa-baixas.
Causa: $|H_a(j\Omega)|\to1$ quando $\Omega\to\infty$, logo as réplicas da
identidade de Poisson não decaem e a soma diverge; o zero em $s=0$ não e
transportado para $z=1$ (zeros não mapeiam por $z=e^{sT_s}$).

## Observação fina (comprovada em verificacao.py)

A identidade de Poisson exata para este problema exige o termo de
meia-descontinuidade, pois $h_a(0^+)=1$ e $h_a(0^-)=0$:

$$H(e^{j\omega})=\sum_k H_a\!\left(j\frac{\omega-2\pi k}{T_s}\right)+\frac{T_s\,h_a(0^+)}{2}$$

com $T_s h_a(0^+)/2 = 1.250000e-05$ (valor conferido numericamente
com erro relativo < 1e-4 na varredura de frequência).
