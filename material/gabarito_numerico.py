#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Autor: Matheus Sousa, LCSI-DES, UFPE.
"""
===============================================================================
 MATERIAL DE ESTUDO (Desafio 08) - GABARITO NUMERICO
===============================================================================
 Todo numero que aparece nos exercicios do material de estudo e calculado aqui.
 Rode:

     /home/matheus/.venv_pds/bin/python material/gabarito_numerico.py

 e confira a saida contra o texto do .tex. Se algo divergir, o texto esta errado.
===============================================================================
"""

from __future__ import annotations

import numpy as np
from scipy import signal

np.set_printoptions(precision=10, suppress=False, linewidth=120)
P = print


def sec(n, t):
    P("\n" + "=" * 78)
    P(f"[{n}] {t}")
    P("=" * 78)


# ===========================================================================
# TOPICO 1 - sinais continuos: energia, potencia, periodicidade
# ===========================================================================
sec("1", "Sinais continuos")
# 1.1  x(t) = 3*cos(200*pi*t + pi/4)
A, f = 3.0, 100.0
P(f"1.1 A={A}, f={f} Hz, T={1/f:.6f} s = 10 ms, potencia media P=A^2/2={A**2/2:.4f} W")
# 1.2  energia de e^{-2t}u(t)  -> 1/(2a) = 1/4
a = 2.0
P(f"1.2 energia = 1/(2a) = {1/(2*a):.6f} J")
# 1.3  x(t) = 2cos(4pi t) + 3sen(6pi t)
P("1.3 f1=2 Hz, f2=3 Hz -> razao 3/2 racional -> T0 = 1/f0 com f0 = mdc(2,3) = 1 Hz")
P(f"    periodo fundamental T0 = {1.0:.4f} s; P = 2^2/2 + 3^2/2 = {4/2 + 9/2:.4f} W")

# ===========================================================================
# TOPICO 2 - convolucao continua e SLIT
# ===========================================================================
sec("2", "Convolucao continua / SLIT")
# 2.1  y = u * (e^{-3t}u)
a3 = 3.0
P(f"2.1 y(t) = (1/{a3:.0f})(1-e^{-a3:.0f}t)u(t); y(inf) = {1/a3:.6f}")
# 2.2  ret * ret (largura 1) -> triangular
t = np.linspace(-1.5, 1.5, 7)
P(f"2.2 y(t) = tri(t) = 1-|t| em |t|<1 -> amostras: {np.round(1-np.abs(t),6)}")
# 2.3  h(t) = e^{2t}u(-t): energia/area
P(f"2.3 integral |h| = 1/2 = {0.5:.6f} -> BIBO estavel (anti-causal)")

# ===========================================================================
# TOPICO 3 - Laplace
# ===========================================================================
sec("3", "Transformada de Laplace")
P("3.1 L{e^{-2t}cos(3t)u(t)} = (s+2)/((s+2)^2+9)  [pode-se conferir ponto a ponto]")
s = 1j * 0.0  # checagem numerica em s = j1: valor esperado (j1+2)/((j1+2)^2+9)
z = (1j + 2) / ((1j + 2) ** 2 + 9)
P(f"    conferencia em s=j: {z:.6f}")
# 3.2 fracoes parciais de (s+4)/((s+1)(s+3))
A1 = (-1 + 4) / ((-1 + 3))
B1 = (-3 + 4) / ((-3 + 1))
P(f"3.2 A={A1:.6f} (em s=-1), B={B1:.6f} (em s=-3) -> x(t)=({A1:.1f}e^-t {B1:+.1f}e^-3t)u(t)")
# 3.3 polos de 1/((s+1)(s+2)): -1, -2 (ambos < 0 -> estavel se causal)
P("3.3 polos s=-1 e s=-2: parte real negativa -> H(s) causal e estavel")

# ===========================================================================
# TOPICO 4 - filtros analogicos de 2a ordem (nucleo do desafio)
# ===========================================================================
sec("4", "Filtros analogicos de 2a ordem")
wn, z = 2 * np.pi * 159.1549, 0.1     # wn = 1000 rad/s
wn = 1000.0
sig = z * wn
wd = wn * np.sqrt(1 - z ** 2)
wr = wn * np.sqrt(1 - 2 * z ** 2)
gmax = 1 / (2 * z * np.sqrt(1 - z ** 2))
P(f"4.1 wn=1000 rad/s, zeta=0.1 -> sigma={sig:.4f}, wd={wd:.4f} rad/s")
P(f"    pico em wr={wr:.4f} rad/s ({wr/(2*np.pi):.4f} Hz), |H|max={gmax:.6f}"
  f" = {20*np.log10(gmax):.4f} dB ; 1/(2zeta)={1/(2*z):.3f}")
P(f"    banda -3 dB = 2*zeta*wn = {2*z*wn:.4f} rad/s = {2*z*wn/(2*np.pi):.6f} Hz")

# 4.2 prototipo do desafio (numeros oficiais)
wnC = 12000 * np.pi
zC = 200 * np.pi / (2 * wnC)
sigC = zC * wnC
wdC = np.sqrt(wnC ** 2 - sigC ** 2)
P(f"4.2 DESAFIO: wn={wnC:.6f} rad/s (fn={wnC/(2*np.pi):.4f} Hz), zeta={zC:.10f}, Q={1/(2*zC):.4f}")
P(f"    sigma={sigC:.6f} rad/s, wd={wdC:.6f} rad/s (fd={wdC/(2*np.pi):.4f} Hz)")
P(f"    |Ha(j*wn)| = 1/(2*zeta*wn) = 1/(200*pi) = {1/(2*zC*wnC):.10e}")
P(f"    pico exato em wr=wn*sqrt(1-2z^2)={wnC*np.sqrt(1-2*zC**2):.4f} rad/s,"
  f" |H|max={1/(2*zC*wnC*np.sqrt(1-zC**2)):.10e}")

# ===========================================================================
# TOPICO 6 - amostragem e aliasing
# ===========================================================================
sec("6", "Amostragem e aliasing")


def alias(f, Fs):
    """Frequencia aparente (0..Fs/2) de uma senoide de frequencia f amostrada a Fs."""
    r = f % Fs
    fa = r if r <= Fs / 2 else Fs - r
    return fa


for f0 in (200, 600, 1300, 2000, 3000):
    P(f"6.1 f={f0:5d} Hz, Fs=4000 Hz, Fs/2=2000 Hz -> aparente {alias(f0,4000):7.2f} Hz")
P("    (600 Hz -> 600; 1300 -> 1300-4000 = -2700 -> 1300; 2000 -> 2000 (Nyquist); 3000 -> 1000)")
P("6.2 decimacao por M=4 (Fs'=1000, Nyquist 500):")
for f0 in (200, 600, 1300):
    ap = f0 % 1000
    ap = ap if ap <= 500 else 1000 - ap
    P(f"    {f0} Hz -> aparente {ap} Hz")
P(f"6.3 Fs=15 kHz, f=9 kHz -> aparente {alias(9000,15000)/1000:.3f} kHz"
  f" (=15-9=6 kHz); f=6 kHz -> {alias(6000,15000)/1000:.3f} kHz")
P(f"    Fs=12 kHz, f=6 kHz -> {alias(6000,12000)/1000:.3f} kHz (exatamente em Nyquist)")

# ===========================================================================
# TOPICO 7 - reconstrucao e Poisson
# ===========================================================================
sec("7", "Reconstrucao, ZOH e identidade de Poisson")
P(f"7.1 ZOH: |H(jw)| = Ts*sinc(f/Fs); em f=Fs/2: sinc(0.5)=sin(pi/2)/(pi/2)={2/np.pi:.6f}"
  f" = {20*np.log10(2/np.pi):.4f} dB")
P(f"    em f=0.4*Fs: sinc(0.4)={np.sinc(0.4):.6f} = {20*np.log10(abs(np.sinc(0.4))):.4f} dB")
# 7.2 Poisson exata com termo de meia-descontinuidade: ha(t)=e^{-1000t}u(t), Ts=1e-4
a = 1000.0
Ts = 1e-4
w = np.array([0.0, 0.2 * np.pi, 0.5 * np.pi, 0.9 * np.pi])
exato = Ts / (1 - np.exp(-a * Ts) * np.exp(-1j * w))
soma = sum(1.0 / (a + 1j * (w - 2 * np.pi * k) / Ts) for k in range(-4000, 4001))
soma_c = soma + Ts * 1.0 / 2.0
for wi, e, sc, s0 in zip(w, exato, soma_c, soma):
    P(f"    w={wi/np.pi:.2f}pi: exato={e.real:.8e}{e.imag:+.3e}j  "
      f"poisson+Ts/2={sc.real:.8e}{sc.imag:+.3e}j  (sem correcao: {s0.real:.8e})")
P(f"    diferenca media (sem correcao - exato) = {np.mean(soma - exato).real:.8e}"
  f"  (== -Ts/2 = {-Ts/2:.8e})")

# ===========================================================================
# TOPICO 8 - sistemas discretos: convolucao, BIBO
# ===========================================================================
sec("8", "Sistemas discretos")
n = np.arange(0, 8)
P(f"8.2 h=(1/2)^n u[n]: degrau y[n] = 2-(1/2)^n -> {np.round(2-0.5**n,6)}")
P(f"8.3 h=2^n u[n]: soma diverge (2^n) -> instavel;"
  f" h=(1/2)^|n|: soma total = {sum(0.5**abs(k) for k in range(-60,61)):.6f} (=3)")

# ===========================================================================
# TOPICO 9 - transformada Z
# ===========================================================================
sec("9", "Transformada Z")
# 9.2 X(z)/z = z/((z-0.5)(z-0.25)) -> x[n] = 2(0.5)^n - (0.25)^n
n = np.arange(0, 6)
P(f"9.2 x[n] = 2(0.5)^n - (0.25)^n -> {np.round(2*0.5**n - 0.25**n, 8)}")
P("    conferencia x[0]=1 (numerador e denominador de mesma ordem em z)")
# 9.1, 9.3 (analiticos)

# ===========================================================================
# TOPICO 10 - funcao de transferencia digital
# ===========================================================================
sec("10", "Funcao de transferencia digital")
b, a = np.array([1.0, -0.5]), np.array([1.0, -0.9])
for wf, nome in ((0.0, "DC"), (np.pi, "Nyquist")):
    H = np.polyval(b[::-1], np.exp(-1j * wf)) / np.polyval(a[::-1], np.exp(-1j * wf))
    P(f"10.1 H(z)=(1-0.5z^-1)/(1-0.9z^-1) em {nome}: |H|={abs(H):.6f}"
      f" = {20*np.log10(abs(H)):.4f} dB")
P(f"     polo z={-a[1]:.4f}, zero z={-b[1]:.4f}, ganho maximo em w=0 (filtro passa-baixas)")
# 10.2 polinomio do Desafio 07
A07 = np.array([1.0, -1.2, 0.61, -0.105])
raizes = np.roots(A07)
P(f"10.2 raizes de A(z)=1-1.2z^-1+0.61z^-2-0.105z^-3: {np.round(raizes,6)}")
P(f"     |z| = {np.round(np.abs(raizes),6)} -> todas < 1 : ESTAVEL (k3=+0.105)")
# 10.3 prototipo do desafio
f = None
bC = np.array([2.5e-5, -1.474755715063e-5])
aC = np.array([1.0, -1.166426294332, 0.9844147633517])
for wf, nome in ((0.0, "DC"), (0.3 * np.pi, "w0=0.3pi"), (np.pi, "Nyquist")):
    H = np.polyval(bC[::-1], np.exp(-1j * wf)) / np.polyval(aC[::-1], np.exp(-1j * wf))
    P(f"10.3 DESAFIO |H| em {nome:10s} = {abs(H):.10e} ({20*np.log10(abs(H)):9.4f} dB)")

# ===========================================================================
# TOPICO 11 - DFT, leakage, janelas
# ===========================================================================
sec("11", "DFT, leakage e janelas")
N = 64
n = np.arange(N)
x_rect = np.ones(N)
X = np.fft.fft(x_rect)
P(f"11.1 retangular N={N}: |X[0]|={abs(X[0]):.4f}, |X[1]|={abs(X[1]):.1e},"
  f" |X[N/2]|={abs(X[N//2]):.3e}")
# senoide com 1.5 ciclos -> leakage
nc = 1.5
xc = np.cos(2 * np.pi * nc * n / N)
Xc = np.abs(np.fft.fft(xc))
k = np.argmax(Xc)
P(f"11.1b senoide de {nc} ciclos em N={N}: pico no bin k={k}"
  f" (freq. real {nc} ciclos); |X[k]|={Xc[k]:.4f}, |X[k-1]|={Xc[k-1]:.4f},"
  f" |X[k+1]|={Xc[k+1]:.4f}  -> leakage")
# resolucao
P("11.2 Fs=1000 Hz, N=1000 -> df = 1 Hz; dois tons com 10 Hz de separacao precisam"
  " df <= 10 Hz -> N >= 100")
# janelas: nivel do primeiro lobo lateral e largura do lobo principal (L=65)
L = 65
for nome, win in (("retangular", np.ones(L)),
                  ("Hamming", np.hamming(L)),
                  ("Hann", np.hanning(L)),
                  ("Blackman", np.blackman(L))):
    NF = 65536
    W = np.abs(np.fft.fft(win, NF))
    W /= W[0]
    d = 20 * np.log10(W + 1e-15)
    i = 1
    while d[i + 1] <= d[i]:
        i += 1
    prim_min = i                      # primeiro nulo (em bins do FFT de NF pontos)
    largura_ciclos = 2 * prim_min / NF           # largura do lobo em ciclos/amostra
    largura_bins_L = largura_ciclos * L          # em bins de uma DFT de L pontos
    pico_lat = d[prim_min:i + 2000].max()
    P(f"    {nome:11s}: 1o lobo lateral = {pico_lat:7.2f} dB,"
      f" largura do lobo principal = {largura_ciclos:.5f} ciclos/amostra"
      f" = {largura_bins_L:4.2f} bins de uma DFT de {L} pontos")

# ===========================================================================
# TOPICO 12 - fase, atraso de grupo, passa-tudo, fase minima
# ===========================================================================
sec("12", "Fase, atraso de grupo, passa-tudo, fase minima")


def atraso_grupo_1a_ordem(b, a, w):
    """
    tau(w) de H(z) = (1-b z^-1)/(1-a z^-1), em amostras.

    Dedução (a usada no material):
      H(e^jw) = (e^jw - b)/(e^jw - a)
      arg H  = arg(e^jw-b) - arg(e^jw-a)
      d/dw arg(e^jw - p) = (1 - p cos w)/(1 - 2p cos w + p^2)
      tau(w) = (1-a cos w)/(1-2a cos w+a^2) - (1-b cos w)/(1-2b cos w+b^2)
    """
    return ((1 - a * np.cos(w)) / (1 - 2 * a * np.cos(w) + a ** 2)
            - (1 - b * np.cos(w)) / (1 - 2 * b * np.cos(w) + b ** 2))


def atraso_grupo_passa_tudo(a, w):
    """tau(w) = (1-a^2)/(1-2a cos w+a^2)  [sempre > 0 para |a|<1]."""
    return (1 - a ** 2) / (1 - 2 * a * np.cos(w) + a ** 2)


b = 0.4     # zero do canal Hc(z) = (1-0.4z^-1)/(1-0.85z^-1)
a = 0.85    # polo do canal


P("    conferencia numerica da derivada (diferenças finitas):")
for w in (0.0, np.pi / 2, np.pi, 0.2 * np.pi):
    z, h = np.exp(1j * w), 1e-8
    H1 = (1 - b / z) / (1 - a / z)
    H2 = (1 - b / np.exp(1j * (w + h))) / (1 - a / np.exp(1j * (w + h)))
    num = -np.unwrap(np.angle([H1, H2]))[1] - (-np.unwrap(np.angle([H1, H2]))[0])
    num = -((-np.unwrap(np.angle([H1, H2]))[1]) - (-np.unwrap(np.angle([H1, H2]))[0]))
    num = -((np.unwrap(np.angle([H1, H2]))[1] - np.unwrap(np.angle([H1, H2]))[0])) / h
    P(f"      w={w/np.pi:.2f}pi: analitico = {atraso_grupo_1a_ordem(b,a,w):.8f},"
      f" numerico = {num:.8f}")
for w in (0.0, np.pi / 2, np.pi):
    P(f"12.1 Hc=(1-{b}z^-1)/(1-{a}z^-1): tau({w/np.pi:.2f}pi) = "
      f"{atraso_grupo_1a_ordem(b, a, w):.8f} amostras")
    P(f"     passa-tudo(a={a}) em omega={w/np.pi:.2f}pi: "
      f"{atraso_grupo_passa_tudo(a, w):.8f} amostras")

# 12.2 equalizador: alvo ALCANCAVEL -> tau_tot(0) = 10 amostras
alvo = 10.0
tau_c0 = atraso_grupo_1a_ordem(b, a, 0.0)
a_eq = ((alvo - tau_c0) - 1) / ((alvo - tau_c0) + 1)
P(f"12.2 alvo tau_tot(0) = {alvo}: tau_c(0) = {tau_c0:.6f} -> "
  f"(1+a)/(1-a) = {alvo - tau_c0:.6f} -> a = {a_eq:.10f} (= 2/3)")
P(f"     conferencia: tau_tot(0) = {tau_c0 + atraso_grupo_passa_tudo(a_eq, 0.0):.10f}")
# e o caso impossivel citado pelo desafio: 3.5 amostras em pi/2
tau_c_pi2 = atraso_grupo_1a_ordem(b, a, np.pi / 2)
max_eq_pi2 = max(atraso_grupo_passa_tudo(x, np.pi / 2) for x in np.linspace(-0.999, 0.999, 4001))
P(f"     em w=pi/2: tau_c = {tau_c_pi2:.6f}; maximo alcancavel do passa-tudo = "
  f"{max_eq_pi2:.6f} -> maximo total = {tau_c_pi2 + max_eq_pi2:.6f} < 3.5 (IMPOSSIVEL)")

# 12.3 decomposicao fase minima / passa-tudo  (H = 2 * Hmin * Hap)
# H(z) = (1-2z^-1)(1+0.3z^-1) / ((1-0.7z^-1)(1+0.6z^-1))
b3 = np.array([1.0, -1.7, -0.6])            # (1-2z^-1)(1+0.3z^-1) = 1 -1.7z^-1 -0.6z^-2
a3 = np.array([1.0, -0.1, -0.42])           # (1-0.7z^-1)(1+0.6z^-1)
P(f"12.3 H(z): zeros = {np.round(np.roots(b3),6)} (2 e -0.3),"
  f" polos = {np.round(np.roots(a3),6)} (0.7 e -0.6)")
h = signal.lfilter(b3, a3, np.eye(1, 8, 0).ravel())
P(f"     h[0..3] = {np.round(h[:4],10)}   soma(h) = {h.sum():.10f} (H(1) = -4)")
bmin = np.polymul(np.array([2.0]), np.polymul([1.0, -0.5], [1.0, 0.3]))  # K=2, zeros 0.5 e -0.3
hap = np.array([0.5, -1.0])                                  # Hap = (0.5-z^-1)/(1-0.5z^-1)
hmin = signal.lfilter(bmin, a3, np.eye(1, 8, 0).ravel())
hap_imp = signal.lfilter(hap, np.array([1.0, -0.5]), np.eye(1, 8, 0).ravel())
prod = np.convolve(hmin, hap_imp)
P(f"     Hmin (K=2): b={bmin}")
P(f"     hmin[0..3] = {np.round(hmin[:4],10)}")
P(f"     Hap: b={hap}, a=[1,-0.5]; |Hap| = "
  f"{np.round([abs(np.polyval(hap[::-1], np.exp(-1j*w))/np.polyval([1,-0.5][::-1], np.exp(-1j*w))) for w in np.linspace(0,np.pi,5)],8)}")
P(f"     h == conv(2*hmin, hap)? erro = {np.max(np.abs(prod[:8] - h)):.3e}")
P(f"     E3(H) = {np.sum(h[:3]**2):.10f}   E3(Hmin) = {np.sum(hmin[:3]**2):.10f}")
E_tot_h = np.sum(signal.lfilter(b3, a3, np.eye(1, 4000, 0).ravel()) ** 2)
E_tot_m = np.sum(signal.lfilter(bmin, a3, np.eye(1, 4000, 0).ravel()) ** 2)
P(f"     energia total (Parseval): H = {E_tot_h:.8f} , Hmin = {E_tot_m:.8f}"
  f"  -> iguais, pois |Hmin| = |H| (Hap e passa-tudo)")
P(f"     fracao da energia nos 3 primeiros: H = {np.sum(h[:3]**2)/E_tot_h:.6f},"
  f" Hmin = {np.sum(hmin[:3]**2)/E_tot_m:.6f}")

# ===========================================================================
# TOPICO 13 - FIR: janelamento, tipos, Kaiser
# ===========================================================================
sec("13", "FIR: janelamento e tipos")
wc = 0.3 * np.pi
M = 20
n = np.arange(M + 1)
hd = np.sin(wc * (n - M / 2)) / (np.pi * (n - M / 2))
hd[M // 2] = wc / np.pi
w_ham = np.hamming(M + 1)
hh = hd * w_ham
P(f"13.1 passa-baixas Hamming wc=0.3pi, M=20 (21 coeficientes), atraso M/2 = {M/2:.1f}")
P(f"     h[0..5] = {np.round(hh[:6],8)}")
P(f"     h[10] (centro) = {hh[10]:.8f} ; soma = {hh.sum():.8f} (ganho DC)")
P(f"     soma alternada (ganho em pi) = {hh.dot((-1.0)**n):.3e}")
# 13.2 Tipo II (M impar, simetrico) -> H(pi) = 0  (NAO ha amostra no centro)
M2 = 21
n2 = np.arange(M2 + 1)
hd2 = np.where(n2 - M2 / 2 == 0, wc,
               np.sin(wc * (n2 - M2 / 2)) / (np.pi * (n2 - M2 / 2)))
h2 = hd2 * np.hamming(M2 + 1)
P(f"13.2 Tipo II (M=21, simetrico, h[n]=h[M-n]): h[0..2]={np.round(h2[:3],8)}"
  f" h[10]={h2[10]:.8f} h[11]={h2[11]:.8f} simetria h[10]-h[11]={h2[10]-h2[11]:.2e}")
P(f"     soma alternada (= H(pi)) = {h2.dot((-1.0)**n2):.3e} ->"
  f" ZERO exato: Tipo II nao serve para passa-altas")
# 13.3 Kaiser: Fs=1000, fp=100, fs=200, As=50 dB
Fs, fp, fs, As = 1000.0, 100.0, 200.0, 50.0
dw = 2 * np.pi * (fs - fp) / Fs
beta = 0.1102 * (As - 8.7)
Mk = (As - 8) / (2.285 * dw)
nc_k = 100 * fp / (Fs / 2)
P(f"13.3 Kaiser: dw={dw:.10f} rad/amostra, beta={beta:.6f}, M~{Mk:.6f} -> M={int(np.ceil(Mk))}"
  f" (ordem), N={int(np.ceil(Mk))+1} taps")
P(f"     atraso = M/2 = {int(np.ceil(Mk))/2:.1f} amostras")
Mk = int(np.ceil(Mk))
nk = np.arange(Mk + 1)
hd_k = np.sin(0.35 * (nk - Mk / 2)) / (np.pi * (nk - Mk / 2))
hd_k[Mk // 2] = 0.35 / np.pi
wk = np.kaiser(Mk + 1, beta)
hk = hd_k * wk
P(f"     h[0..3] = {np.round(hk[:4],8)} ; ganho DC = {hk.sum():.6f}")

# ===========================================================================
# TOPICO 14 - quantizacao de coeficientes (Desafio 06)
# ===========================================================================
sec("14", "Quantizacao de coeficientes de IIR")


def quant12(x, lim=8.0, bits=12):
    """Quantiza para ponto fixo de 12 bits na faixa [-lim, +lim] (truncamento)."""
    q = 2 * lim / 2 ** bits
    return np.trunc(x / q) * q


polos = np.array([0.92 * np.exp(1j * 0.1 * np.pi), 0.92 * np.exp(-1j * 0.1 * np.pi),
                  0.95 * np.exp(1j * 0.25 * np.pi), 0.95 * np.exp(-1j * 0.25 * np.pi),
                  0.98 * np.exp(1j * 0.4 * np.pi), 0.98 * np.exp(-1j * 0.4 * np.pi)])
aden = np.real(np.poly(polos))
P(f"14.1 polos ideais: {np.round(polos,6)}")
P(f"     denominador (forma direta): a1..a6 = {np.round(aden[1:],12)}")
P(f"     passo de quantizacao (12 bits, [-8,8]) = {2*8/2**12:.10f}")
aden_q = quant12(aden)
P(f"     quantizado: {np.round(aden_q[1:],12)}")
rq = np.roots(aden_q)
P(f"     raizes quantizadas: {np.round(rq,8)}")
P(f"     |z| = {np.round(np.abs(rq),8)} -> max = {np.abs(rq).max():.8f}"
  f"  INSTAVEL? {np.abs(rq).max() >= 1}")
# em cascata (3 secoes de 2a ordem)
secoes = [(0.92, 0.1 * np.pi), (0.95, 0.25 * np.pi), (0.98, 0.4 * np.pi)]
P("14.3 cascata SOS (cada secao quantizada independentemente):")
maxmod = 0.0
for r, th in secoes:
    a1 = -2 * r * np.cos(th)
    a2 = r ** 2
    a1q, a2q = quant12(np.array([a1]))[0], quant12(np.array([a2]))[0]
    pq = np.roots([1.0, a1q, a2q])
    maxmod = max(maxmod, np.abs(pq).max())
    P(f"     r={r:.2f} theta={th/np.pi:.2f}pi: a1={a1:.10f}->{a1q:.10f},"
      f" a2={a2:.10f}->{a2q:.10f}, polos |z|={np.round(np.abs(pq),8)}")
P(f"     max |polo| na cascata = {maxmod:.8f} (< 1 -> ESTAVEL)")

# 14.4 varredura de bits: quando cada estrutura fica instavel
P("14.4 varredura do numero de bits (faixa [-8,+8]) - max |polo|:")
P("     bits   passo        Forma Direta        Cascata (SOS)")
for bits in (6, 7, 8, 9, 10, 11, 12, 14, 16):
    q = 2 * 8 / 2 ** bits
    adq = np.trunc(aden / q) * q
    md = np.abs(np.roots(adq)).max()
    mc = 0.0
    for r, th in secoes:
        a1q = np.trunc(-2 * r * np.cos(th) / q) * q
        a2q = np.trunc(r ** 2 / q) * q
        mc = max(mc, np.abs(np.roots([1.0, a1q, a2q])).max())
    flag_d = "INSTAVEL" if md >= 1 else "estavel"
    flag_c = "INSTAVEL" if mc >= 1 else "estavel"
    P(f"     {bits:3d}   {q:10.8f}   {md:.6f} {flag_d:9s}  {mc:.6f} {flag_c}")
# bits minimos para cada estrutura ficar ESTAVEL de forma monotonica
def max_polo_direta(b):
    q = 2 * 8 / 2 ** b
    return float(np.abs(np.roots(np.trunc(aden / q) * q)).max())


def max_polo_cascata(b):
    q = 2 * 8 / 2 ** b
    return float(max(np.abs(np.roots([1.0,
                                      np.trunc(-2 * r * np.cos(th) / q) * q,
                                      np.trunc(r ** 2 / q) * q])).max()
                     for r, th in secoes))


for nome, fn in (("forma direta", max_polo_direta), ("cascata SOS", max_polo_cascata)):
    est = {b: fn(b) < 1 for b in range(5, 25)}
    limiar = next((b for b in range(5, 25) if all(est[k] for k in range(b, 25))), None)
    P(f"     {nome}: estavel para todo b >= {limiar} bits"
      f" (em {limiar - 1} bits: max|z| = {fn(limiar - 1):.6f})")

# ===========================================================================
# TOPICO 15 - invariancia ao impulso
# ===========================================================================
sec("15", "Invariancia ao impulso")
# 15.1 1a ordem
a = 1000.0
Ts = 1e-4
P(f"15.1 ha=e^-1000t u(t), Ts=0.1 ms: r=exp(-aTs)={np.exp(-a*Ts):.10f},"
  f" H(z)={Ts:.1e}/(1-{np.exp(-a*Ts):.6f}z^-1); ganho DC = {Ts/(1-np.exp(-a*Ts)):.10f}"
  f" (analitico 1/a = {1/a:.10f})")
# 15.2 2a ordem do desafio
P(f"15.2 DESAFIO: r={np.exp(-sigC*25e-6):.12f}, theta={wdC*25e-6:.12f} rad"
  f" = {wdC*25e-6/np.pi:.12f}pi")
b0 = 25e-6
b1 = -25e-6 * np.exp(-sigC * 25e-6) * (np.cos(wdC * 25e-6) + (sigC / wdC) * np.sin(wdC * 25e-6))
a1 = -2 * np.exp(-sigC * 25e-6) * np.cos(wdC * 25e-6)
a2 = np.exp(-2 * sigC * 25e-6)
P(f"     b0={b0:.12e} b1={b1:.12e} a1={a1:.12e} a2={a2:.12e}")
# 15.3 passa-altas degenerado
ahp = 2000 * np.pi
Ts = 25e-6
al = np.exp(-ahp * Ts)
P(f"15.3 passa-altas a=2000pi: alpha={al:.10f}, |H(1)| = {Ts - ahp*Ts/(1-al):.6f}"
  f" ({20*np.log10(abs(Ts - ahp*Ts/(1-al))):.4f} dB)")
# em z = -1: H = (b0 - b1)/(1 - a1)  com b0 = Ts - a*Ts, b1 = -Ts*alpha, a1 = -alpha
H_ny = (Ts - ahp * Ts + Ts * al) / (1 + al)
P(f"     |H(-1)| = {abs(H_ny):.6f} ({20*np.log10(abs(H_ny)):.4f} dB)"
  f"  [confere com desafio08.py: 0.084671 / -21.4453 dB]")

# ===========================================================================
# TOPICO 16 - bilinear e pre-warping
# ===========================================================================
sec("16", "Transformacao bilinear e pre-deformacao")
Fs = 1000.0
for f in (5.0, 20.0, 150.0, 300.0, 100.0):
    Om = 2 * Fs * np.tan(np.pi * f / Fs)
    P(f"16.1 f={f:6.1f} Hz -> Omega = 2*Fs*tan(pi*f/Fs) = {Om:.6f} rad/s"
      f"   (sem warp: 2*pi*f = {2*np.pi*f:.6f})")
# 16.2 1a ordem: Ha(s) = Om_p/(s+Om_p), Om_p = 2Fs tan(pi*100/Fs)
Omp = 2 * Fs * np.tan(np.pi * 100 / Fs)
b = np.array([Omp, Omp])
a = np.array([2 * Fs + Omp, Omp - 2 * Fs])
b = b / a[0]
a = a / a[0]
P(f"16.2 bilinear 1a ordem (fc=100 Hz, Fs=1000): b={np.round(b,10)}, a={np.round(a,10)}")
for f, nome in ((0.0, "DC"), (100.0, "-3dB nominal"), (500.0, "Nyquist")):
    w = 2 * np.pi * f / Fs
    H = np.polyval(b[::-1], np.exp(-1j * w)) / np.polyval(a[::-1], np.exp(-1j * w))
    P(f"     |H| em {nome:14s} = {abs(H):.10f} ({20*np.log10(abs(H)):9.4f} dB)")
# 16.3 ordem Butterworth passa-faixa (Desafio 09)
Fp1, Fp2, Fs1, Fs2, Fs = 20.0, 150.0, 5.0, 300.0, 1000.0
Ap, As = 1.0, 30.0
Sp1 = 2 * Fs * np.tan(np.pi * Fp1 / Fs)
Sp2 = 2 * Fs * np.tan(np.pi * Fp2 / Fs)
Ss1 = 2 * Fs * np.tan(np.pi * Fs1 / Fs)
Ss2 = 2 * Fs * np.tan(np.pi * Fs2 / Fs)
Bp = (Sp2 - Sp1) / (Sp1 * Sp2)
Bs = (Ss2 - Ss1) / (Ss1 * Ss2)
P(f"16.3 pre-warp: Op1={Sp1:.4f}, Op2={Sp2:.4f}, Os1={Ss1:.4f}, Os2={Ss2:.4f} rad/s")
P(f"     Bp={Bp:.10f}, Bs={Bs:.10f}")
eps_p = np.sqrt(10 ** (0.1 * Ap) - 1)
eps_s = np.sqrt(10 ** (0.1 * As) - 1)
P(f"     eps_p={eps_p:.6f}, eps_s={eps_s:.6f}")
r1 = Bs / Bp
r2 = (Ss1 / Sp1) * (Sp2 / Ss1)  # razao de frequencias para o prototipo LP equivalente
P(f"     Bs/Bp = {r1:.6f}")
N_lp = np.log10(eps_s / eps_p) / np.log10(r1)
P(f"     N >= log10(eps_s/eps_p)/log10(Bs/Bp) = {N_lp:.6f} -> N = {int(np.ceil(N_lp))}"
  f" -> para passa-faixa: 2N = {2*int(np.ceil(N_lp))} pólos")
# ordem via prototipo passa-baixas equivalente (mascara geometrica)
Om0 = np.sqrt(Sp1 * Sp2)
P(f"     centro geometrico do passa-faixa Om0 = {Om0:.6f} rad/s"
  f" ({Om0/2/np.pi:.4f} Hz)")
B = Sp2 - Sp1
Om_stop = []
for Os in (Ss1, Ss2):
    v = abs(Os ** 2 - Om0 ** 2) / (B * Os)
    Om_stop.append(v)
P(f"     Om_stop equivalentes = {np.round(Om_stop,6)} rad/s"
  f" (usar o MENOR)")
Om_pb = 1.0
N_pb = np.log10(eps_s / eps_p) / np.log10(max(Om_stop))
P(f"     N >= log10(eps_s/eps_p)/log10(Omega_s,min) = {N_pb:.6f} -> N={int(np.ceil(N_pb))}")

P("\n" + "=" * 78)
P("FIM DO GABARITO")
P("=" * 78)
