#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Autor: Matheus Sousa, LCSI-DES, UFPE.
"""
===============================================================================
 PDS 2026.2 - DUPLA: DESAFIO 08
 Discretizacao de Filtros IIR por Invariancia ao Impulso
 e Analise de Aliasing Espectral
===============================================================================

 SCRIPT PRINCIPAL - rode este arquivo:

     python src/desafio08.py

 Ele faz, em ordem:
   1) a resolucao analitica/numerica dos itens 1, 2 e 3 do enunciado;
   2) a resolucao do Desafio de Reflexao da Dupla (falha no passa-altas);
   3) a Parte Computacional completa -> 7 figuras em figuras/ (.png e .pdf);
   4) o relatorio numerico em resultados/relatorio_resultados.md.

 Dependencias: numpy, matplotlib, sympy  (ver README.md).
===============================================================================
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator

RAIZ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))

import pds08_core as core   # noqa: E402

FIG = RAIZ / "figuras"
RES = RAIZ / "resultados"
FIG.mkdir(exist_ok=True)
RES.mkdir(exist_ok=True)

PI = np.pi
COR_ANALOG = "#1f4e79"      # azul escuro  - filtro analogico
COR_DIGITAL = "#c00000"     # vermelho     - filtro digital
COR_REPLICA = "#7f7f7f"     # cinza        - replicas / aliasing
COR_EXTRA = "#2e8b57"       # verde        - casos extras
COR_HL = "#e8a33d"          # laranja      - destaques/marcadores

CABECALHO_LOG: list[str] = []


# ---------------------------------------------------------------------------
# utilidades de saida
# ---------------------------------------------------------------------------
def log(texto: str = "") -> None:
    """Imprime na tela e guarda o texto para o relatorio .md."""
    print(texto)
    CABECALHO_LOG.append(texto)


def titulo(n: int, texto: str) -> None:
    log()
    log("=" * 78)
    log(f"({n}) {texto}")
    log("=" * 78)


def estilo() -> None:
    plt.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": 160,
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "axes.grid": True,
        "grid.alpha": 0.30,
        "grid.linestyle": ":",
        "legend.fontsize": 8.5,
        "legend.framealpha": 0.92,
        "figure.autolayout": False,
        "axes.axisbelow": True,
        "mathtext.fontset": "dejavusans",
    })


def salva(fig, nome: str) -> None:
    """Salva a figura em PNG (slides) e PDF (vetorial)."""
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{nome}.{ext}", bbox_inches="tight")
    plt.close(fig)
    log(f"    [figura] figuras/{nome}.png | figuras/{nome}.pdf")


def dB(x):
    """Magnitude -> dB, com piso para não gerar -inf nos graficos."""
    x = np.asarray(x, dtype=float)
    return 20.0 * np.log10(np.maximum(x, 1e-12))


# ===========================================================================
# ITEM 1 - ha(t) = L^-1{Ha(s)}
# ===========================================================================
def item1(p) -> None:
    titulo(1, "ITEM 1 - Resposta ao impulso analógica ha(t) = L^-1{Ha(s)}")

    log("  Ha(s) = s / [ s^2 + 200*pi*s + (12000*pi)^2 ]")
    log("        = s / [ s^2 + 2*zeta*wn*s + wn^2 ]")
    log("")
    log(f"  wn     = 12000*pi        = {p.wn:.6f} rad/s      "
        f"(fn = {p.fn:.4f} Hz)")
    log(f"  2*zeta*wn = 200*pi       => zeta = 200/(2*12000) = 1/120 = "
        f"{p.zeta:.10f}")
    log(f"  Q = 1/(2*zeta)           = {p.Q:.6f}")
    log(f"  sigma  = zeta*wn = 100*pi = {p.sigma:.6f} rad/s")
    log(f"  wd = wn*sqrt(1-zeta^2)   = 100*pi*sqrt(14399) = {p.wd:.6f} rad/s")
    log("")
    log("  Fracoes parciais:")
    log("        s                (s+sigma)             sigma")
    log("  ---------------  =  ----------------  -  ----------------")
    log("  (s+sigma)^2+wd^2   (s+sigma)^2+wd^2     (s+sigma)^2+wd^2")
    log("")
    log("  Logo:")
    log("    ha(t) = e^(-100*pi*t) * [ cos(wd*t) - (100/11999.583)*sen(wd*t) ], t >= 0")
    log("    ha(t) = e^(-100*pi*t) * [ cos(wd*t) - (1/sqrt(14399))*sen(wd*t) ]")
    log("")
    log(f"  com wd = 100*pi*sqrt(14399) = {p.wd:.6f} rad/s")
    log(f"  Forma exata (SymPy): (-sqrt(119)*sen(1100*sqrt(119)*pi*t)/1309 + "
        f"cos(1100*sqrt(119)*pi*t))*exp(-100*pi*t)")
    log("")
    log(f"  ha(0+) = {core.ha_analogica(0.0):.6f}   (salto em t = 0 -> relevante "
        f"para o item de aliasing)")
    log(f"  |Ha(j*wn)| = 1/(2*zeta*wn) = 1/(200*pi) = "
        f"{1.0 / core.TWO_ZETA_WN:.10e}")


# ===========================================================================
# ITEM 2 - H(z) por invariancia ao impulso
# ===========================================================================
def item2(f: core.FiltroDigital, p) -> None:
    titulo(2, "ITEM 2 - Discretizacao por invariância ao impulso: coeficientes de H(z)")

    log(f"  Ts = 25 us  (Fs = 40 kHz)")
    log(f"  h[n] = Ts*ha(n*Ts)")
    log("")
    log("  Deduzindo com  r = exp(-sigma*Ts)  e  theta = wd*Ts :")
    log()
    log("         Ts*[ 1 - r*(cos(theta) + (sigma/wd)*sen(theta)) * z^-1 ]")
    log("  H(z) = ---------------------------------------------------------")
    log("                1 - 2*r*cos(theta)*z^-1 + r^2*z^-2")
    log()
    log(f"  r     = exp(-100*pi*25e-6) = {f.r:.12f}")
    log(f"  theta = wd*Ts              = {f.theta:.12f} rad "
        f"= {f.theta / PI:.12f}*pi rad/amostra")
    log(f"  cos(theta) = {np.cos(f.theta):.12f}   sen(theta) = {np.sin(f.theta):.12f}")
    log("")
    log("  COEFICIENTES NUMERICOS EXATOS:")
    log(f"    b0 = Ts                 = {f.b[0]:.12e}")
    log(f"    b1 = -Ts*r*(cos+ (sig/wd)sen) = {f.b[1]:.12e}")
    log(f"    a1 = -2*r*cos(theta)    = {f.a[1]:.12e}")
    log(f"    a2 = r^2                = {f.a[2]:.12e}")
    log("")
    log(f"  H(z) = ({f.b[0]:.6e} {f.b[1]:+.6e} z^-1) / "
        f"(1 {f.a[1]:+.6e} z^-1 {f.a[2]:+.6e} z^-2)")
    log("")
    polos = f.polos
    z = f.zeros[0]
    log(f"  Polos (complexos conjugados): z = {np.abs(polos[0]):.9f} * "
        f"exp(+/-j*{np.abs(np.angle(polos[0])):.9f})  ->  |z| < 1  (ESTAVEL)")
    log(f"     |z| = r = exp(-sigma*Ts) = {f.r:.9f}  (|z| < 1 confirma estabilidade)")
    log(f"     angulo = theta = {f.theta / PI:.9f}*pi rad/amostra  "
        f"-> f = {f.theta / (2 * PI * f.Ts):.6f} Hz")
    log(f"  Zero: z = -b1/b0 = {z:.9f}  (inventado pelo metodo; zeros NAO mapeiam "
        f"por z = e^(s*Ts))")
    log("")
    log(f"  Ganho em DC: |H(e^j0)| = {np.abs(f.reais_em(0.0)).item():.6e} "
        f"(~ Ts/2 = {f.Ts / 2:.6e}, residuo do termo de meia-descontinuidade)")
    log(f"  Frequência pedida: w0 = 2*pi*6000*Ts = {f.w0_raio / PI:.10f}*pi "
        f"= 0.3*pi rad/amostra (6000 Hz)")


# ===========================================================================
# ITEM 3 - ganho de pico e aliasing
# ===========================================================================
def item3(f, p) -> dict:
    titulo(3, "ITEM 3 - Ganho de pico em w0 = 0.3*pi e analise de aliasing")

    pa = core.pico_analogico()
    pd = core.pico_digital(f)
    banda = core.banda_meia_potencia(f)

    log("  (a) Filtro ANALOGICO")
    log(f"      pico |Ha| = {pa['ganho']:.10e}  = {pa['ganho_dB']:.4f} dB")
    log(f"      em Omega_r = wn*sqrt(1-2*zeta^2) = {pa['Omega_r']:.4f} rad/s "
        f"= {pa['f_r']:.4f} Hz")
    log(f"      em Omega = wn: |Ha(j*wn)| = {pa['ganho_em_wn']:.10e} "
        f"= 1/(2*zeta*wn) = 1/(200*pi)")
    log("")
    log("  (b) Filtro DIGITAL obtido por invariância ao impulso (Fs = 40 kHz)")
    log(f"      |H(e^j*w0)| com w0 = 0.3*pi = {pd['ganho_em_w0']:.10e} "
        f"= {pd['ganho_em_w0_dB']:.4f} dB")
    log(f"      pico varrido  = {pd['ganho_pico']:.10e} = {pd['ganho_pico_dB']:.4f} dB "
        f"em w = {pd['w_pico'] / PI:.8f}*pi  (f = {pd['f_pico']:.4f} Hz)")
    log(f"      banda de -3 dB: {banda['df_3dB_Hz']:.4f} Hz  "
        f"(w de {banda['w_3dB_baixo'] / PI:.6f}*pi a {banda['w_3dB_alto'] / PI:.6f}*pi)")
    log(f"      ganho em DC   = {pd['ganho_DC']:.6e} (~0, coerente com Ha(j0) = 0)")
    log("")
    dif_db = pd["ganho_em_w0_dB"] - pa["ganho_dB"]
    dif_pct = 100.0 * (pd["ganho_em_w0"] / pa["ganho_em_wn"] - 1.0)
    log("  (c) Comparacao analógico x digital em 6000 Hz")
    log(f"      analógico: {pa['ganho_em_wn']:.10e} "
        f"({20 * np.log10(pa['ganho_em_wn']):.4f} dB)")
    log(f"      digital  : {pd['ganho_em_w0']:.10e} ({pd['ganho_em_w0_dB']:.4f} dB)")
    log(f"      diferença: {dif_db:+.4f} dB  ({dif_pct:+.4f} %)")
    log("")
    log("  RESPOSTA: a diferença e de ~0,07 dB (0,8 %), ou seja, DESPREZIVEL.")
    log("  Não ha alteracao significativa do pico: as réplicas do espectro caem em")
    log("  34 kHz e 46 kHz (= Fs -/+ 6 kHz), onde |Ha| e cerca de 50 dB abaixo do")
    log("  pico:")
    log(f"      |Ha(j*2*pi*34kHz)| = {core.ganho_analogico_em(2 * PI * 34000)[0]:.3e}"
        f" = {dB(core.ganho_analogico_em(2 * PI * 34000)[0]):.1f} dB absolutos"
        f" (pico = {20 * np.log10(pa['ganho_em_wn']):.1f} dB)")
    log(f"      |Ha(j*2*pi*46kHz)| = {core.ganho_analogico_em(2 * PI * 46000)[0]:.3e}"
        f" = {dB(core.ganho_analogico_em(2 * PI * 46000)[0]):.1f} dB absolutos")
    log("  A condição de Nyquist e folgadamente satisfeita: Fs/2 = 20 kHz >> 6 kHz.")
    log("")
    log("  (d) Contraprova - reducao de Fs (aliasing frequêncial):")
    for Fs in (20_000.0, 15_000.0, 12_000.0, 8_000.0):
        fx = core.discretiza_invariancia_impulso(Fs)
        px = core.pico_digital(fx, n_pontos=200_001)
        razao = px["ganho_pico"] / pa["ganho_em_wn"]
        log(f"      Fs = {Fs / 1000:5.1f} kHz -> pico {px['ganho_pico_dB']:8.3f} dB "
            f"em {px['f_pico']:8.2f} Hz   ({razao:6.3f}x o pico analógico)")
    log("")
    log("  Com Fs = 15 kHz a imagem da componente de -6 kHz (em 9 kHz) dobra")
    log("  exatamente sobre 6 kHz -> pico sobe de +0,19 dB. Com Fs = 12 kHz a")
    log("  ressonância cai em Nyquist e a imagem se soma EM FASE: o pico dobra")
    log("  (+6,1 dB). Abaixo de 12 kHz a própria ressonância e subamostrada.")
    return {"analógico": pa, "digital": pd, "banda": banda}


# ===========================================================================
# DESAFIO DE REFLEXAO - falha do metodo no passa-altas
# ===========================================================================
def reflexao() -> dict:
    titulo(4, "DESAFIO DE REFLEXAO DA DUPLA - Falha no filtro PASSA-ALTAS")

    fhp = core.filtro_passa_altas_invariancia()
    a = core.HP_A_CONST
    alpha = np.exp(-a * fhp.Ts)
    gDC = np.abs(fhp.reais_em(0.0)).item()
    gNy = np.abs(fhp.reais_em(PI)).item()

    log("  Ha_hp(s) = s/(s + 2000*pi)   ->  fc = 1 kHz, |Ha(j*0)| = 0 e")
    log("  |Ha(j*Omega)| -> 1 quando Omega -> infinito (passa-altas ideal em HF).")
    log("")
    log("  ha_hp(t) = delta(t) - 2000*pi*e^(-2000*pi*t)*u(t)")
    log("")
    log("  Amostrando:  h[n] = Ts*delta[n] - Ts*2000*pi*alpha^n*u[n],")
    log(f"               com alpha = exp(-2000*pi*Ts) = {alpha:.9f}")
    log("")
    log(f"  H(z) = Ts - Ts*2000*pi/(1 - alpha*z^-1)")
    log(f"       = [ (Ts - Ts*2000*pi) - Ts*alpha*z^-1 ] / (1 - alpha*z^-1)")
    log("")
    log(f"    b = [{fhp.b[0]:.10e}, {fhp.b[1]:.10e}]")
    log(f"    a = [1, {fhp.a[1]:.10e}]")
    log(f"    zero em z = -b1/b0 = {-fhp.b[1] / fhp.b[0]:.10e}  (praticamente z = 0)")
    log(f"    polo em z = alpha   = {alpha:.10f}")
    log("")
    log(f"  |H(e^j0)|   = {gDC:.6f}  ({20 * np.log10(gDC):+.4f} dB)  <-- DC")
    log(f"  |H(e^jpi)|  = {gNy:.6f}  ({20 * np.log10(gNy):+.4f} dB)  <-- Nyquist")
    log("")
    log("  CONCLUSAO: o ganho agora DECRESCE com a frequência (1,0806 -> 0,0847):")
    log("  o 'passa-altas' virou um PASSA-BAIXAS degenerado, com ganho > 1 em DC.")
    log("  O zero do eixo j*Omega em s = 0 (que definia o passa-altas) não tem")
    log("  nenhum zero correspondente em z = 1: os zeros não são mapeados por")
    log("  z = e^(s*Ts) pelo metodo da invariância ao impulso.")
    log("")
    log("  Causa de fundo: |Ha(j*Omega)| NAO tende a zero quando Omega -> oo.")
    log("  A identidade de Poisson")
    log("      H(e^jw) = sum_k Ha(j*(w-2*pi*k)/Ts) + Ts*ha(0+)/2")
    log("  tem, para o passa-altas, réplica de amplitude ~1 em TODO k, ou seja, o")
    log("  somatorio das réplicas NAO converge (energia de todo o eixo de frequência")
    log("  e dobrada para dentro da banda digital). Nenhuma taxa de amostragem")
    log("  corrige isso - o problema e do metodo, não de Fs.")
    log("")
    log("  RESSALVA PROVENIENTE DA MEDICAO NUMERICA (verificacao.py):")
    log(f"    |Ha(j*34kHz)| = {core.resposta_analogica_passa_altas(2 * PI * 34000)[0]:.3f} "
        f"-> a réplica em 34 kHz tem ganho praticamente 1.")
    return {"fhp": fhp, "gDC": gDC, "gNy": gNy, "alpha": alpha}


# ===========================================================================
# FIGURA 1 - ha(t) + amostras
# ===========================================================================
def fig01_ha_amostras(f, p) -> None:
    fig, ax = plt.subplots(1, 2, figsize=(11.2, 4.1))

    t = np.linspace(0.0, 250e-6, 4000)
    n = np.arange(0, 11)
    # ATENCAO: h[n]/Ts = ha(n*Ts) - e a amplitude que se sobrepoe a curva ha(t)
    hn_amostras = core.ha_analogica(n * f.Ts)

    ax[0].plot(1e6 * t, core.ha_analogica(t), color=COR_ANALOG, lw=1.8,
               label=r"$h_a(t)$ analitica")
    ax[0].stem(1e6 * n * f.Ts, hn_amostras, linefmt=":", markerfmt="o", basefmt=" ",
               label=r"amostras $h[n]/T_s = h_a(nT_s)$")
    ax[0].set_xlabel(r"tempo $t$ [$\mu$s]")
    ax[0].set_ylabel(r"$h_a(t)$  [1/s]")
    ax[0].set_title(r"(a) $h_a(t)$ e amostras ($T_s=25\ \mu$s) - detalhe inicial")
    ax[0].set_xlim(-5, 250)
    ax[0].legend(loc="upper right")
    ax[0].annotate(r"$h[0]/T_s = h_a(0^+) = 1$", xy=(0, 1.0),
                   xytext=(55, 0.80), fontsize=8.5,
                   arrowprops=dict(arrowstyle="->", lw=0.9, color="k"))

    t2 = np.linspace(0.0, 1.0e-3, 6000)
    n2 = np.arange(0, 41)
    ax[1].plot(1e6 * t2, core.ha_analogica(t2), color=COR_ANALOG, lw=1.5,
               label=r"$h_a(t)$")
    ax[1].stem(1e6 * n2 * f.Ts, core.ha_analogica(n2 * f.Ts),
               linefmt=":", markerfmt="o", basefmt=" ",
               label=r"$h[n]/T_s = h_a(nT_s)$")
    ax[1].plot(1e6 * t2, np.exp(-p.sigma * t2), color=COR_HL, ls="--", lw=1.3,
               label=r"envoltória $\pm e^{-\sigma t}$")
    ax[1].plot(1e6 * t2, -np.exp(-p.sigma * t2), color=COR_HL, ls="--", lw=1.3)
    ax[1].set_xlabel(r"tempo $t$ [$\mu$s]")
    ax[1].set_ylabel(r"$h_a(t)$  [1/s]")
    ax[1].set_title(r"(b) $h_a(t)$ oscila em $f_d \approx 5999{,}8$ Hz (Q = 60)")
    ax[1].legend(loc="upper right")

    fig.suptitle("Desafio 08 - Item 1/Parte Computacional: $h_a(t)$ analitica vs. "
                 r"amostras $h[n]/T_s = h_a(nT_s)$  ($F_s=40$ kHz)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    salva(fig, "fig01_ha_t_e_amostras")


# ===========================================================================
# FIGURA 2 - resposta ao impulso digital
# ===========================================================================
def fig02_impulso_digital(f) -> None:
    n = np.arange(0, 200)
    hn = core.resposta_impulso_digital(n_amostras=200)

    fig, ax = plt.subplots(figsize=(9.2, 3.9))
    ax.stem(n, hn, linefmt="C0-", markerfmt="C0.", basefmt="k-")
    env = f.Ts * np.abs(core.ha_analogica(np.arange(250) * f.Ts))
    ax.plot(np.arange(250), env, color=COR_HL, ls="--", lw=1.3,
            label=r"envoltória $T_s e^{-\sigma n T_s}$")
    ax.plot(np.arange(250), -env, color=COR_HL, ls="--", lw=1.3)
    ax.set_xlabel(r"$n$ [amostras]")
    ax.set_ylabel(r"$h[n]$")
    ax.set_title(r"Resposta ao impulso do filtro digital $h[n]=T_s h_a(nT_s)$ "
                 r"($F_s=40$ kHz, $\omega_0=0{,}3\pi$ rad/amostra)")
    ax.set_xlim(-2, 202)
    ax.legend(loc="upper right")
    ax.annotate(r"$h[0]=T_s=25\ \mu$s", xy=(0, f.Ts), xytext=(25, 1.55e-5),
                arrowprops=dict(arrowstyle="->", lw=0.9, color="k"), fontsize=9)
    fig.tight_layout()
    salva(fig, "fig02_resposta_impulso_digital")


# ===========================================================================
# FIGURA 3 - magnitude analogico vs digital (40 kHz e 15 kHz)
# ===========================================================================
def fig03_magnitude(f40, dados3) -> None:
    pa = dados3["analógico"]
    fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.4))

    # ---- (a) Fs = 40 kHz
    ff = np.linspace(0.0, 21e3, 21001)
    ax[0].plot(ff / 1e3, dB(core.ganho_analogico_em(2 * PI * ff)),
               color=COR_ANALOG, lw=2.4, alpha=0.85, ls="-",
               label=r"analógico $|H_a(j\Omega)|$")
    w = np.linspace(0, PI, 8001)
    ax[0].plot(w * f40.Fs / (2 * PI) / 1e3, dB(np.abs(f40.reais_em(w))),
               color=COR_DIGITAL, lw=1.4, ls="--",
               label=r"digital $|H(e^{j\omega})|$, $F_s=40$ kHz")
    ax[0].axvline(6.0, color=COR_HL, ls=":", lw=1.4)
    ax[0].axvline(f40.Fs / 2e3, color="k", ls="-.", lw=1.2, alpha=0.7,
                  label=r"Nyquist $F_s/2=20$ kHz")
    ax[0].annotate("pico em 6 kHz\nanalógico: %.2f dB\ndigital:   %.2f dB"
                   % (20 * np.log10(pa["ganho_em_wn"]),
                      20 * np.log10(dados3["digital"]["ganho_em_w0"])),
                   xy=(6, 20 * np.log10(pa["ganho_em_wn"])),
                   xytext=(8.6, -32), fontsize=8.5,
                   arrowprops=dict(arrowstyle="->", lw=0.9, color="k"))
    ax[0].set_ylabel("magnitude  [dB]")
    ax[0].set_xlabel("frequência  [kHz]")
    ax[0].set_title(r"(a) $F_s = 40$ kHz: as duas curvas praticamente coincidem")
    ax[0].set_ylim(-130, -45)
    ax[0].set_xlim(0, 21)
    ax[0].legend(loc="lower right")

    # detalhe do pico: evidencia a diferenca de 0,07 dB
    wz = np.linspace(0.28 * PI, 0.32 * PI, 2001)
    ins = ax[0].inset_axes([0.10, 0.16, 0.42, 0.40])
    ins.plot(wz / PI, dB(np.abs(f40.reais_em(wz))), color=COR_DIGITAL,
             lw=1.6, ls="--", label="digital")
    ins.plot(wz / PI, dB(core.ganho_analogico_em(wz / f40.Ts)), color=COR_ANALOG,
             lw=1.6, label="analógico")
    ins.set_title("detalhe do pico (0,01 dB/div)", fontsize=7.5)
    ins.tick_params(labelsize=7)
    ins.set_xlabel(r"$\omega/\pi$", fontsize=7.5)
    ins.set_ylabel("dB", fontsize=7.5)
    ins.legend(fontsize=6.5, loc="lower right")
    ins.grid(alpha=0.3, ls=":")

    # ---- (b) Fs = 15 kHz
    f15 = core.discretiza_invariancia_impulso(15_000.0)
    ax[1].plot(ff / 1e3, dB(core.ganho_analogico_em(2 * PI * ff)),
               color=COR_ANALOG, lw=2.4, alpha=0.85, ls="-",
               label=r"analógico $|H_a(j\Omega)|$")
    w15 = np.linspace(0, PI, 8001)
    ax[1].plot(w15 * f15.Fs / (2 * PI) / 1e3, dB(np.abs(f15.reais_em(w15))),
               color=COR_DIGITAL, lw=1.6, ls="--",
               label=r"digital $|H(e^{j\omega})|$, $F_s=15$ kHz")
    ax[1].axvline(6.0, color=COR_HL, ls=":", lw=1.4)
    ax[1].axvline(9.0, color=COR_REPLICA, ls="-.", lw=1.4,
                  label="imagem de -6 kHz: 15-6 = 9 kHz")
    ax[1].axvline(f15.Fs / 2e3, color="k", ls="-.", lw=1.2, alpha=0.7,
                  label=r"Nyquist $F_s/2=7{,}5$ kHz")
    ax[1].annotate("", xy=(6.0, -60), xytext=(9.0, -60),
                   arrowprops=dict(arrowstyle="->", lw=1.5, color=COR_REPLICA))
    ax[1].text(7.5, -58.5, "dobra sobre 6 kHz", fontsize=8.5,
               color=COR_REPLICA, ha="center")
    ax[1].annotate("pico: 0,8$\\pi$ = 6 kHz\n+0,19 dB por aliasing",
                   xy=(6, 20 * np.log10(pa["ganho_em_wn"])),
                   xytext=(0.8, -80), fontsize=8.5,
                   arrowprops=dict(arrowstyle="->", lw=0.9, color="k"))
    ax[1].set_title(r"(b) $F_s = 15$ kHz: imagem em 9 kHz dobra para 6 kHz")
    ax[1].set_xlabel("frequência  [kHz]")
    ax[1].set_ylabel("magnitude  [dB]")
    ax[1].set_ylim(-130, -45)
    ax[1].set_xlim(0, 10)
    ax[1].legend(loc="lower right", fontsize=7.5)

    fig.suptitle("Desafio 08 - Item 3 / Parte Computacional: resposta em magnitude "
                 "do protótipo analógico e do filtro digital IIR (invariância ao impulso)",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    salva(fig, "fig03_magnitude_analogico_vs_digital")


# ===========================================================================
# FIGURA 4 - mecanismo de aliasing: replicas do espectro
# ===========================================================================
def fig04_replicas() -> None:
    figs = [(40_000.0, "(a) $F_s = 40$ kHz: réplicas em 34 kHz e 46 kHz"),
            (15_000.0, "(b) $F_s = 15$ kHz: réplica em 9 kHz invade a banda"),
            (8_000.0, "(c) $F_s = 8$ kHz: a ressonância de 6 kHz dobra para 2 kHz")]
    fig, ax = plt.subplots(1, 3, figsize=(13.6, 4.2))

    for a, (Fs, ttl) in zip(ax, figs):
        f = core.discretiza_invariancia_impulso(Fs)
        base = np.linspace(0.0, 1.6 * Fs, 40001)
        for k in range(-2, 3):
            a.plot(base / 1e3, dB(core.ganho_analogico_em(2 * PI * (base - k * Fs))),
                   color=COR_REPLICA if k else COR_ANALOG, lw=1.9 if k == 0 else 1.0,
                   ls="-" if k == 0 else "--", alpha=0.95 if k == 0 else 0.70,
                   label=(r"analógico $|H_a(j2\pi f)|$" if k == 0 else
                          r"réplicas $|H_a(f-kF_s)|$") if k in (0, 1) else None)
        a.axvspan(0, Fs / 2e3, color=COR_EXTRA, alpha=0.10,
                  label=r"banda-base digital $0\ldots F_s/2$")
        wg = np.linspace(0, PI, 4001)
        a.plot(wg * Fs / (2 * PI) / 1e3, dB(np.abs(f.reais_em(wg))),
               color="k", lw=1.8, label=r"digital exato $|H(e^{j\omega})|$")
        a.axvline(Fs / 2e3, color="k", ls=":", lw=1.1)
        a.text(Fs / 2e3, -44, r"$F_s/2$", fontsize=8, ha="center", va="top")
        pk = np.argmax(np.abs(f.reais_em(wg)))
        a.plot(wg[pk] * Fs / (2 * PI) / 1e3, dB(np.abs(f.reais_em(wg)))[pk], "o",
               ms=6, mfc=COR_HL, mec="k", zorder=5,
               label="pico do filtro digital")
        a.set_title(ttl, fontsize=9.5)
        a.set_xlabel("frequência  [kHz]")
        a.set_ylim(-130, -40)
        a.set_xlim(0, 1.6 * Fs / 1e3)
    ax[0].set_ylabel("magnitude  [dB] (valor absoluto)")
    ax[0].legend(loc="upper right", fontsize=7)
    ax[1].annotate("", xy=(6.0, -50), xytext=(9.0, -50),
                   arrowprops=dict(arrowstyle="->", lw=1.6, color=COR_DIGITAL))
    ax[1].text(7.5, -47, "9 kHz $\\to$ 6 kHz", fontsize=8.5, color=COR_DIGITAL,
               ha="center")
    ax[2].annotate("", xy=(2.0, -50), xytext=(6.0, -50),
                   arrowprops=dict(arrowstyle="->", lw=1.6, color=COR_DIGITAL))
    ax[2].text(4.0, -47, "6 kHz $\\to$ 2 kHz", fontsize=8.5, color=COR_DIGITAL,
               ha="center")
    fig.suptitle("Mecanismo do aliasing: a TDFT do filtro digital e a soma das "
                 "réplicas deslocadas de $kF_s$ (+ termo de meia-descontinuidade)",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    salva(fig, "fig04_aliasing_replicas")


# ===========================================================================
# FIGURA 5 - mapa de polos e zeros
# ===========================================================================
def fig05_polos_zeros(f40, p) -> None:
    f15 = core.discretiza_invariancia_impulso(15_000.0)
    fig, ax = plt.subplots(1, 3, figsize=(12.4, 4.3))

    ang = np.linspace(0, 2 * PI, 400)
    # (a) plano s
    ax[0].axhline(0, color="k", lw=0.8)
    ax[0].axvline(0, color="k", lw=0.8)
    ax[0].plot([0], [p.wd], "x", ms=11, mew=2.2, color=COR_ANALOG,
               label=r"polos $s=-\sigma\pm j\omega_d$")
    ax[0].plot([0], [-p.wd], "x", ms=11, mew=2.2, color=COR_ANALOG)
    ax[0].plot([0], [0], "o", ms=9, mfc="w", mec=COR_DIGITAL, mew=2,
               label="zero $s=0$ (DC)")
    ax[0].set_xlabel(r"$\sigma$  [rad/s]")
    ax[0].set_ylabel(r"$j\Omega$  [rad/s]")
    ax[0].set_title("(a) Plano $s$: protótipo analógico")
    ax[0].set_xlim(-900, 200)
    ax[0].set_ylim(-45e3, 45e3)
    ax[0].legend(loc="lower right", fontsize=8)

    # (b) plano z, Fs = 40 kHz
    ax[1].plot(np.cos(ang), np.sin(ang), color="k", lw=1.0)
    pol = f40.polos
    zer = f40.zeros
    ax[1].plot(pol.real, pol.imag, "x", ms=12, mew=2.2, color=COR_ANALOG,
               label=r"polos: $e^{sT_s}$, $|z|=r=%.6f$" % f40.r)
    ax[1].plot(zer.real, zer.imag, "o", ms=9, mfc="w", mec=COR_DIGITAL, mew=2,
               label="zero em $z=%.4f$" % zer[0].real)
    ax[1].plot([1], [0], "s", ms=6, color=COR_HL, label="$z=1$ (DC, $\\omega=0$)")
    ax[1].set_xlabel(r"$\mathrm{Re}\{z\}$")
    ax[1].set_ylabel(r"$\mathrm{Im}\{z\}$")
    ax[1].set_title(r"(b) Plano $z$: $F_s=40$ kHz  ($\theta=0{,}29998\pi$)")
    ax[1].set_xlim(-1.15, 1.15)
    ax[1].set_ylim(-1.15, 1.15)
    ax[1].set_aspect("equal")
    ax[1].legend(loc="lower left", fontsize=7.5)

    # (c) comparacao dos dois filtros (disco unitario completo)
    ax[2].plot(np.cos(ang), np.sin(ang), color="k", lw=1.0)
    pol15 = f15.polos
    ax[2].plot(np.real(pol), np.imag(pol), "x", ms=12, mew=2.4, color=COR_ANALOG,
               label=r"polos $F_s=40$ kHz ($\theta=0{,}29998\pi$, $r=%.6f$)" % f40.r)
    ax[2].plot(np.real(pol15), np.imag(pol15), "x", ms=12, mew=2.4,
               color=COR_EXTRA,
               label=r"polos $F_s=15$ kHz ($\theta=0{,}79988\pi$, $r=%.6f$)" % f15.r)
    ax[2].plot(np.real(zer), np.imag(zer), "o", ms=9, mfc="w", mec=COR_DIGITAL,
               mew=2, label="zeros $z=%.4f$ ($F_s=40$ kHz)" % zer[0].real)
    ax[2].plot(np.real(f15.zeros), np.imag(f15.zeros), "o", ms=9, mfc="w",
               mec=COR_HL, mew=2, label="zero $z=%.4f$ ($F_s=15$ kHz)"
               % f15.zeros[0].real)
    ax[2].plot([1], [0], "s", ms=6, color=COR_HL, label="$z=1$: eixo $\\omega=0$ / DC")
    for pz in (pol[0], pol15[0]):
        ax[2].plot([0, pz.real], [0, pz.imag], color="k", lw=0.6, ls=":", alpha=0.7)
    ax[2].annotate(r"$\theta=0{,}3\pi$", xy=(pol.real[0], pol.imag[0]),
                   xytext=(0.10, 0.86), fontsize=8,
                   arrowprops=dict(arrowstyle="->", lw=0.8))
    ax[2].annotate(r"$\theta=0{,}8\pi$", xy=(pol15.real[0], pol15.imag[0]),
                   xytext=(-0.95, 0.60), fontsize=8,
                   arrowprops=dict(arrowstyle="->", lw=0.8))
    ax[2].set_xlim(-1.08, 1.08)
    ax[2].set_ylim(-1.08, 1.08)
    ax[2].set_aspect("equal")
    ax[2].set_title("(c) Comparacao direta: os 4 polos e os zeros com $|z|<1$")
    ax[2].set_xlabel(r"$\mathrm{Re}\{z\}$")
    ax[2].set_ylabel(r"$\mathrm{Im}\{z\}$")
    ax[2].legend(loc="lower left", fontsize=6.8)

    fig.suptitle("Desafio 08 - Mapa de polos e zeros: polos analógicos mapeiam por "
                 r"$z=e^{sT_s}$ (estabilidade preservada, $|z|<1$)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    salva(fig, "fig05_polos_zeros")


# ===========================================================================
# FIGURA 6 - falha no passa-altas
# ===========================================================================
def fig06_passa_altas(dados) -> None:
    fhp = dados["fhp"]
    a = core.HP_A_CONST
    fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.3))

    ff = np.linspace(1.0, 21e3, 40001)
    ax[0].plot(ff / 1e3, dB(core.resposta_analogica_passa_altas(2 * PI * ff)),
               color=COR_ANALOG, lw=2.4, label="analógico $|H_a(j2\\pi f)|$ (passa-altas)")
    w = np.linspace(0, PI, 12001)
    ax[0].plot(w * fhp.Fs / (2 * PI) / 1e3, dB(np.abs(fhp.reais_em(w))),
               color=COR_DIGITAL, lw=1.7, ls="--",
               label="digital $|H(e^{j\\omega})|$ (invariância ao impulso)")
    # variante ignorando o termo delta[n] (implementacao pratica do metodo)
    w_sem = w
    H_sem = -a * fhp.Ts / (1 - dados["alpha"] * np.exp(-1j * w_sem))
    ax[0].plot(w_sem * fhp.Fs / (2 * PI) / 1e3, dB(np.abs(H_sem)),
               color=COR_EXTRA, lw=1.1, ls=":",
               label="digital sem o termo $T_s\\delta[n]$ (coincide)")
    ax[0].axvline(1.0, color=COR_HL, ls=":", lw=1.3)
    ax[0].text(1.2, -25, "corte analógico\n1 kHz", fontsize=8, color=COR_HL)
    ax[0].set_xlabel("frequência  [kHz]")
    ax[0].set_ylabel("magnitude  [dB]")
    ax[0].set_xlim(0, 21)
    ax[0].set_ylim(-30, 3)
    ax[0].set_title("(a) O passa-altas analógico vira PASSA-BAIXAS digital "
                    "($F_s = 40$ kHz)")
    ax[0].legend(loc="upper right", fontsize=7.5)
    ax[0].annotate("DC: %.3f (%.2f dB)" % (dados["gDC"], 20 * np.log10(dados["gDC"])),
                   xy=(0.05, 20 * np.log10(dados["gDC"])), xytext=(3.0, 1.6),
                   arrowprops=dict(arrowstyle="->", lw=0.9), fontsize=8.5)
    ax[0].annotate("Nyquist (20 kHz): %.3f (%.2f dB)"
                   % (dados["gNy"], 20 * np.log10(dados["gNy"])),
                   xy=(20.0, 20 * np.log10(dados["gNy"])), xytext=(8.0, -21),
                   arrowprops=dict(arrowstyle="->", lw=0.9), fontsize=8.5)

    # (b) replicas nao decaem -> somatorio divergente
    base = np.linspace(0.0, 1.5e5, 40001)
    for k in range(-1, 3):
        ax[1].plot(
            base / 1e3,
            dB(core.resposta_analogica_passa_altas(2 * PI * np.abs(base - k * 40e3))),
            color=COR_REPLICA, lw=1.9 if k == 0 else 1.1,
            ls="-" if k == 0 else "--",
            label=("resposta analógica $|H_a(j\\Omega)|$" if k == 0 else
                   r"réplicas deslocadas de $kF_s$" if k == 1 else None),
        )
    ax[1].axvspan(0, 20, color=COR_EXTRA, alpha=0.10,
                  label=r"banda-base digital $0\ldots F_s/2$")
    ax[1].set_xlabel("frequência  [kHz]")
    ax[1].set_ylabel("magnitude  [dB]")
    ax[1].set_ylim(-30, 3)
    ax[1].set_title("(b) Replicas de ganho $\\approx 1$ em TODO $k$: "
                    "soma das réplicas não converge")
    ax[1].legend(loc="lower right", fontsize=8)
    ax[1].text(21, -3.0, "o metodo so funciona se $|H_a(j\\Omega)|\\to 0$\n"
               "quando $\\Omega\\to\\infty$ (resposta limitada em banda)",
               fontsize=8.5, color=COR_DIGITAL,
               bbox=dict(fc="white", ec=COR_DIGITAL, alpha=0.9))

    fig.suptitle("Desafio de Reflexão da Dupla - Por que a invariância ao impulso "
                 "falha miseravelmente no filtro passa-altas", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    salva(fig, "fig06_falha_passa_altas")


# ===========================================================================
# FIGURA 7 - erro de aliasing em funcao de Fs
# ===========================================================================
def fig07_erro_vs_Fs(dados) -> None:
    pa = dados["analógico"]
    Fss = np.linspace(12_050.0, 60_000.0, 120)
    erros = []
    for Fs in Fss:
        f = core.discretiza_invariancia_impulso(float(Fs))
        w = np.linspace(0, PI, 30001)
        m = np.abs(f.reais_em(w))
        erros.append(20 * np.log10(m.max() / pa["ganho_em_wn"]))
    erros = np.array(erros)

    fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.3))

    ax[0].plot(Fss / 1e3, erros, color=COR_DIGITAL, lw=1.9,
               label="erro do pico do filtro digital (vs. analógico)")
    ax[0].axhline(0, color="k", lw=0.9)
    for Fs, cor, txt, dxy in ((12.0, COR_HL, "+6,1 dB:\nressonância em Nyquist\n"
                               "(imagem em fase)", (2.0, -0.5)),
                              (15.0, COR_EXTRA, "+0,19 dB", (2.6, 0.30)),
                              (40.0, COR_ANALOG, "+0,07 dB (nominal)", (3.0, 0.50))):
        f = core.discretiza_invariancia_impulso(Fs * 1e3)
        w = np.linspace(0, PI, 100001)
        m = np.abs(f.reais_em(w))
        e = 20 * np.log10(m.max() / pa["ganho_em_wn"])
        ax[0].plot([Fs], [e], "o", ms=7, mfc=cor, mec="k", zorder=5)
        ax[0].annotate(txt, xy=(Fs, e), xytext=(Fs + dxy[0], e + dxy[1]),
                       fontsize=8, color="k",
                       arrowprops=dict(arrowstyle="->", lw=0.8))
    ax[0].set_xlabel(r"$F_s$  [kHz]")
    ax[0].set_ylabel("erro do pico em relação ao analógico  [dB]")
    ax[0].set_title(r"(a) Erro de aliasing no pico ($F_s \geq 12$ kHz)")
    ax[0].set_xlim(12, 60)
    ax[0].set_xticks([12, 15, 20, 25, 30, 40, 50, 60])
    ax[0].set_ylim(-0.6, 5.4)
    ax[0].legend(loc="upper right", fontsize=8)

    Fss2 = np.linspace(4_100.0, 40_000.0, 400)
    f_pk2 = []
    for Fs in Fss2:
        f = core.discretiza_invariancia_impulso(float(Fs))
        w = np.linspace(0, PI, 20001)
        m = np.abs(f.reais_em(w))
        f_pk2.append(w[int(np.argmax(m))] * Fs / (2 * PI) / 1e3)
    ax[1].plot(Fss2 / 1e3, f_pk2, color=COR_DIGITAL, lw=2.0, zorder=4,
               label="frequência do pico do filtro digital")
    ax[1].axhline(6.0, color=COR_ANALOG, lw=2.2, ls="--", alpha=0.75, zorder=3,
                  label="máscara: ressonância em 6 kHz")
    ax[1].axvline(12.0, color="k", ls=":", lw=1.2)
    ax[1].axvspan(4.1, 12.0, color=COR_HL, alpha=0.10)
    ax[1].text(12.4, 7.2, r"$F_s = 12$ kHz: $F_s/2 = f_0$", fontsize=8.5)
    ax[1].text(5.0, 4.4, "ressonância\nsubamostrada\n(dobra para $F_s-f_0$)",
               fontsize=8, color=COR_HL)
    ax[1].set_xlabel(r"$F_s$  [kHz]")
    ax[1].set_ylabel("frequência de pico  [kHz]")
    ax[1].set_title(r"(b) Abaixo de $F_s=12$ kHz a ressonância e subamostrada")
    ax[1].legend(loc="center right", fontsize=8)
    ax[1].set_ylim(0, 8)
    ax[1].set_xlim(4, 40)

    fig.suptitle("Desafio 08 - Item 3: quando o aliasing passa a ser significativo "
                 "(varredura da taxa de amostragem)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    salva(fig, "fig07_erro_aliasing_vs_Fs")


# ===========================================================================
# RELATORIO
# ===========================================================================
def escreve_relatorio(dados3, dados_ref) -> None:
    f40 = core.discretiza_invariancia_impulso()
    pa, pd, banda = dados3["analógico"], dados3["digital"], dados3["banda"]

    linhas = [
        "# Desafio 08 - PDS 2026.2 - Resultados numéricos",
        "",
        "Gerado automaticamente por `src/desafio08.py`.",
        "",
        "## Parâmetros do protótipo analógico",
        "",
        "| Grandeza | Valor |",
        "|---|---|",
        f"| $\\omega_n$ | 12000π = {core.WN:.6f} rad/s (f_n = 6000 Hz) |",
        f"| $\\zeta$ | 1/120 = {core.parametros_analogicos().zeta:.10f} |",
        f"| $Q$ | 60 |",
        f"| $\\sigma=\\zeta\\omega_n$ | 100π = {core.parametros_analogicos().sigma:.6f} rad/s |",
        f"| $\\omega_d$ | 100π√14399 = {core.parametros_analogicos().wd:.6f} rad/s |",
        f"| $T_s$ | 25 µs (F_s = 40 kHz) |",
        "",
        "## Item 1 - ha(t)",
        "",
        "$$h_a(t)=e^{-100\\pi t}\\left[\\cos(\\omega_d t)"
        "-\\frac{1}{\\sqrt{14399}}\\sin(\\omega_d t)\\right]u(t),\\quad"
        "\\omega_d=100\\pi\\sqrt{14399}$$",
        "",
        f"Forma exata (SymPy): $(-\\sqrt{{119}}\\sin(1100\\sqrt{{119}}\\pi t)/1309"
        f"+\\cos(1100\\sqrt{{119}}\\pi t))e^{{-100\\pi t}}$",
        "",
        "## Item 2 - Coeficientes de H(z)",
        "",
        f"$r=e^{{-\\sigma T_s}}={f40.r:.12f}$, "
        f"$\\theta=\\omega_dT_s={f40.theta:.12f}$ rad "
        f"(${f40.theta / PI:.10f}\\pi$ rad/amostra)",
        "",
        "$$H(z)=\\frac{b_0+b_1z^{-1}}{1+a_1z^{-1}+a_2z^{-2}}$$",
        "",
        "| Coeficiente | Expressão | Valor |",
        "|---|---|---|",
        f"| $b_0$ | $T_s$ | {f40.b[0]:.12e} |",
        f"| $b_1$ | $-T_s r(\\cos\\theta+\\frac{{\\sigma}}{{\\omega_d}}\\sin\\theta)$ | {f40.b[1]:.12e} |",
        f"| $a_1$ | $-2r\\cos\\theta$ | {f40.a[1]:.12e} |",
        f"| $a_2$ | $r^2$ | {f40.a[2]:.12e} |",
        "",
        f"Polos: $z={f40.r:.9f}e^{{\\pm j{f40.theta:.9f}}}$ (|z| < 1 : estável). "
        f"Zero: $z={f40.zeros[0].real:.9f}$.",
        "",
        "## Item 3 - Ganho de pico",
        "",
        "| Filtro | Ganho de pico | dB |",
        "|---|---|---|",
        f"| Analógico (em $\\Omega=\\omega_n$) | {pa['ganho_em_wn']:.10e} | {20 * np.log10(pa['ganho_em_wn']):.4f} |",
        f"| Analógico (máximo exato) | {pa['ganho']:.10e} | {pa['ganho_dB']:.4f} |",
        f"| Digital em $\\omega_0=0{{,}}3\\pi$ | {pd['ganho_em_w0']:.10e} | {pd['ganho_em_w0_dB']:.4f} |",
        "",
        f"Diferença em relação ao pico analógico ({pa['ganho_dB']:.4f} dB): "
        f"**{pd['ganho_em_w0_dB'] - pa['ganho_dB']:+.4f} dB** "
        f"({100 * (pd['ganho_em_w0'] / pa['ganho'] - 1):+.4f} %) ; em relação a "
        f"$|H_a(j\\omega_n)|$ ({20 * np.log10(pa['ganho_em_wn']):.4f} dB): "
        f"{pd['ganho_em_w0_dB'] - 20 * np.log10(pa['ganho_em_wn']):+.4f} dB "
        "-> aliasing desprezivel em F_s = 40 kHz.",
        "",
        f"Banda de -3 dB do filtro digital: {banda['df_3dB_Hz']:.4f} Hz.",
        "",
        "### Varredura da taxa de amostragem (aliasing)",
        "",
        "| $F_s$ [kHz] | pico [dB] | freq. do pico [Hz] | razão vs. analógico |",
        "|---|---|---|---|",
    ]
    for Fs in (40_000.0, 20_000.0, 15_000.0, 12_000.0, 8_000.0):
        fx = core.discretiza_invariancia_impulso(Fs)
        px = core.pico_digital(fx, n_pontos=200_001)
        linhas.append(f"| {Fs / 1000:.1f} | {px['ganho_pico_dB']:.3f} | "
                      f"{px['f_pico']:.2f} | {px['ganho_pico'] / pa['ganho_em_wn']:.3f} |")
    fhp = dados_ref["fhp"]
    linhas += [
        "",
        "## Desafio de Reflexão - passa-altas",
        "",
        f"$H_{{hp}}(z)=\\frac{{({fhp.b[0]:.9e}){fhp.b[1]:+.9e}z^{{-1}}}}"
        f"{{1{fhp.a[1]:+.9f}z^{{-1}}}}$",
        "",
        "| Frequência | Ganho | dB |",
        "|---|---|---|",
        f"| DC (z=1) | {dados_ref['gDC']:.6f} | {20 * np.log10(dados_ref['gDC']):+.4f} |",
        f"| Nyquist (z=-1) | {dados_ref['gNy']:.6f} | {20 * np.log10(dados_ref['gNy']):+.4f} |",
        "",
        "O ganho decresce com a frequência: o passa-altas degenerou em passa-baixas.",
        "Causa: $|H_a(j\\Omega)|\\to1$ quando $\\Omega\\to\\infty$, logo as réplicas da",
        "identidade de Poisson não decaem e a soma diverge; o zero em $s=0$ não e",
        "transportado para $z=1$ (zeros não mapeiam por $z=e^{sT_s}$).",
        "",
        "## Observação fina (comprovada em verificacao.py)",
        "",
        "A identidade de Poisson exata para este problema exige o termo de",
        "meia-descontinuidade, pois $h_a(0^+)=1$ e $h_a(0^-)=0$:",
        "",
        "$$H(e^{j\\omega})=\\sum_k H_a\\!\\left(j\\frac{\\omega-2\\pi k}{T_s}\\right)"
        "+\\frac{T_s\\,h_a(0^+)}{2}$$",
        "",
        f"com $T_s h_a(0^+)/2 = {f40.Ts / 2:.6e}$ (valor conferido numericamente",
        "com erro relativo < 1e-4 na varredura de frequência).",
        "",
    ]
    (RES / "relatorio_resultados.md").write_text("\n".join(linhas), encoding="utf-8")
    log(f"\n  [relatorio] resultados/relatorio_resultados.md")


# ===========================================================================
# MAIN
# ===========================================================================
def main() -> int:
    estilo()
    p = core.parametros_analogicos()
    f40 = core.discretiza_invariancia_impulso(core.FS_NOMINAL)

    log("=" * 78)
    log(" PDS 2026.2 - DESAFIO 08")
    log(" Discretizacao de Filtros IIR por Invariância ao Impulso")
    log(" e Analise de Aliasing Espectral")
    log("=" * 78)

    item1(p)
    item2(f40, p)
    dados3 = item3(f40, p)
    dados_ref = reflexao()

    titulo(5, "PARTE COMPUTACIONAL - geracao das figuras")
    fig01_ha_amostras(f40, p)
    fig02_impulso_digital(f40)
    fig03_magnitude(f40, dados3)
    fig04_replicas()
    fig05_polos_zeros(f40, p)
    fig06_passa_altas(dados_ref)
    fig07_erro_vs_Fs(dados3)

    titulo(6, "RELATORIO")
    escreve_relatorio(dados3, dados_ref)

    log("")
    log("=" * 78)
    log(" FIM - 7 figuras em figuras/ e relatorio em resultados/")
    log("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
