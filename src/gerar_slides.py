#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Autor: Matheus Sousa, LCSI-DES, UFPE.
"""
===============================================================================
 PDS 2026.2 - DESAFIO 08 - GERADOR DOS SLIDES (PDF)
===============================================================================

 Monta o arquivo  slides_desafio08.pdf  (formato 16:9) com o texto da
 apresentacao oral e as figuras geradas por src/desafio08.py.

 Nao depende de nenhuma biblioteca extra: usa matplotlib.backends.backend_pdf,
 portanto o mesmo ambiente que roda a simulacao gera os slides.

 Uso:
     python src/desafio08.py        # gera as figuras (figuras/*.png)
     python src/gerar_slides.py     # gera slides_desafio08.pdf
===============================================================================
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

RAIZ = Path(__file__).resolve().parent.parent
FIG = RAIZ / "figuras"
SAIDA = RAIZ / "slides_desafio08.pdf"

AZUL = "#1f4e79"
VERM = "#c00000"
CINZA = "#404040"
LAR, ALT = 13.333, 7.5          # polegadas (16:9)


def slide_base(titulo: str, subtitulo: str | None = None):
    fig = plt.figure(figsize=(LAR, ALT))
    fig.patch.set_facecolor("white")
    fig.text(0.045, 0.935, titulo, fontsize=23, weight="bold", color=AZUL,
             va="top")
    fig.patches.append(plt.Rectangle((0.045, 0.892), 0.91, 0.004,
                                     transform=fig.transFigure,
                                     facecolor=VERM))
    if subtitulo:
        fig.text(0.045, 0.868, subtitulo, fontsize=13, color=CINZA, va="top")
    return fig


def slide_bullets(titulo, itens, subtitulo=None, rodape=None):
    fig = slide_base(titulo, subtitulo)
    y = 0.815
    for it in itens:
        if it.startswith("##"):
            y -= 0.030
            fig.text(0.055, y, it[2:].strip(), fontsize=15, weight="bold",
                     color=VERM, va="top")
            y -= 0.075
            continue
        if it.startswith("$"):       # linha de formula (mathtext)
            fig.text(0.095, y, it, fontsize=16, color="black", va="top")
            y -= 0.085
            continue
        fig.text(0.062, y, "\u2022", fontsize=14, color=AZUL, va="top")
        fig.text(0.088, y, it, fontsize=14, color="black", va="top",
                 wrap=True)
        y -= 0.070 + 0.045 * (len(it) // 95)
    if rodape:
        fig.text(0.045, 0.035, rodape, fontsize=10.5, color=CINZA, va="bottom")
    return fig


def slide_figura(titulo, imagem: str, texto: str | None = None,
                 subtitulo=None, rodape=None):
    fig = slide_base(titulo, subtitulo)
    ax = fig.add_axes([0.045, 0.10, 0.91, 0.74])
    ax.imshow(mpimg.imread(FIG / imagem))
    ax.axis("off")
    if texto:
        fig.text(0.045, 0.070, texto, fontsize=12, color=CINZA, va="bottom")
    if rodape:
        fig.text(0.955, 0.035, rodape, fontsize=10.5, color=CINZA, va="bottom",
                 ha="right")
    return fig


def main() -> int:
    n_slides = 0
    with PdfPages(SAIDA) as pdf:
        # ------------------------------------------------------------ 1
        f = plt.figure(figsize=(LAR, ALT))
        f.patch.set_facecolor("white")
        f.patches.append(plt.Rectangle((0, 0), 1, 0.16,
                                       transform=f.transFigure,
                                       facecolor=AZUL))
        f.text(0.5, 0.76, "Desafio 08", fontsize=40, weight="bold",
               color=AZUL, ha="center")
        f.text(0.5, 0.63, "Discretização de Filtros IIR por Invariância ao "
               "Impulso\ne Análise de Aliasing Espectral", fontsize=22,
               color=CINZA, ha="center")
        f.text(0.5, 0.44, "PDS 2026.2 — Departamento de Eletrônica e Sistemas "
               "(UFPE)", fontsize=15, ha="center", color="black")
        f.text(0.5, 0.38, "Avaliação em duplas — 10 minutos "
               "(7–8 min de exposição + 2–3 min de arguição)",
               fontsize=13, ha="center", color=VERM)
        f.text(0.5, 0.27, "Dupla: ______________________  |  "
               "______________________", fontsize=14, ha="center",
               color=CINZA)
        f.text(0.5, 0.13, "F_s = 40 kHz  (T_s = 25 µs)  •  f_0 = 6 kHz  •  "
               "Q = 60", fontsize=13, ha="center", color=AZUL)
        pdf.savefig(f)
        plt.close(f)

        # ------------------------------------------------------------ 2
        f = slide_bullets(
            "Contexto e protótipo analógico",
            ["O enunciado pede um filtro digital passa-faixa obtido a partir "
             "do protótipo analógico:",
             "$H_a(s) = \\frac{s}{s^2 + 200\\pi s + (12000\\pi)^2} = "
             "\\frac{s}{s^2+2\\zeta\\omega_n s+\\omega_n^2}$",
             "## Parâmetros exatos",
             "$\\omega_n = 12000\\pi = 37699{,}1118$ rad/s  ($f_n$ = 6000 Hz)",
             "$2\\zeta\\omega_n = 200\\pi \\Rightarrow \\zeta = 1/120 = "
             "8{,}333\\times10^{-3} \\Rightarrow Q = 60$",
             "$\\sigma = \\zeta\\omega_n = 100\\pi$, $\\omega_d = 100\\pi"
             "\\sqrt{14399} = 37697{,}8028$ rad/s",
             "Método exigido: INVARIÂNCIA AO IMPULSO com $T_s = 25\\ \\mu$s "
             "($F_s$ = 40 kHz)",
             "Definição:  $h[n] = T_s\\,h_a(nT_s)$, ou seja, a sequência é a "
             "resposta ao impulso analógica amostrada e escalada por $T_s$."],
            rodape="Desafio 08 — PDS 2026.2")
        pdf.savefig(f)
        plt.close(f)

        # ------------------------------------------------------------ 3
        f = slide_bullets(
            "Item 1 — resposta ao impulso analógica $h_a(t)$",
            ["Frações parciais:",
             "$\\frac{s}{(s+\\sigma)^2+\\omega_d^2} = "
             "\\frac{s+\\sigma}{(s+\\sigma)^2+\\omega_d^2} - "
             "\\frac{\\sigma}{(s+\\sigma)^2+\\omega_d^2}$",
             "## Resultado (tabelas de Laplace)",
             "$h_a(t) = e^{-100\\pi t}\\left[\\cos(\\omega_d t) - "
             "\\frac{1}{\\sqrt{14399}}\\sin(\\omega_d t)\\right]u(t)$",
             "$\\omega_d = 100\\pi\\sqrt{14399} = 37697{,}8028$ rad/s  "
             "($f_d$ = 5999,79 Hz)",
             "## Pontos de atenção",
             "$h_a(0^+) = 1$ e $h_a(0^-) = 0$: existe um SALTO em $t=0$ "
             "(será usado na análise de aliasing)",
             "$|H_a(j\\omega_n)| = \\frac{1}{2\\zeta\\omega_n} = "
             "\\frac{1}{200\\pi} = 1{,}59155\\times10^{-3}$"],
            rodape="Verificado com SymPy (forma exata com $\\sqrt{119}$) e "
                   "numericamente — ver src/verificacao.py")
        pdf.savefig(f)
        plt.close(f)

        # ------------------------------------------------------------ 4
        f = slide_figura(
            "Item 1 — $h_a(t)$ e as amostras $h[n]/T_s = h_a(nT_s)$",
            "fig01_ha_t_e_amostras.png",
            "As amostras caem exatamente sobre a curva analítica: é a "
            "definição do método (não há aproximação no tempo).",
            subtitulo="Oscilação em ~5999,8 Hz com decaimento "
                      "$e^{-\\sigma t}$ (τ = 1/σ = 3,18 ms; Q = 60)")
        pdf.savefig(f)
        plt.close(f)

        # ------------------------------------------------------------ 5
        f = slide_bullets(
            "Item 2 — coeficientes exatos de $H(z)$ (2ª ordem)",
            ["Com $r = e^{-\\sigma T_s}$ e $\\theta = \\omega_d T_s$:",
             "$H(z) = \\frac{T_s\\left[1 - r\\left(\\cos\\theta + "
             "\\frac{\\sigma}{\\omega_d}\\sin\\theta\\right)z^{-1}\\right]}"
             "{1 - 2r\\cos\\theta\\,z^{-1} + r^2z^{-2}}$",
             "$r = e^{-100\\pi\\cdot25\\mu} = 0{,}992176780$   |   "
             "$\\theta = 0{,}942445070$ rad = $0{,}29998958\\pi$ rad/amostra",
             "## Coeficientes numéricos",
             "$b_0 = T_s = 2{,}500000000\\times10^{-5}$   |   "
             "$b_1 = -1{,}474755715\\times10^{-5}$",
             "$a_1 = -1{,}166426294$   |   "
             "$a_2 = 0{,}9844147634$",
             "Polos: $z = 0{,}992176780\\,e^{\\pm j0{,}29999\\pi}$ "
             "$\\Rightarrow |z| < 1$ (ESTÁVEL).  Zero do método em "
             "$z = 0{,}58990$.",
             "O ganho em DC vale $|H(e^{j0})| = 1{,}2534\\times10^{-5} "
             "\\approx T_s/2$ — resíduo do termo de meia-descontinuidade, "
             "e não exatamente zero."],
            rodape="Conferido por dois caminhos independentes: resíduos "
                   "complexos dos polos e soma direta (TDFT) de h[n]")
        pdf.savefig(f)
        plt.close(f)

        # ------------------------------------------------------------ 6
        f = slide_figura(
            "Item 2 — polos e zeros no plano z",
            "fig05_polos_zeros.png",
            "Os polos analógicos mapeiam por $z=e^{sT_s}$: raio "
            "$r=e^{-\\sigma T_s}<1$ garante estabilidade. Os ZEROS não são "
            "mapeados — o zero em $s=0$ (DC) não vira zero em $z=1$.",
            subtitulo="Comparação $F_s$ = 40 kHz ($\\theta=0{,}3\\pi$) e "
                      "$F_s$ = 15 kHz ($\\theta=0{,}8\\pi$)")
        pdf.savefig(f)
        plt.close(f)

        # ------------------------------------------------------------ 7
        f = slide_figura(
            "Item 3 — ganho de pico em $\\omega_0 = 2\\pi\\cdot6000\\cdot T_s "
            "= 0{,}3\\pi$",
            "fig03_magnitude_analogico_vs_digital.png",
            "ANALÓGICO: $|H_a(j\\omega_n)| = 1{,}59155\\times10^{-3}$ "
            "(−55,9636 dB)   |   DIGITAL: $|H(e^{j0{,}3\\pi})| = "
            "1{,}60409\\times10^{-3}$ (−55,8954 dB)",
            subtitulo="Diferença de apenas +0,0685 dB (+0,79 %): "
                      "aliasing desprezível com $F_s$ = 40 kHz")
        pdf.savefig(f)
        plt.close(f)

        # ------------------------------------------------------------ 8
        f = slide_figura(
            "Mecanismo do aliasing: réplicas deslocadas de $kF_s$",
            "fig04_aliasing_replicas.png",
            "A TDFT do filtro digital é a soma das réplicas: quanto menor "
            "$F_s$, mais perto da banda-base caem as réplicas.",
            subtitulo="$F_s$ = 40 kHz: réplicas em 34/46 kHz  |  "
                      "$F_s$ = 15 kHz: imagem em 9 kHz dobra para 6 kHz  |  "
                      "$F_s$ = 8 kHz: a ressonância inteira dobra para 2 kHz")
        pdf.savefig(f)
        plt.close(f)

        # ------------------------------------------------------------ 9
        f = slide_figura(
            "Item 3 — quando o aliasing passa a importar",
            "fig07_erro_aliasing_vs_Fs.png",
            "Erro do pico em dB relativo ao protótipo analógico: "
            "$F_s$ = 40 kHz → +0,07 dB (nominal); $F_s$ = 15 kHz → +0,19 dB; "
            "$F_s$ = 12 kHz → +6,1 dB; abaixo de 12 kHz a ressonância de "
            "6 kHz é subamostrada e reaparece em $F_s - 6$ kHz.",
            subtitulo="Resposta ao item 3: NÃO houve alteração significativa "
                      "do pico em relação ao filtro analógico original")
        pdf.savefig(f)
        plt.close(f)

        # ----------------------------------------------------------- 10
        f = slide_figura(
            "Desafio de reflexão — por que o método falha no passa-altas",
            "fig06_falha_passa_altas.png",
            "Passa-altas analógico $H_a(s)=\\frac{s}{s+2000\\pi}$ → "
            "$H(z)$ obtido: ganho em DC = 1,0806 (+0,67 dB) e em Nyquist = "
            "0,0847 (−21,45 dB): o filtro ficou PASSA-BAIXAS.",
            subtitulo="Causa: $|H_a(j\\Omega)| \\to 1$ quando "
                      "$\\Omega\\to\\infty$ — a soma das réplicas não converge")
        pdf.savefig(f)
        plt.close(f)

        # ----------------------------------------------------------- 11
        f = slide_bullets(
            "Desafio de reflexão — a explicação matemática",
            ["1) A identidade exata do método é uma soma periódica "
             "(Poisson):",
             "$H(e^{j\\omega}) = \\sum_{k} H_a\\!\\left(j\\frac{\\omega-2\\pi "
             "k}{T_s}\\right) + \\frac{T_s h_a(0^+)}{2}$",
             "2) Para o passa-altas cada réplica vale ~1 em TODO $k$: "
             "a série NÃO converge (energia de todo o eixo Ω entra na "
             "banda digital).",
             "3) Nenhum $F_s$ resolve: o problema é do MÉTODO, não da taxa.",
             "4) O zero em $s=0$ que definia o passa-altas NÃO é transportado "
             "para $z=1$ — na invariancia ao impulso só os POLOS mapeiam, por "
             "$z=e^{sT_s}$.",
             "5) Consequência: o termo $T_s\\delta[n]$ (de $\\delta(t)$) "
             "somado ao passa-baixas do polo produz a curva decrescente "
             "observada.",
             "Para passa-altas usa-se TRANSFORMAÇÃO BILINEAR (com "
             "pré-deformação) — é o método do Desafio 09."],
            rodape="Conclusão: a invariancia ao impulso só é válida para "
                   "respostas limitadas em banda ($|H_a| \\to 0$ em HF)")
        pdf.savefig(f)
        plt.close(f)

        # ----------------------------------------------------------- 12
        f = slide_bullets(
            "Rigor extra e verificação numérica",
            ["Termo de meia-descontinuidade: como $h_a(0^+)=1$, a soma "
             "clássica de réplicas fica exatamente $T_s/2 = 12{,}5\\ \\mu$s "
             "abaixo da resposta exata.",
             "Comprovado na simulação: desvio medido "
             "$1{,}2500\\times10^{-5}$ (previsto $1{,}2500\\times10^{-5}$); "
             "erro relativo final $6{,}5\\times10^{-5}$.",
             "23 verificações independentes (src/verificacao.py), incluindo:",
             "## Três caminhos para $H(z)$",
             "$\\bullet$ forma fechada (tabelas)  $\\bullet$ resíduos "
             "complexos dos polos  $\\bullet$ TDFT/FFT de $h[n]$",
             "## Mapa de polos e escala de $F_s$",
             "$\\bullet$ polos digitais $= e^{p_aT_s}$ (exato)  $\\bullet$ "
             "$F_s$=12 kHz: pico dobra (2,026×) como previsto pelo "
             "dobramento em Nyquist"],
            rodape="Tudo reproduzível: python src/verificacao.py → "
                   "“TODAS AS VERIFICAÇÕES PASSARAM” (23/23)")
        pdf.savefig(f)
        plt.close(f)

        # ----------------------------------------------------------- 13
        f = slide_figura(
            "Resposta ao impulso do filtro digital",
            "fig02_resposta_impulso_digital.png",
            "$h[n]$ decai lentamente porque os polos estão muito próximos do "
            "círculo unitário ($r = 0{,}9922$): o filtro tem banda estreita "
            "(≈100 Hz de banda de −3 dB) e tempo de acomodação longo.",
            subtitulo="Estabilidade garantida por $|z| = r = e^{-\\sigma T_s} "
                      "< 1$")
        pdf.savefig(f)
        plt.close(f)

        # ----------------------------------------------------------- 14
        f = slide_bullets(
            "Conclusões",
            ["1) O protótipo tem $f_n$ = 6 kHz e Q = 60; sua resposta ao "
             "impulso e sua amostragem por $T_s$ = 25 µs foram obtidas de "
             "forma fechada e verificadas por três caminhos independentes.",
             "2) A invariancia ao impulso gerou um IIR de 2ª ordem exato "
             "($b_0,b_1,a_1,a_2$ acima), estável, com polos em "
             "$0{,}99218\\,e^{\\pm j0{,}3\\pi}$.",
             "3) O ganho de pico em $\\omega_0 = 0{,}3\\pi$ difere do "
             "analógico em +0,07 dB: NÃO houve alteração significativa, "
             "porque $F_s/2 = 20$ kHz ≫ 6 kHz (réplicas 50 dB abaixo do pico).",
             "4) O aliasing cresce ao reduzir $F_s$: +0,19 dB (15 kHz), "
             "+6,1 dB (12 kHz) e subamostragem real abaixo de 12 kHz — "
             "a fronteira é $F_s = 2f_0$.",
             "5) No passa-altas o método DEGENERA (resposta decrescente) "
             "porque as réplicas não decaem; a solução é a transformação "
             "bilinear.",
             "Próximos passos naturais: pré-deformação e bilinear (Desafio 09) "
             "e banda estreita em cascata (SOS, Desafio 06)."],
            rodape="Obrigado! Perguntas?")
        pdf.savefig(f)
        plt.close(f)

    print(f"[slides] {SAIDA} ({SAIDA.stat().st_size / 1024:.0f} kB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
