#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Autor: Matheus Sousa, LCSI-DES, UFPE.
"""
===============================================================================
 PDS 2026.2 - DESAFIO 08 - VERIFICACAO INDEPENDENTE DOS RESULTADOS
===============================================================================

 Cada bloco abaixo confere um resultado da resolucao POR UM CAMINHO DIFERENTE
 do usado no script principal. Se todas as checagens passarem, o script
 termina com "TODAS AS VERIFICACOES PASSARAM" (exit code 0).

 Uso:  python src/verificacao.py
===============================================================================
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pds08_core as core   # noqa: E402

TOL = 1e-6
LINHA = "-" * 78


def checa(nome: str, valor: float, esperado: float, tol: float = TOL) -> bool:
    ok = abs(valor - esperado) <= tol * max(1.0, abs(esperado))
    print(f"  [{'OK ' if ok else 'FALHA'}] {nome:<58s} "
          f"{valor:.10g}  (ref {esperado:.10g})")
    return ok


def main() -> int:
    resultados: list[bool] = []

    # ------------------------------------------------------------------
    print(LINHA)
    print("1) ha(t) simbolica (SymPy) x ha(t) fechada (core)")
    print(LINHA)
    s, t = sp.symbols("s t", positive=True)
    wn = 12000 * sp.pi
    Ha = s / (s ** 2 + 200 * sp.pi * s + wn ** 2)
    ha_simbolica = sp.simplify(sp.inverse_laplace_transform(Ha, s, t))
    print("  ha(t) obtida pelo SymPy:")
    print("   ", ha_simbolica)
    # comparacao numerica ponto a ponto com a forma fechada usada no codigo
    ha_sim_num = sp.lambdify(t, ha_simbolica, "numpy")
    tt = np.array([0.0, 1e-5, 2.5e-5, 1e-4, 5e-4, 1e-3])
    ref = np.array([float(ha_sim_num(x)) for x in tt])
    calc = core.ha_analogica(tt)
    err = np.max(np.abs(ref - calc) / np.max(np.abs(ref)))
    resultados.append(checa("erro relativo maximo ha(t) simbolico x fechado", err, 0.0, 1e-9))

    # ------------------------------------------------------------------
    print(LINHA)
    print("2) Coeficientes de H(z): formula fechada x ajuste exato (algebra)")
    print(LINHA)
    f = core.discretiza_invariancia_impulso()
    # Residuo: H(z) = Ts*[p1/(2j wd)/(1-exp(p1*Ts) z^-1) + conjugado]
    p = core.parametros_analogicos()
    p1 = -p.sigma + 1j * p.wd
    Ts = f.Ts
    # racional obtido somando os dois polos complexos conjugados
    # montagem exata: A/(1-p1d z^-1) + A*/(1-p2d z^-1)
    p1d, p2d = np.exp(p1 * Ts), np.exp(np.conj(p1) * Ts)
    A = p1 / (p1 - np.conj(p1))
    # (A(1-p2d z^-1) + A*(1-p1d z^-1)) / ((1-p1d z^-1)(1-p2d z^-1))
    num_res = Ts * np.array([A + np.conj(A), -(A * p2d + np.conj(A) * p1d)])
    den_res = np.poly([p1d, p2d])                     # em ordem decrescente
    den_res = np.real(den_res) / den_res[0]           # normaliza a0 = 1
    print("  residuos   b =", np.round(np.real(num_res), 12),
          " a =", np.round(den_res, 12))
    print("  fechada    b =", np.round(f.b, 12), " a =", np.round(f.a, 12))
    resultados.append(checa("b0 residuos x fechado", np.real(num_res[0]), f.b[0]))
    resultados.append(checa("b1 residuos x fechado", np.real(num_res[1]), f.b[1]))
    resultados.append(checa("a1 residuos x fechado", den_res[1], f.a[1]))
    resultados.append(checa("a2 residuos x fechado", den_res[2], f.a[2]))

    # ------------------------------------------------------------------
    print(LINHA)
    print("3) H(e^jw) racional x TDFT (soma direta) x FFT de h[n]")
    print(LINHA)
    w = np.array([0.0, 0.1, 0.3 * np.pi, 0.5 * np.pi, 0.9 * np.pi, np.pi])
    h_exata = core.resposta_impulso_digital(n_amostras=400_000)
    H_tdft = core.dtft_de_amostras(h_exata, w)
    H_rac = f.reais_em(w)
    erro = np.max(np.abs(H_tdft - H_rac) / np.abs(H_rac))
    print(f"  |H| racional : {np.abs(H_rac)}")
    print(f"  |H| TDFT     : {np.abs(H_tdft)}")
    resultados.append(checa("erro relativo maximo racional x TDFT truncada",
                            erro, 0.0, 1e-6))

    # FFT de N amostras == amostras da DTFT em w = 2*pi*k/N (conferencia)
    N = 2_000_000
    hN = core.resposta_impulso_digital(n_amostras=N)
    H_fft = np.fft.fft(hN)
    k = int(round(0.3 * N / 2))                 # w = 0.3*pi
    w_k = 2 * np.pi * k / N
    resultados.append(checa("|H| FFT(k) x formula em w=0.3*pi",
                            abs(H_fft[k]), abs(f.reais_em(w_k).item()), 1e-9))

    # ------------------------------------------------------------------
    print(LINHA)
    print("4) Soma de Poisson (replicas) x resposta racional")
    print(LINHA)
    # A soma de Poisson converge apenas CONDICIONALMENTE, pois |Ha(j*Omega)|
    # decai como 1/Omega (passa-faixa de 2a ordem). Alem disso ha(t) tem um
    # salto em t = 0 (ha(0+) = 1, ha(0-) = 0), o que exige o termo de
    # meia-descontinuidade Ts*ha(0+)/2 na identidade de Poisson.
    w_ok = np.linspace(0.05 * np.pi, 0.6 * np.pi, 11)
    H_sp = core.soma_de_replicas(f, w_ok, k_max=2000,
                                 correcao_meia_descontinuidade=False)
    H_sc = core.soma_de_replicas(f, w_ok, k_max=2000,
                                 correcao_meia_descontinuidade=True)
    H_r = f.reais_em(w_ok)
    desvio = H_r - H_sp
    print(f"  desvio medio replicas-sem-correcao x exato = "
          f"{desvio.real.mean():.6e} (deve ser Ts*ha(0+)/2 = {f.Ts / 2:.6e})")
    resultados.append(checa("desvio das replicas == Ts*ha(0+)/2",
                            float(np.abs(desvio.real.mean())), f.Ts / 2, 1e-3))
    erro_poisson = np.max(np.abs(H_sc - H_r) / np.abs(H_r))
    print(f"  erro relativo maximo das replicas CORRIGIDAS: {erro_poisson:.3e}")
    resultados.append(checa("erro relativo |H| Poisson corrigido x racional",
                            erro_poisson, 0.0, 2e-3))

    # ------------------------------------------------------------------
    print(LINHA)
    print("5) Ganho de pico: item 3 do desafio x pico varrido x analogico")
    print(LINHA)
    pa = core.pico_analogico()
    pd = core.pico_digital(f)
    print(f"  analogico: pico {pa['ganho']:.10g} em f = {pa['f_r']:.3f} Hz "
          f"({pa['ganho_dB']:.4f} dB)")
    print(f"  analogico em wn: {pa['ganho_em_wn']:.10g}")
    print(f"  digital  : pico {pd['ganho_pico']:.10g} em w = "
          f"{pd['w_pico'] / np.pi:.8f} pi (f = {pd['f_pico']:.4f} Hz)")
    print(f"  digital em w0 = 0.3 pi: {pd['ganho_em_w0']:.10g} "
          f"({pd['ganho_em_w0_dB']:.4f} dB)")
    # valor teorico exato: |Ha(j*wn)| = 1/(2*zeta*wn) = 1/(200*pi)
    teo = 1.0 / core.TWO_ZETA_WN
    resultados.append(checa("|H(e^jw0)| x 1/(2*zeta*wn)", pd['ganho_em_w0'], teo, 1e-3))
    resultados.append(checa("ganho DC do filtro digital (deve ser ~0)",
                            pd['ganho_DC'], 0.0, 1e-2))

    # ------------------------------------------------------------------
    print(LINHA)
    print("6) Aliasing em Fs = 15 kHz (imagem em 9 kHz dobra para 6 kHz)")
    print(LINHA)
    f15 = core.discretiza_invariancia_impulso(15_000.0)
    pd15 = core.pico_digital(f15, n_pontos=600_001)
    print(f"  pico digital Fs=15k: {pd15['ganho_pico']:.10g} em w = "
          f"{pd15['w_pico'] / np.pi:.6f} pi  (f = {pd15['f_pico']:.3f} Hz)")
    # posicao da imagem da componente de -6000 Hz apos amostragem em 15 kHz
    f_img = 15_000.0 - 6000.0
    f_fold = abs(15_000.0 - f_img) if f_img > 7500.0 else f_img
    print(f"  imagem: 15 kHz - 6 kHz = {f_img / 1e3:.3f} kHz -> dobra para "
          f"{f_fold / 1e3:.3f} kHz (== 6 kHz, coincide com o pico)")
    resultados.append(checa("frequencia de dobramento da imagem",
                            f_fold, 6000.0, 1e-12))
    # aumento previsto: 1 + |Ha(j*2*pi*9000)| / |Ha(j*wn)|
    ganho_img = core.ganho_analogico_em(2 * np.pi * f_img).item()
    razao_prev = 1.0 + ganho_img / pa["ganho_em_wn"]
    razao_sim = pd15["ganho_pico"] / pa["ganho_em_wn"]
    print(f"  |Ha(j*2*pi*9kHz)| = {ganho_img:.6e} -> razao prevista "
          f"{razao_prev:.6f} ; razao simulada {razao_sim:.6f}")
    resultados.append(checa("razao pico(Fs=15k)/pico analogico x previsto",
                            razao_sim, razao_prev, 2e-2))
    # caso critico: Fs = 12 kHz -> ressonancia exatamente em Nyquist: pico DOBRA
    f12 = core.discretiza_invariancia_impulso(12_000.0)
    ganho_pi = np.abs(f12.reais_em(np.pi)).item()
    print(f"  Fs = 12 kHz (ressonancia em Nyquist): |H(e^j*pi)| = {ganho_pi:.6e}"
          f"  (2x analogico = {2 * pa['ganho_em_wn']:.6e})")
    resultados.append(checa("pico em Fs=12k == 2x pico analogico",
                            ganho_pi, 2 * pa["ganho_em_wn"], 5e-2))

    # ------------------------------------------------------------------
    print(LINHA)
    print("7) Passa-altas degenerado: H(z) x soma de residuos; ganhos DC/Nyquist")
    print(LINHA)
    fhp = core.filtro_passa_altas_invariancia()
    Ts = fhp.Ts
    a = core.HP_A_CONST
    alpha = np.exp(-a * Ts)
    # h[n] = Ts*delta[n] - (a*Ts)*alpha^n*u[n]  ->  H(z) = [Ts - a*Ts - Ts*alpha*z^-1]/(1 - alpha*z^-1)
    b_ref = np.array([Ts - a * Ts, -Ts * alpha])
    resultados.append(checa("b0 passa-altas", fhp.b[0], b_ref[0]))
    resultados.append(checa("b1 passa-altas", fhp.b[1], b_ref[1]))
    resultados.append(checa("a1 passa-altas", fhp.a[1], -alpha))
    gDC = np.abs(fhp.reais_em(0.0)).item()
    gNy = np.abs(fhp.reais_em(np.pi)).item()
    print(f"  |H(1)| = {gDC:.6f}   |H(-1)| = {gNy:.6f}")
    print(f"  analogico: |Ha(j0)| = {core.resposta_analogica_passa_altas(1e-9)[0]:.3e} "
          f"(zero em DC)  vs  |Ha(j*inf)| -> 1")
    resultados.append(checa("ganho DC > ganho Nyquist (filtro virou passa-baixas)",
                            float(gDC > gNy), 1.0, 1e-12))
    # limite exato: |H(1)| = Ts - a*Ts/(1-alpha)
    gDC_teo = Ts - a * Ts / (1 - alpha)
    resultados.append(checa("|H(1)| x expressao algebrica", gDC, abs(gDC_teo), 1e-9))

    # ------------------------------------------------------------------
    print(LINHA)
    print("8) Propriedades da invariancia ao impulso (conferencias extras)")
    print(LINHA)
    # h[0] = Ts*ha(0) = Ts = b0
    resultados.append(checa("h[0] = Ts*ha(0) = b0", core.ha_analogica(0.0) * Ts, f.b[0]))
    # polos digitais = exp(polos analogicos * Ts)
    polos_d = f.polos
    polos_a = -p.sigma + 1j * p.wd, -p.sigma - 1j * p.wd
    esperados = sorted(np.exp(np.array(polos_a) * Ts), key=lambda z: np.angle(z))
    obtidos = sorted(polos_d, key=lambda z: np.angle(z))
    for i, (zo, ze) in enumerate(zip(obtidos, esperados)):
        resultados.append(checa(f"polo digital {i + 1} = exp(p_a*Ts)", abs(zo), abs(ze), 1e-12))
    # soma dos h[n] deve ser ~0 (Ha(j0) = 0) - Riemann nao exato
    soma = core.resposta_impulso_digital(n_amostras=2_000_000).sum()
    print(f"  sum h[n] = {soma:.3e}  (esperado ~0, pois Ha(j0)=0)")
    resultados.append(checa("sum h[n] proximo de zero", float(abs(soma) < 1e-4), 1.0, 1e-12))

    # ------------------------------------------------------------------
    print(LINHA)
    n_ok = sum(resultados)
    print(f"RESULTADO: {n_ok}/{len(resultados)} verificacoes passaram")
    print(LINHA)
    if n_ok == len(resultados):
        print("TODAS AS VERIFICACOES PASSARAM")
        return 0
    print("EXISTEM VERIFICACOES QUE FALHARAM - revise o modelo")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
