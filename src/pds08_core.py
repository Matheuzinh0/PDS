#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Autor: Matheus Sousa, LCSI-DES, UFPE.
"""
===============================================================================
 PDS 2026.2 - DESAFIO 08
 Discretizacao de Filtros IIR por Invariancia ao Impulso e
 Analise de Aliasing Espectral
===============================================================================

 Nucleo de calculo (funcoes puras, SEM graficos e SEM I/O).

 Prototipo analogico passa-faixa do enunciado:

                 s                         s
  Ha(s) = -------------------  =  ---------------------------
          s^2 + 200*pi*s + (12000*pi)^2    s^2 + 2*zeta*wn*s + wn^2

  com  wn = 12000*pi rad/s  (f_n = 6000 Hz)
       2*zeta*wn = 200*pi  ->  zeta = 1/120  (Q = 1/(2*zeta) = 60)

  Metodo da invariancia ao impulso, Ts = 25 us (Fs = 40 kHz):
       h[n] = Ts * ha(n*Ts)
       H(z) = sum_n h[n] z^-n

  Item 1 do desafio  -> ha(t)  = L^-1{Ha(s)}
  Item 2 do desafio  -> coeficientes exatos de H(z) (2a ordem)
  Item 3 do desafio  -> ganho de pico em omega0 = 2*pi*6000*Ts

 Referencias das equacoes: Oppenheim & Schafer, "Discrete-Time Signal
 Processing", cap. 7 (projeto por invariancia ao impulso).
===============================================================================
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# --------------------------------------------------------------------------
# 1. PARAMETROS DO PROTOTIPO ANALOGICO
# --------------------------------------------------------------------------

PI = np.pi

WN = 12000.0 * PI            # frequencia natural do prototipo, rad/s
TWO_ZETA_WN = 200.0 * PI     # coeficiente do termo de 1a ordem, rad/s

FS_NOMINAL = 40_000.0        # Hz  -> Ts = 25 us  (enunciado)
TS_NOMINAL = 1.0 / FS_NOMINAL


@dataclass(frozen=True)
class ParametrosAnalogicos:
    """Parametros do prototipo analogico de 2a ordem (polo complexo conjugado)."""
    wn: float        # frequencia natural, rad/s
    zeta: float      # fator de amortecimento
    Q: float         # fator de qualidade = 1/(2*zeta)
    sigma: float     # parte real dos polos = zeta*wn, rad/s
    wd: float        # parte imaginaria dos polos = wn*sqrt(1-zeta^2), rad/s
    fn: float        # frequencia natural, Hz
    bw: float        # banda de -3 dB aproximada, Hz


def parametros_analogicos() -> ParametrosAnalogicos:
    """Calcula wn, zeta, sigma e wd a partir de Ha(s) do enunciado."""
    zeta = TWO_ZETA_WN / (2.0 * WN)
    sigma = zeta * WN
    wd = np.sqrt(WN ** 2 - sigma ** 2)          # = 100*pi*sqrt(14399)
    return ParametrosAnalogicos(
        wn=WN,
        zeta=zeta,
        Q=1.0 / (2.0 * zeta),
        sigma=sigma,
        wd=wd,
        fn=WN / (2.0 * PI),
        bw=WN / (2.0 * PI) / (1.0 / (2.0 * zeta)),
    )


# --------------------------------------------------------------------------
# 2. ITEM 1 - RESPOSTA AO IMPULSO ANALOGICA ha(t)
# --------------------------------------------------------------------------

def ha_analogica(t: np.ndarray | float) -> np.ndarray:
    """
    ha(t) = L^-1{ Ha(s) } = e^(-sigma*t) * [ cos(wd*t) - (sigma/wd)*sen(wd*t) ]

    Deduzida por decomposicao em fracoes parciais:

        s                       (s + sigma)              sigma
    -------------  =  ----------------------------  -  --------------
    (s+sigma)^2+wd^2        (s+sigma)^2 + wd^2      (s+sigma)^2 + wd^2

    cujas inversas sao  e^(-sigma t) cos(wd t)  e  e^(-sigma t) sen(wd t).

    Validada simbolicamente com SymPy em ``verificacao.py`` e comparada
    numericamente com a transformada inversa de Laplace numerica.
    """
    p = parametros_analogicos()
    t = np.asarray(t, dtype=float)
    return np.exp(-p.sigma * t) * (
        np.cos(p.wd * t) - (p.sigma / p.wd) * np.sin(p.wd * t)
    )


def ha_analogica_tempo_continuo_equivalente(t: np.ndarray | float) -> np.ndarray:
    """Mesma ha(t), escrita na forma A*e^(p1 t) + A* e^(p2 t) (util p/ conferencia)."""
    p = parametros_analogicos()
    t = np.asarray(t, dtype=float)
    p1 = -p.sigma + 1j * p.wd
    A = p1 / (p1 - np.conj(p1))
    return 2.0 * np.real(A * np.exp(p1 * t))


# --------------------------------------------------------------------------
# 3. ITEM 2 - DISCRETIZACAO POR INVARIANCIA AO IMPULSO  (h[n] = Ts*ha(nTs))
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class FiltroDigital:
    """Filtro digital de 2a ordem H(z) = (b0 + b1 z^-1)/(1 + a1 z^-1 + a2 z^-2)."""
    Fs: float
    Ts: float
    b: np.ndarray        # [b0, b1]
    a: np.ndarray        # [1, a1, a2]
    r: float             # raio dos polos  = exp(-sigma*Ts)
    theta: float         # angulo dos polos = wd*Ts  [rad]
    w0_raio: float       # omega0 = 2*pi*6000*Ts [rad/amostra]
    f0: float            # 6000 Hz (frequencia do prototipo)

    @property
    def polos(self) -> np.ndarray:
        return np.roots(self.a)

    @property
    def zeros(self) -> np.ndarray:
        return np.roots(self.b)

    def reais_em(self, omega: np.ndarray | float) -> np.ndarray:
        """
        H(e^jw) exatamente, avaliando o polinomio racional em z = e^(jw).

        ATENCAO ao ordenamento: H(z) e um polinomio em z^-1, logo
        np.polyval precisa dos coeficientes em ordem decrescente de grau em
        z^-1, isto e, os coeficientes invertidos:  b0 + b1*z^-1  ->
        polyval([b1, b0], z^-1).
        """
        w = np.atleast_1d(np.asarray(omega, dtype=float))
        zinv = np.exp(-1j * w)
        return np.polyval(self.b[::-1], zinv) / np.polyval(self.a[::-1], zinv)


def cria_zeros(n: int) -> np.ndarray:
    """Vetor de zeros (alocacao explicita - exigida pelo numpy moderno)."""
    return np.zeros(n, dtype=float)


def discretiza_invariancia_impulso(Fs: float = FS_NOMINAL) -> FiltroDigital:
    """
    Metodo da invariancia ao impulso aplicado a Ha(s) do enunciado.

    Deduzcao (usada no codigo e nos slides)
    ---------------------------------------
    Com sigma = zeta*wn, wd = wn*sqrt(1-zeta^2) e Ts = 1/Fs:

        h[n] = Ts*ha(nTs)
             = Ts*e^(-sigma*n*Ts)*cos(wd*n*Ts)
               - Ts*(sigma/wd)*e^(-sigma*n*Ts)*sen(wd*n*Ts)

    Definindo  r = e^(-sigma*Ts)  e  theta = wd*Ts, as transformadas Z
    tabeladas dos dois termos levam a:

                   Ts*[ 1 - r*(cos(theta) + (sigma/wd)*sen(theta)) * z^-1 ]
        H(z) = ---------------------------------------------------------------
                      1 - 2*r*cos(theta) * z^-1 + r^2 * z^-2

    Portanto:   b0 = Ts
                b1 = -Ts*r*( cos(theta) + (sigma/wd)*sen(theta) )
                a1 = -2*r*cos(theta)
                a2 = r^2
    """
    p = parametros_analogicos()
    Ts = 1.0 / Fs
    r = np.exp(-p.sigma * Ts)
    theta = p.wd * Ts

    b = np.array(
        [
            Ts,
            -Ts * r * (np.cos(theta) + (p.sigma / p.wd) * np.sin(theta)),
        ]
    )
    a = np.array([1.0, -2.0 * r * np.cos(theta), r * r])

    return FiltroDigital(
        Fs=Fs, Ts=Ts, b=b, a=a, r=r, theta=theta,
        w0_raio=2.0 * PI * 6000.0 * Ts, f0=6000.0,
    )


def resposta_impulso_digital(Fs: float = FS_NOMINAL, n_amostras: int = 2000):
    """
    h[n] = Ts*ha(n*Ts) amostrado diretamente da expressao analitica
    (usado para conferir os coeficientes de H(z) e para os graficos).
    """
    Ts = 1.0 / Fs
    n = np.arange(n_amostras)
    return Ts * ha_analogica(n * Ts)


def dtft_de_amostras(h: np.ndarray, omega: np.ndarray,
                     Ts: float = 1.0) -> np.ndarray:
    """
    TDFT (DTFT) numerica da sequencia h[n] por soma direta truncada:

        H(e^jw) ~= sum_{n=0}^{N-1} h[n] e^(-jwn)

    Usada como verificacao INDEPENDENTE da formula racional fechada.
    (Ts serve apenas para o rotulo; a escala ja esta em h[n].)
    """
    omega = np.atleast_1d(np.asarray(omega, dtype=float))
    # matriz (len(omega), len(h)) - limita memoria em malhas grandes
    return (h[None, :] * np.exp(-1j * omega[:, None] * np.arange(h.size))).sum(axis=1)


# --------------------------------------------------------------------------
# 4. ITEM 3 - GANHO DE PICO E ALIASING
# --------------------------------------------------------------------------

def ganho_analogico_em(omega_rad_s: np.ndarray | float) -> np.ndarray:
    """|Ha(j*Omega)| para Omega em rad/s."""
    p = parametros_analogicos()
    Om = np.atleast_1d(np.asarray(omega_rad_s, dtype=float))
    num = np.abs(Om)
    den = np.sqrt((p.wn ** 2 - Om ** 2) ** 2 + (TWO_ZETA_WN * Om) ** 2)
    return num / den


def pico_analogico() -> dict:
    """Pico exato de |Ha(j*Omega)|: em Omega_r = wn*sqrt(1-2*zeta^2)."""
    p = parametros_analogicos()
    Om_r = p.wn * np.sqrt(max(1.0 - 2.0 * p.zeta ** 2, 0.0))
    gmax = ganho_analogico_em(Om_r).item()
    return {"Omega_r": Om_r, "f_r": Om_r / (2 * PI), "ganho": gmax,
            "ganho_dB": 20 * np.log10(gmax),
            "ganho_em_wn": ganho_analogico_em(p.wn).item()}


def pico_digital(f: FiltroDigital, n_pontos: int = 4_000_001) -> dict:
    """
    Pico de |H(e^jw)| varrendo a banda digital [0, pi] com malha fina,
    alem do valor exato na frequencia pedida w0 = 2*pi*6000*Ts.
    """
    w = np.linspace(0.0, PI, n_pontos)
    mag = np.abs(f.reais_em(w))
    k = int(np.argmax(mag))
    return {
        "w_pico": w[k],
        "f_pico": w[k] * f.Fs / (2 * PI),
        "ganho_pico": mag[k],
        "ganho_pico_dB": 20 * np.log10(mag[k]),
        "ganho_em_w0": np.abs(f.reais_em(f.w0_raio)).item(),
        "ganho_em_w0_dB": (20 * np.log10(np.abs(f.reais_em(f.w0_raio)))).item(),
        "ganho_DC": np.abs(f.reais_em(0.0)).item(),
    }


def soma_de_replicas(f: FiltroDigital, omega: np.ndarray,
                     k_max: int = 12, correcao_meia_descontinuidade: bool = True
                     ) -> np.ndarray:
    """
    Soma periodica (Poisson) que explica o aliasing frequencial:

        H(e^jw)  =  sum_k  Ha( j*(w - 2*pi*k)/Ts )  +  Ts*ha(0+)/2
                    \\_______________________________/     \\_________/
                       replicas deslocadas de k*Fs      termo de meia-
                                                        descontinuidade

    ATENCAO (resultado fino, comprovado numericamente em verificacao.py):
    a soma classica de replicas SOZINHA nao reproduz H(e^jw) quando ha(t)
    possui descontinuidade de salto em t = 0. Como
    ha(0+) = 1 e ha(0-) = 0, a formula de Poisson aplicada a uma funcao
    de variacao limitada com salto valem as replicas MAIS metade do salto
    amostrado, Ts*ha(0+)/2 = Ts/2. Sem esse termo a soma das replicas
    fica exatamente Ts/2 abaixo da resposta exata (constante em w).
    Isso decorre do decaimento lento |Ha(j*Omega)| ~ 1/Omega (mesma razao
    pela qual a identidade falha de forma catastrofica no passa-altas).

    A convergencia e apenas CONDICIONAL (termos ~1/k com cancelamento
    par a par), por isso k_max deve ser grande para boa precisao.
    """
    w = np.atleast_1d(np.asarray(omega, dtype=float))
    total = np.zeros_like(w, dtype=complex)
    p = parametros_analogicos()
    for k in range(-k_max, k_max + 1):
        Om = (w - 2.0 * PI * k) / f.Ts
        # resposta complexa exata de Ha(s) em s = j*Om
        Ha = 1j * Om / ((1j * Om) ** 2 + TWO_ZETA_WN * (1j * Om) + p.wn ** 2)
        total += Ha
    if correcao_meia_descontinuidade:
        total = total + f.Ts * float(ha_analogica(0.0)) / 2.0
    return total


def banda_meia_potencia(f: FiltroDigital) -> dict:
    """Banda de -3 dB do filtro digital (a partir do pico) em Hz e rad/amostra."""
    w = np.linspace(0.0, PI, 2_000_001)
    mag = np.abs(f.reais_em(w))
    kp = int(np.argmax(mag))
    lim = mag[kp] / np.sqrt(2.0)
    # procura de cada lado do pico
    esq = np.where(mag[:kp] <= lim)[0]
    dir_ = np.where(mag[kp:] <= lim)[0]
    w_esq = w[esq[-1]] if esq.size else 0.0
    w_dir = w[kp + dir_[0]] if dir_.size else PI
    return {
        "w_pico": w[kp],
        "w_3dB_baixo": w_esq,
        "w_3dB_alto": w_dir,
        "df_3dB_Hz": (w_dir - w_esq) * f.Fs / (2 * PI),
    }


# --------------------------------------------------------------------------
# 5. PARTE COMPUTACIONAL EXTRA - FALHA DO METODO NO PASSA-ALTAS
# --------------------------------------------------------------------------

HP_A_CONST = 2000.0 * PI     # Ha_hp(s) = s/(s + 2000*pi)  -> fc = 1 kHz


def ha_passa_altas(t: np.ndarray | float):
    """
    Ha_hp(s) = s/(s + a),  a = 2000*pi rad/s:

        ha_hp(t) = delta(t) - a*e^(-a*t)*u(t)

    O termo delta(t) e a origem da degeneracao: ele NAO tem versao
    amostrada que preserve a resposta em frequencia do passa-altas.
    """
    a = HP_A_CONST
    t = np.asarray(t, dtype=float)
    return -a * np.exp(-a * t)


def filtro_passa_altas_invariancia(Fs: float = FS_NOMINAL) -> FiltroDigital:
    """
    Invariancia ao impulso aplicada ao passa-altas.

    h[n] = Ts*ha_hp(nTs) = Ts*delta[n] - Ts*a*e^(-a*n*Ts)*u[n]
         = Ts*delta[n] - A*alpha^n*u[n],      A = a*Ts,  alpha = e^(-a*Ts)

    H(z) = Ts - A/(1 - alpha*z^-1)
         = [ (Ts - A) - Ts*alpha*z^-1 ] / (1 - alpha*z^-1)
    """
    Ts = 1.0 / Fs
    a = HP_A_CONST
    A = a * Ts
    alpha = np.exp(-a * Ts)
    b = np.array([Ts - A, -Ts * alpha])
    a_ = np.array([1.0, -alpha])
    omega1 = a * Ts                        # ... reaproveitado apenas p/ rotulo
    return FiltroDigital(Fs=Fs, Ts=Ts, b=b, a=a_, r=alpha, theta=omega1,
                         w0_raio=2 * PI * 1000.0 * Ts, f0=1000.0)


def resposta_analogica_passa_altas(omega_rad_s: np.ndarray) -> np.ndarray:
    """|Ha_hp(j*Omega)| = Omega/sqrt(Omega^2 + a^2)  (passa-altas, ->1 em HF)."""
    Om = np.atleast_1d(np.asarray(omega_rad_s, dtype=float))
    a = HP_A_CONST
    return Om / np.sqrt(Om ** 2 + a ** 2)
