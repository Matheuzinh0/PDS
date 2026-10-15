#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Autor: Matheus Sousa, LCSI-DES, UFPE.
"""
Aplica acentuacao correta do portugues SOMENTE no texto visivel de literais
de string de um arquivo .py (titulos de figura, legendas, mensagens, relatorio).

Salvaguardas:
  * preserva os trechos entre { } de f-strings (codigo Python: por exemplo
    {core.ganho_analogico_em(x)} nao pode virar ganho_analógico_em);
  * ignora literais "tecnicos" (nomes de arquivo/pasta/chave/estilo: contem
    apenas minusculas, digitos, '_', '.' ou '/');
  * ignora comentarios e identificadores.

Uso: python tools/acentuar_strings.py src/desafio08.py
"""
from __future__ import annotations

import re
import sys

MAPA = [
    ("prototipo", "protótipo"),
    ("Parametros", "Parâmetros"),
    ("parametros", "parâmetros"),
    ("Numericos", "Numéricos"),
    ("numericos", "numéricos"),
    ("Reflexao", "Reflexão"),
    ("reflexao", "reflexão"),
    ("Observacao", "Observação"),
    ("observacao", "observação"),
    ("Analogico", "Analógico"),
    ("analogico", "analógico"),
    ("analogica", "analógica"),
    ("maximo", "máximo"),
    ("minimo", "mínimo"),
    ("posicao", "posição"),
    ("Posicao", "Posição"),
    ("razao", "razão"),
    ("Razao", "Razão"),
    ("condicao", "condição"),
    ("frequencia", "frequência"),
    ("Frequencia", "Frequência"),
    ("replicas", "réplicas"),
    ("replica", "réplica"),
    ("ressonancia", "ressonância"),
    ("Ressonancia", "Ressonância"),
    ("invariancia", "invariância"),
    ("Invariancia", "Invariância"),
    ("estavel", "estável"),
    ("instavel", "instável"),
    ("unitario", "unitário"),
    ("propria", "própria"),
    ("proximo", "próximo"),
    ("diferenca", "diferença"),
    ("Diferenca", "Diferença"),
    ("mascara", "máscara"),
    ("relacao", "relação"),
    ("envoltoria", "envoltória"),
    ("logaritmica", "logarítmica"),
    ("discretizacao", "discretização"),
    ("Nao ", "Não "),
    ("nao ", "não "),
    ("Sao ", "São "),
    ("sao ", "são "),
    ("Analogo", "Análogo"),
    ("analogo", "análogo"),
]

TECNICO = re.compile(r"^[a-z0-9_./\-]*[._/][a-z0-9_./\-]*$")
LITERAL = re.compile(r'("(?:[^"\\]|\\.)*"|\'(?:[^\'\\]|\\.)*\')')


def acentua(texto: str) -> str:
    """Acentua o texto mantendo intactos os blocos {...} (expressoes)."""
    partes: list[tuple[str, str]] = []
    buf, prof = "", 0
    for ch in texto:
        if ch == "{":
            if prof == 0:
                partes.append(("txt", buf))
                buf = ""
            prof += 1
            buf += ch
        elif ch == "}":
            prof -= 1
            buf += ch
            if prof <= 0:
                prof = 0
                partes.append(("cod", buf))
                buf = ""
        else:
            buf += ch
    if buf:
        partes.append(("txt" if prof == 0 else "cod", buf))
    saida = []
    for tipo, s in partes:
        if tipo == "txt":
            for de, para in MAPA:
                s = s.replace(de, para)
        saida.append(s)
    return "".join(saida)


def processa_literal(lit: str) -> str:
    aspas = lit[0]
    corpo = lit[1:-1]
    if TECNICO.match(corpo):
        return lit
    if corpo.count(aspas) == 0:
        return aspas + acentua(corpo) + aspas
    return lit


def main(caminho: str) -> int:
    linhas = open(caminho, encoding="utf-8").read().splitlines(keepends=True)
    saida = []
    for linha in linhas:
        if linha.lstrip().startswith("#") or "MAPA = [" in linha:
            saida.append(linha)
            continue
        saida.append(LITERAL.sub(lambda m: processa_literal(m.group(0)), linha))
    open(caminho, "w", encoding="utf-8").write("".join(saida))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1]))
