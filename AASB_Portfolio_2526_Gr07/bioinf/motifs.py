# motifs.py
import re
import math


# =========================================================
# IUPAC (DNA) -> regex + procura
# =========================================================

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]",
    "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]",
    "N": "[ACGT]"
}

def iupac_para_regex(padrao):
    """Converte um padrão IUPAC (DNA) para regex."""
    if padrao is None:
        raise TypeError("padrao não pode ser None")

    padrao = padrao.upper().strip()
    if not padrao:
        raise ValueError("Padrão vazio")

    for c in padrao:
        if c not in IUPAC:
            raise ValueError(f"Símbolo IUPAC inválido: {c}")

    return "".join(IUPAC[c] for c in padrao)

def procura_iupac(sequencia, padrao):
    """Procura um padrão IUPAC numa sequência."""
    if sequencia is None:
        raise TypeError("sequencia não pode ser None")

    sequencia = sequencia.upper().strip()
    regex = iupac_para_regex(padrao)
    return [m.start() for m in re.finditer(regex, sequencia)]


# =========================================================
# PROSITE -> regex + procura
# =========================================================

def prosite_para_regex(prosite):
    """Converte PROSITE (mínimo) para regex."""
    if prosite is None:
        raise TypeError("prosite não pode ser None")

    prosite = prosite.strip()
    if not prosite:
        raise ValueError("PROSITE vazio")

    p = prosite.replace("-", "")
    p = p.replace("x", ".").replace("X", ".")
    p = re.sub(r"\{([A-Z]+)\}", r"[^\1]", p)
    p = re.sub(r"\((\d+)\)", r"{\1}", p)
    p = re.sub(r"\((\d+),(\d+)\)", r"{\1,\2}", p)
    return p

def procura_prosite(sequencia, prosite):
    """Procura um padrão PROSITE numa sequência."""
    if sequencia is None:
        raise TypeError("sequencia não pode ser None")

    sequencia = sequencia.upper().strip()
    regex = prosite_para_regex(prosite)
    return [m.start() for m in re.finditer(regex, sequencia)]


# =========================================================
# ENZIMAS DE RESTRIÇÃO (DNA)
# =========================================================

def digestao_dna(sequencia, sitio_restricao):
    """Faz digestão de DNA com sítio com '^' e devolve cortes e fragmentos."""
    if sequencia is None or sitio_restricao is None:
        raise TypeError("sequencia/sitio_restricao não pode ser None")

    sequencia = sequencia.upper().strip()
    sitio_restricao = sitio_restricao.upper().strip()

    if "^" not in sitio_restricao:
        raise ValueError("O sítio de restrição deve conter '^'")

    pos_corte = sitio_restricao.index("^")
    motivo = sitio_restricao.replace("^", "")

    if not motivo:
        raise ValueError("Motivo de restrição vazio")
    if pos_corte < 0 or pos_corte > len(motivo):
        raise ValueError("Posição de corte fora do intervalo")

    regex = iupac_para_regex(motivo)

    inicios = [m.start() for m in re.finditer(regex, sequencia)]
    cortes = sorted(set([i + pos_corte for i in inicios]))

    fragmentos = []
    ultimo = 0
    for c in cortes:
        fragmentos.append(sequencia[ultimo:c])
        ultimo = c
    fragmentos.append(sequencia[ultimo:])

    return cortes, fragmentos


# =========================================================
# PWM (com pseudocontagem=1)
# =========================================================

def tabela_contagens(seqs, alfabeto="ACGT", pseudocontagem=1):
    """Conta símbolos por posição (com pseudocontagem)."""
    if not seqs:
        raise ValueError("Lista de sequências vazia")
    if any(len(seqs[0]) != len(s) for s in seqs):
        raise ValueError("As sequências não têm todas o mesmo tamanho!")

    return [
        {b: ocorrencias.count(b) + pseudocontagem for b in alfabeto}
        for ocorrencias in zip(*seqs)
    ]

def pwm(seqs, tipo="DNA", pseudocontagem=1):
    """Cria uma PWM (probabilidades) a partir de sequências."""
    if not seqs:
        raise ValueError("Lista de sequências vazia")
    if any(len(seqs[0]) != len(s) for s in seqs):
        raise ValueError("As sequências não têm todas o mesmo tamanho!")

    tipo = tipo.upper()
    if tipo not in ("DNA", "PROTEIN"):
        raise ValueError(f"Tipo inválido: {tipo}")

    alfabeto = "ACGT" if tipo == "DNA" else "ARNDCQEGHILKMFPSTWYVBZX_"
    tabela = tabela_contagens(seqs, alfabeto=alfabeto, pseudocontagem=pseudocontagem)

    n = len(seqs)
    a = len(alfabeto)
    denominador = n + a * pseudocontagem

    return [{k: v / denominador for k, v in coluna.items()} for coluna in tabela]

def imprime_pwm(matriz_pwm, casas=2):
    """Arredonda valores da PWM."""
    return [{k: round(v, casas) for k, v in col.items()} for col in matriz_pwm]

def prob_gerar_sequencia(sequencia, matriz_pwm):
    """Probabilidade de uma sequência segundo a PWM."""
    if len(sequencia) != len(matriz_pwm):
        raise ValueError("Tamanho da sequência e do motif não são iguais!")

    prob = 1.0
    for letra, coluna in zip(sequencia.upper(), matriz_pwm):
        prob *= coluna[letra]
    return prob

def seq_mais_provavel(sequencia, matriz_pwm):
    """Devolve a(s) janela(s) mais provável(eis) segundo a PWM."""
    sequencia = sequencia.upper().strip()
    l = len(matriz_pwm)

    if len(sequencia) < l:
        raise ValueError("A sequência tem que ser >= ao motif!")

    janelas = [sequencia[i:i+l] for i in range(0, len(sequencia) - l + 1)]
    probs = {s: prob_gerar_sequencia(s, matriz_pwm) for s in janelas}
    maior = max(probs.values())

    return sorted(set([s for s, p in probs.items() if p == maior]))


# =========================================================
# PSSM + scoring + melhor subsequência
# =========================================================

def pssm_de_pwm(matriz_pwm, alfabeto="ACGT"):
    """Converte PWM em PSSM (log2(p/bg))."""
    bg = 1.0 / len(alfabeto)
    pssm = []
    for coluna in matriz_pwm:
        d = {}
        for b in alfabeto:
            p = coluna.get(b, 0.0)
            d[b] = float("-inf") if p == 0 else math.log2(p / bg)
        pssm.append(d)
    return pssm

def score_kmer(pssm, kmer):
    """Score de um kmer segundo a PSSM."""
    if len(kmer) != len(pssm):
        raise ValueError("O comprimento do kmer deve coincidir com a PSSM")

    s = 0.0
    for i, b in enumerate(kmer.upper()):
        s += pssm[i].get(b, float("-inf"))
    return s

def melhor_subsequencia(pssm, sequencia):
    """Melhor subsequência segundo a PSSM (kmer, posição, score)."""
    sequencia = sequencia.upper().strip()
    k = len(pssm)

    if k == 0:
        raise ValueError("PSSM vazia")
    if len(sequencia) < k:
        raise ValueError("Sequência mais curta do que o motif")

    melhor = ""
    melhor_pos = -1
    melhor_score = float("-inf")

    for i in range(0, len(sequencia) - k + 1):
        kmer = sequencia[i:i+k]
        sc = score_kmer(pssm, kmer)
        if sc > melhor_score:
            melhor_score = sc
            melhor = kmer
            melhor_pos = i

    return melhor, melhor_pos, melhor_score
