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
    """
    Converte um padrão IUPAC (DNA) para expressão regular (regex).
    
    Suporta códigos IUPAC DNA (por ex. R=[AG], Y=[CT], N=[ACGT]).
    
    Args:
        padrao (str): Padrão IUPAC.
    
    Returns:
        str: Regex equivalente.
    
    Raises:
        TypeError: Se `padrao` for None.
        ValueError: Se o padrão for vazio ou contiver símbolos não suportados.
    
    Example:
        >>> iupac_para_regex("ATN")
        'AT[ACGT]'
    """
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
    """
    Procura ocorrências (índices iniciais) de um padrão IUPAC numa sequência.
    
    Args:
        sequencia (str): Sequência onde procurar.
        padrao (str): Padrão IUPAC.
    
    Returns:
        List[int]: Lista de posições iniciais (0-based) onde o padrão ocorre.
    
    Raises:
        TypeError: Se `sequencia` for None.
        ValueError: Propagado de `iupac_para_regex` se o padrão for inválido.
    
    Example:
        >>> procura_iupac("ATGCAATG", "ATG")
        [0, 5]
    """
    if sequencia is None:
        raise TypeError("sequencia não pode ser None")

    sequencia = sequencia.upper().strip()
    regex = iupac_para_regex(padrao)
    return [m.start() for m in re.finditer(regex, sequencia)]


# =========================================================
# PROSITE -> regex + procura
# =========================================================

def prosite_para_regex(prosite):
    """
    Converte um padrão PROSITE (subconjunto) para regex.
    
    Implementa um conjunto mínimo típico:
    - Remove '-' (separadores)
    - x/X -> '.' (qualquer)
    - {ABC} -> [^ABC] (negação)
    - (n) e (n,m) -> repetição {n} ou {n,m}
    
    Args:
        prosite (str): Padrão PROSITE.
    
    Returns:
        str: Regex equivalente.
    
    Raises:
        TypeError: Se `prosite` for None.
        ValueError: Se o padrão for vazio.
    
    Example:
        >>> prosite_para_regex("C-x(2)-C")
        'C.{2}C'
    """
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
    """
    Procura ocorrências (índices iniciais) de um padrão PROSITE numa sequência.
    
    Args:
        sequencia (str): Sequência onde procurar.
        prosite (str): Padrão PROSITE.
    
    Returns:
        List[int]: Lista de posições iniciais (0-based).
    
    Raises:
        TypeError: Se `sequencia` for None.
        ValueError: Propagado de `prosite_para_regex` se o padrão for inválido.
    
    Example:
        >>> procura_prosite("ACCC", "C-x(2)-C")
        [1]
    """
    if sequencia is None:
        raise TypeError("sequencia não pode ser None")

    sequencia = sequencia.upper().strip()
    regex = prosite_para_regex(prosite)
    return [m.start() for m in re.finditer(regex, sequencia)]


# =========================================================
# ENZIMAS DE RESTRIÇÃO (DNA)
# =========================================================

def digestao_dna(sequencia, sitio_restricao):
    """
    Executa digestão de DNA para um sítio de restrição com posição de corte.
    
    O sítio deve conter '^' para indicar o ponto de corte (ex.: 'G^AATTC').
    O motivo pode conter códigos IUPAC (DNA).
    
    Args:
        sequencia (str): Sequência de DNA a digerir.
        sitio_restricao (str): Sítio de restrição com '^'.
    
    Returns:
        Tuple[List[int], List[str]]: (cortes, fragmentos), onde `cortes` são as posições
        de corte (0-based) e `fragmentos` é a lista de fragmentos resultantes.
    
    Raises:
        TypeError: Se `sequencia` ou `sitio_restricao` forem None.
        ValueError: Se não existir '^', se o motivo for vazio, se a posição for inválida,
        ou se o motivo contiver símbolos IUPAC inválidos.
    
    Example:
        >>> digestao_dna("GAATTCGAATTC", "G^AATTC")[0]
        [1, 7]
    """
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
    """
    Calcula contagens por coluna para um conjunto de sequências alinhadas.
    
    As sequências têm de ter o mesmo comprimento. Aplica pseudocontagens
    (Laplace) adicionando `pseudocontagem` a cada símbolo do alfabeto em cada coluna.
    
    Args:
        seqs (List[str]): Lista de sequências (motifs) do mesmo tamanho.
        alfabeto (str): Símbolos permitidos. Default: "ACGT".
        pseudocontagem (int): Valor a adicionar a cada contagem. Default: 1.
    
    Returns:
        List[dict]: Lista de dicionários (um por posição), {símbolo: contagem}.
    
    Raises:
        ValueError: Se `seqs` for vazio ou se houver comprimentos diferentes.
    
    Example:
        >>> tabela_contagens(["AA", "AT"], alfabeto="AT", pseudocontagem=1)
        [{'A': 3, 'T': 1}, {'A': 2, 'T': 2}]
    """
    if not seqs:
        raise ValueError("Lista de sequências vazia")
    if any(len(seqs[0]) != len(s) for s in seqs):
        raise ValueError("As sequências não têm todas o mesmo tamanho!")

    return [
        {b: ocorrencias.count(b) + pseudocontagem for b in alfabeto}
        for ocorrencias in zip(*seqs)
    ]

def pwm(seqs, tipo="DNA", pseudocontagem=1):
    """
    Constrói uma PWM (Position Weight Matrix) de probabilidades.
    
    As probabilidades são calculadas por coluna usando pseudocontagens:
        p = (contagem + pseudo) / (N + \|alfabeto\| * pseudo)
    
    Args:
        seqs (List[str]): Lista de sequências do mesmo comprimento.
        tipo (str): "DNA" ou "PROTEIN". Default: "DNA".
        pseudocontagem (int): Pseudocontagem por símbolo e por coluna. Default: 1.
    
    Returns:
        List[dict]: PWM como lista de dicionários (um por coluna), {símbolo: prob}.
    
    Raises:
        ValueError: Se `seqs` for vazio, comprimentos diferentes, ou `tipo` inválido.
    
    Example:
        >>> m = pwm(["AT", "AA"], tipo="DNA", pseudocontagem=1)
        >>> round(m[0]["A"], 2) >= 0.5
        True
    """
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
    """
    Arredonda os valores de uma PWM para apresentação.
    
    Args:
        matriz_pwm (List[dict]): PWM (lista de colunas).
        casas (int): Número de casas decimais. Default: 2.
    
    Returns:
        List[dict]: PWM com valores arredondados.
    
    Example:
        >>> imprime_pwm([{'A': 0.3333, 'C': 0.6666}], casas=2)
        [{'A': 0.33, 'C': 0.67}]
    """
    return [{k: round(v, casas) for k, v in col.items()} for col in matriz_pwm]

def prob_gerar_sequencia(sequencia, matriz_pwm):
    """
    Calcula a probabilidade de uma sequência segundo uma PWM.
    
    Multiplica as probabilidades coluna a coluna.
    
    Args:
        sequencia (str): Sequência com comprimento igual ao motif (len(PWM)).
        matriz_pwm (List[dict]): PWM (lista de colunas).
    
    Returns:
        float: Probabilidade da sequência.
    
    Raises:
        ValueError: Se o comprimento da sequência não coincidir com a PWM.
    
    Example:
        >>> m = [{'A': 1.0, 'C': 0.0}, {'A': 0.5, 'C': 0.5}]
        >>> prob_gerar_sequencia('AA', m)
        0.5
    """
    if len(sequencia) != len(matriz_pwm):
        raise ValueError("Tamanho da sequência e do motif não são iguais!")

    prob = 1.0
    for letra, coluna in zip(sequencia.upper(), matriz_pwm):
        prob *= coluna[letra]
    return prob

def seq_mais_provavel(sequencia, matriz_pwm):
    """
    Encontra a(s) subsequência(s) mais provável(eis) segundo uma PWM.
    
    Usa janelas deslizantes de tamanho len(PWM) e devolve todas as janelas
    com probabilidade máxima (sem duplicados).
    
    Args:
        sequencia (str): Sequência maior/igual ao tamanho do motif.
        matriz_pwm (List[dict]): PWM.
    
    Returns:
        List[str]: Lista ordenada de janelas mais prováveis.
    
    Raises:
        ValueError: Se a sequência for mais curta do que o motif.
    
    Example:
        >>> m = pwm(['AAA', 'AAT'], tipo='DNA', pseudocontagem=1)
        >>> seq_mais_provavel('CAAAT', m)[0]
        'AAA'
    """
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
    """
    Converte uma PWM em PSSM (log-odds, base 2).
    
    Assume background uniforme: bg = 1/\|alfabeto\|.
    
    Args:
        matriz_pwm (List[dict]): PWM.
        alfabeto (str): Alfabeto considerado no background. Default: "ACGT".
    
    Returns:
        List[dict]: PSSM como lista de colunas {símbolo: score_log2}.
    
    Example:
        >>> pssm = pssm_de_pwm([{'A': 0.5, 'C': 0.5}], alfabeto='AC')
        >>> round(pssm[0]['A'], 3)
        0.0
    """
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
    """
    Calcula o score (soma) de um k-mer segundo uma PSSM.
    
    Args:
        pssm (List[dict]): PSSM (lista de colunas).
        kmer (str): K-mer com comprimento igual ao número de colunas.
    
    Returns:
        float: Score total (soma dos scores por posição).
    
    Raises:
        ValueError: Se o comprimento do k-mer não coincidir com a PSSM.
    
    Example:
        >>> score_kmer([{'A': 1.0}, {'A': 2.0}], 'AA')
        3.0
    """
    if len(kmer) != len(pssm):
        raise ValueError("O comprimento do kmer deve coincidir com a PSSM")

    s = 0.0
    for i, b in enumerate(kmer.upper()):
        s += pssm[i].get(b, float("-inf"))
    return s

def melhor_subsequencia(pssm, sequencia):
    """
    Devolve a melhor subsequência segundo uma PSSM.
    
    Percorre todas as janelas de tamanho len(PSSM) e escolhe a de score máximo.
    
    Args:
        pssm (List[dict]): PSSM.
        sequencia (str): Sequência onde procurar.
    
    Returns:
        Tuple[str, int, float]: (melhor_kmer, posição_inicial, score).
    
    Raises:
        ValueError: Se a PSSM for vazia ou se a sequência for mais curta que o motif.
    
    Example:
        >>> melhor_subsequencia([{'A': 1.0}, {'A': 1.0}], 'CAAAT')
        ('AA', 2, 2.0)
    """
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
