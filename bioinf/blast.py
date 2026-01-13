

# blast.py
from typing import List, Dict, Tuple

# ---------- Helpers de extensão/score (baixa complexidade) ----------

def _score_seed(query: str, seq: str, i_query: int, i_seq: int,
                k: int, match: int, mismatch: int) -> int:
    """Pontua o k-mer seed."""
    score = 0
    for t in range(k):
        score += match if query[i_query + t] == seq[i_seq + t] else mismatch
    return score

def _extend_left(query: str, seq: str, i_query: int, i_seq: int,
                 match: int, mismatch: int) -> Tuple[int, int]:
    """Extensão à esquerda: devolve (delta_score, extensão_em_caracteres)."""
    limite = min(i_query, i_seq)
    delta = 0
    left = 0
    # percorre o número máximo seguro de passos para a esquerda
    for off in range(1, limite + 1):
        delta += match if query[i_query - off] == seq[i_seq - off] else mismatch
        left = off
    return delta, left

def _extend_right(query: str, seq: str, i_query: int, i_seq: int, k: int,
                  match: int, mismatch: int) -> Tuple[int, int]:
    """Extensão à direita: devolve (delta_score, extensão_em_caracteres)."""
    # número máximo de passos sem sair da string
    limite = min(len(query) - (i_query + k), len(seq) - (i_seq + k))
    delta = 0
    right = 0
    for off in range(limite):
        delta += match if query[i_query + k + off] == seq[i_seq + k + off] else mismatch
        right = off + 1
    return delta, right

# ---------- Funções públicas ----------

def gerar_palavras(seq: str, k: int = 3) -> List[str]:
    """
    Divide uma sequência em palavras (k-mers).
    
    Args:
        seq (str): Sequência (query ou subject) a partir da qual são gerados k-mers.
        k (int): Tamanho do k-mer. Default: 3.
    
    Returns:
        List[str]: Lista de k-mers consecutivos (janela deslizante).
    
    Example:
        >>> gerar_palavras("ATGCA", k=3)
        ['ATG', 'TGC', 'GCA']
    """
    return [seq[i:i+k] for i in range(len(seq) - k + 1)]

def indexar_banco(banco: List[str], k: int = 3) -> Dict[str, List[Tuple[int, int]]]:
    """
    Cria um índice de k-mers para um banco de sequências.
    
    O índice mapeia cada k-mer para uma lista de ocorrências (id da sequência no banco,
    posição inicial do k-mer nessa sequência).
    
    Args:
        banco (List[str]): Lista de sequências (subjects).
        k (int): Tamanho do k-mer. Default: 3.
    
    Returns:
        Dict[str, List[Tuple[int, int]]]: Dicionário {k-mer: [(id_seq, pos), ...]}.
    
    Example:
        >>> indexar_banco(["ATGCA"], k=3)["TGC"]
        [(0, 1)]
    """
    indice: Dict[str, List[Tuple[int, int]]] = {}
    for id_seq, seq in enumerate(banco):
        for i in range(len(seq) - k + 1):
            palavra = seq[i:i+k]
            indice.setdefault(palavra, []).append((id_seq, i))
    return indice

def encontrar_hits(query: str, indice: Dict[str, List[Tuple[int, int]]], k: int = 3) -> List[Tuple[int, int, int]]:
    """
    Encontra hits de k-mers da query num índice de banco.
    
    Para cada k-mer da query, devolve todas as ocorrências no banco.
    
    Args:
        query (str): Sequência query.
        indice (Dict[str, List[Tuple[int, int]]]): Índice produzido por `indexar_banco`.
        k (int): Tamanho do k-mer. Default: 3.
    
    Returns:
        List[Tuple[int, int, int]]: Lista de hits (pos_query, id_seq, pos_seq).
    
    Example:
        >>> ind = indexar_banco(["ATGCA"], k=3)
        >>> encontrar_hits("TTGC", ind, k=3)
        [(1, 0, 1)]
    """
    hits: List[Tuple[int, int, int]] = []
    for i, palavra in enumerate(gerar_palavras(query, k)):
        for id_seq, pos in indice.get(palavra, []):
            hits.append((i, id_seq, pos))
    return hits

def estender_hit(query: str, seq: str, i_query: int, i_seq: int,
                 k: int = 3, match: int = 2, mismatch: int = -1) -> Tuple[int, str, str]:
    """
    Estende um hit (seed) para ambos os lados e calcula o score.
    
    A extensão é feita com base em match/mismatch (não usa gaps). Devolve o alinhamento
    local resultante como duas strings (substrings alinhadas por posição).
    
    Args:
        query (str): Sequência query.
        seq (str): Sequência do banco (subject).
        i_query (int): Índice inicial do hit na query.
        i_seq (int): Índice inicial do hit no subject.
        k (int): Tamanho do seed (k-mer). Default: 3.
        match (int): Score para match. Default: 2.
        mismatch (int): Score para mismatch. Default: -1.
    
    Returns:
        Tuple[int, str, str]: (score_total, query_alinhada, subject_alinhada).
    
    Example:
        >>> estender_hit("ATGCA", "TTGCA", i_query=1, i_seq=1, k=3)
        (6, 'TGCA', 'TGCA')
    """
    # pontuar o seed k-mer
    score = _score_seed(query, seq, i_query, i_seq, k, match, mismatch)
    # extensões esquerda e direita com loops de limites fixos (sem while)
    delta_esq, left = _extend_left(query, seq, i_query, i_seq, match, mismatch)
    delta_dir, right = _extend_right(query, seq, i_query, i_seq, k, match, mismatch)
    score += (delta_esq + delta_dir)

    q_aln = query[i_query - left : i_query + k + right]
    s_aln = seq[i_seq   - left : i_seq   + k + right]
    return score, q_aln, s_aln

def blast_simplificado(query: str, banco: List[str], k: int = 3) -> List[Tuple[int, int, str, str]]:
    """
    Executa um BLAST simplificado baseado em seeds (k-mers) e extensão sem gaps.
    
    Pipeline:
    1) Indexa o banco por k-mers.
    2) Encontra hits de k-mers da query no índice.
    3) Estende cada hit e calcula o score.
    4) Ordena resultados por score (descendente).
    
    Args:
        query (str): Sequência query.
        banco (List[str]): Lista de sequências do banco (subjects).
        k (int): Tamanho do k-mer. Default: 3.
    
    Returns:
        List[Tuple[int, int, str, str]]: Lista de resultados (score, id_seq, q_aln, s_aln)
        ordenada por score desc.
    
    Example:
        >>> blast_simplificado("ATGCA", ["TTGCA", "AAAAA"], k=3)[0][1]
        0
    """
    indice = indexar_banco(banco, k)
    hits = encontrar_hits(query, indice, k)
    melhores: List[Tuple[int, int, str, str]] = []
    for i_query, id_seq, i_seq in hits:
        score, q_aln, s_aln = estender_hit(query, banco[id_seq], i_query, i_seq, k)
        melhores.append((score, id_seq, q_aln, s_aln))
    return sorted(melhores, key=lambda x: x[0], reverse=True)
