

# filo.py
from __future__ import annotations
from typing import Dict, List, Tuple
from itertools import combinations

__all__ = [
    "calcular_matriz_distancias",
    "upgma",
]

# ---------------------------------------------------------------------
# Helpers de distância (baixam complexidade das funções públicas)
# ---------------------------------------------------------------------

def _pdist(a: str, b: str) -> float:
    """
    p-distância: proporção de posições diferentes entre duas sequências.

    Args:
        a (str): Sequência 1.
        b (str): Sequência 2.

    Returns:
        float: Proporção de diferenças em [0, 1].
    """
    n = len(a)
    if n == 0:
        return 0.0
    dif = sum(c1 != c2 for c1, c2 in zip(a, b))
    return dif / n


def _validate_equal_lengths(sequencias: List[str]) -> None:
    """
    Valida que a lista não é vazia e que todas as sequências têm o mesmo tamanho.

    Args:
        sequencias (List[str]): Lista de sequências.

    Raises:
        ValueError: Se a lista for vazia ou se houver tamanhos distintos.
    """
    if not sequencias:
        raise ValueError("A lista de sequências não pode ser vazia.")
    size = len(sequencias[0])
    if any(len(s) != size for s in sequencias):
        raise ValueError("Todas as sequências devem ter o mesmo tamanho.")


# ---------------------------------------------------------------------
# API pública: Matriz de Distâncias
# ---------------------------------------------------------------------

def calcular_matriz_distancias(sequencias: List[str]) -> Dict[str, Dict[str, float]]:
    """
    Calcula a matriz de distâncias usando **p-distância**.

    Usa cada sequência como etiqueta (chave do dicionário). Para projetos
    reais em que possa haver sequências iguais com nomes diferentes, é
    recomendável passar nomes explicitamente — aqui segue o enunciado.

    Args:
        sequencias (List[str]): Lista de sequências (todas com o mesmo tamanho).

    Returns:
        Dict[str, Dict[str, float]]: Matriz de distâncias simétrica (0 na diagonal).

    Raises:
        ValueError: Se a lista for vazia ou se os comprimentos diferirem.
    """
    _validate_equal_lengths(sequencias)
    nomes = list(sequencias)

    matriz: Dict[str, Dict[str, float]] = {n: {} for n in nomes}
    # Diagonal
    for a in nomes:
        matriz[a][a] = 0.0
    # Metade superior (e espelhar para a inferior)
    for a, b in combinations(nomes, 2):
        d = _pdist(a, b)
        matriz[a][b] = d
        matriz[b][a] = d
    return matriz


# ---------------------------------------------------------------------
# Helpers de UPGMA
# ---------------------------------------------------------------------

def _closest_pair(matriz: Dict[str, Dict[str, float]],
                  clusters: Dict[str, Tuple[str, float]]) -> Tuple[str, str, float]:
    """
    Encontra o par de clusters com menor distância.

    Returns:
        Tuple[str, str, float]: (a, b, dmin) com a <-> b o par mais próximo.
    """
    best_a = best_b = ""
    best = float("inf")
    for a, b in combinations(clusters.keys(), 2):
        d = matriz[a][b]
        if d < best:
            best = d
            best_a, best_b = a, b
    return best_a, best_b, best


def _newick_join(a: str, b: str, dist_ab: float,
                 clusters: Dict[str, Tuple[str, float]]) -> Tuple[str, float]:
    """
    Une dois clusters em Newick e calcula a nova altura (UPGMA).

    Args:
        a, b (str): Clusters a unir.
        dist_ab (float): Distância entre a e b.
        clusters: Mapa etiqueta -> (newick, altura).

    Returns:
        Tuple[str, float]: (novo_newick, nova_altura).
    """
    altura = dist_ab / 2.0
    na, ha = clusters[a]
    nb, hb = clusters[b]
    ramo_a = max(0.0, altura - ha)
    ramo_b = max(0.0, altura - hb)
    return f"({na}:{ramo_a:.6f},{nb}:{ramo_b:.6f})", altura


def _average_linkage_row(matriz: Dict[str, Dict[str, float]],
                         a: str, b: str, novo: str,
                         restantes: List[str]) -> None:
    """
    Cria as distâncias do novo cluster para os restantes via média aritmética.
    """
    matriz[novo] = {}
    for k in restantes:
        d = (matriz[a][k] + matriz[b][k]) / 2.0
        matriz[novo][k] = d
        matriz[k][novo] = d
    matriz[novo][novo] = 0.0


def _remove_old_clusters(matriz: Dict[str, Dict[str, float]],
                         a: str, b: str) -> None:
    """
    Remove linhas/colunas dos clusters a e b da matriz.
    """
    for x in (a, b):
        matriz.pop(x, None)
        for linha in matriz.values():
            linha.pop(x, None)


# ---------------------------------------------------------------------
# API pública: UPGMA
# ---------------------------------------------------------------------

def upgma(matriz: Dict[str, Dict[str, float]]) -> str:
    """
    Constrói uma árvore filogenética com **UPGMA** e retorna em **Newick**.

    Args:
        matriz (Dict[str, Dict[str, float]]): Matriz de distâncias simétrica
            com 0 na diagonal (ex.: saída de `calcular_matriz_distancias`).

    Returns:
        str: Árvore em formato Newick terminada com `;`.

    Notes:
        - Caso trivial: 0 clusters → retorna apenas `;`.
        - Caso trivial: 1 cluster  → retorna a etiqueta + `;`.
        - As etiquetas usadas são as chaves da matriz (neste projeto, as próprias sequências).
    """
    # clusters: etiqueta -> (newick, altura)
    clusters: Dict[str, Tuple[str, float]] = {name: (name, 0.0) for name in list(matriz.keys())}

    # Casos triviais
    if not clusters:
        return ";"
    if len(clusters) == 1:
        return f"{next(iter(clusters))};"

    # Itera até restar um cluster (a raiz)
    while len(clusters) > 1:
        a, b, dmin = _closest_pair(matriz, clusters)

        novo_newick, altura = _newick_join(a, b, dmin, clusters)

        # Prepara chaves restantes (atuais na matriz) sem (a, b)
        restantes = [k for k in list(matriz.keys()) if k not in (a, b)]

        # Cria linha/coluna do novo cluster antes de remover a/b
        _average_linkage_row(matriz, a, b, novo_newick, restantes)

        # Remove antigos e atualiza clusters
        _remove_old_clusters(matriz, a, b)
        clusters.pop(a, None)
        clusters.pop(b, None)
        clusters[novo_newick] = (novo_newick, altura)

    # Único elemento remanescente é a raiz em Newick
    return f"{next(iter(clusters))};"
