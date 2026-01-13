"""
Algoritmos de Alinhamento de Sequências

Inclui:
- Matrizes de substituição simples (PAM/BLOSUM-like reduzidas)
- Needleman-Wunsch (alinhamento global): score + reconstrução
- Smith-Waterman (alinhamento local): score + reconstrução
- Alinhamento progressivo (múltiplo): consenso + alinhamento
"""

import math
from collections import defaultdict

# -------------------- BLOSUM62 --------------------
class Blosum62:
    """
    Matriz de substituição BLOSUM62 (acesso via método `subst`).
    
    A tabela é carregada internamente e disponibiliza scores inteiros para pares
    de aminoácidos.
    
    Example:
        >>> b = Blosum62()
        >>> b.subst("A", "A")
        4
    """
    def __init__(self):
        tabela = """
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  -
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
- -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
"""
        tabela = [l.split() for l in tabela.splitlines() if l]
        headers, *rows = tabela
        self.tab = defaultdict(dict)
        for row in rows:
            aa, *vals = row
            for h, v in zip(headers, vals):
                self.tab[aa][h] = int(v)

    def subst(self, x, y):
        """
        Devolve o score de substituição entre dois símbolos.
        
        Args:
            x (str): Símbolo 1 (ex.: aminoácido).
            y (str): Símbolo 2 (ex.: aminoácido).
        
        Returns:
            int: Score BLOSUM62 para o par (x, y).
        
        Example:
            >>> Blosum62().subst("W", "W")
            11
        """
        return self.tab[x][y]


# -------------------- Utilitários --------------------
def _clean_seq(seq):
    if seq is None:
        raise TypeError("seq cannot be None")
    return seq.upper().strip()

def _argmax(items):
    return max(range(len(items)), key=lambda i: items[i])

# -------------------- Matrizes simples --------------------
def simple_substitution_matrix(match=1, mismatch=-1, alphabet="ACGT"):
    """
    Cria uma matriz de substituição simples (match/mismatch) para um alfabeto.
    
    Nota: esta função devolve um dicionário indexado por tuplos (a, b).
    
    Args:
        match (int): Score para símbolos iguais. Default: 1.
        mismatch (int): Score para símbolos diferentes. Default: -1.
        alphabet (str): Alfabeto permitido. Default: "ACGT".
    
    Returns:
        dict: Dicionário {(a, b): score}.
    
    Example:
        >>> m = simple_substitution_matrix(match=2, mismatch=-1, alphabet="AT")
        >>> m[("A", "T")]
        -1
    """
    mat = {}
    for a in alphabet:
        for b in alphabet:
            mat[(a,b)] = match if a==b else mismatch
    return mat

# -------------------- Needleman-Wunsch --------------------
def _nw_traceback(bt, s1, s2):
    a1, a2 = [], []
    i, j = len(s1), len(s2)
    while i > 0 or j > 0:
        move = bt[i][j]
        if move == 'D':
            a1.append(s1[i-1]); a2.append(s2[j-1])
            i -= 1; j -= 1
        elif move == 'U':
            a1.append(s1[i-1]); a2.append('-')
            i -= 1
        else:
            a1.append('-'); a2.append(s2[j-1])
            j -= 1
    return ''.join(reversed(a1)), ''.join(reversed(a2))

def needleman_wunsch(seq1, seq2, submat, gap=-1):
    """
    Alinhamento global usando Needleman-Wunsch.
    
    Implementa programação dinâmica e faz traceback para reconstruir as sequências
    alinhadas.
    
    Args:
        seq1 (str): Primeira sequência.
        seq2 (str): Segunda sequência.
        submat: Objeto com método `subst(x, y)` que devolve o score de substituição.
        gap (int): Penalização de gap. Default: -1.
    
    Returns:
        Tuple[int, str, str]: (score, seq1_alinhada, seq2_alinhada).
    
    Raises:
        TypeError: Se alguma sequência for None (propagado de `_clean_seq`).
    
    Example:
        >>> b = Blosum62()
        >>> score, a1, a2 = needleman_wunsch("A", "A", b, gap=-1)
        >>> (score, a1, a2)
        (4, 'A', 'A')
    """
    s1, s2 = _clean_seq(seq1), _clean_seq(seq2)
    n, m = len(s1), len(s2)
    dp = [[0]*(m+1) for _ in range(n+1)]
    bt = [[None]*(m+1) for _ in range(n+1)]
    for i in range(1,n+1): dp[i][0]=i*gap; bt[i][0]='U'
    for j in range(1,m+1): dp[0][j]=j*gap; bt[0][j]='L'
    for i in range(1,n+1):
        for j in range(1,m+1):
            scores=[dp[i-1][j-1]+submat.subst(s1[i-1],s2[j-1]),
                    dp[i-1][j]+gap,
                    dp[i][j-1]+gap]
            k=_argmax(scores)
            dp[i][j]=scores[k]
            bt[i][j]=['D','U','L'][k]
    return dp[n][m], *_nw_traceback(bt,s1,s2)

# -------------------- Smith-Waterman --------------------
def _sw_traceback(dp, bt, s1, s2, start_i, start_j):
    a1,a2=[],[]
    i,j=start_i,start_j
    while i>0 and j>0 and dp[i][j]>0:
        move=bt[i][j]
        if move=='D': a1.append(s1[i-1]); a2.append(s2[j-1]); i-=1;j-=1
        elif move=='U': a1.append(s1[i-1]); a2.append('-'); i-=1
        elif move=='L': a1.append('-'); a2.append(s2[j-1]); j-=1
        else: break
    return ''.join(reversed(a1)), ''.join(reversed(a2))

def smith_waterman(seq1, seq2, submat, gap=-1):
    """
    Alinhamento local usando Smith-Waterman.
    
    Usa programação dinâmica com piso a 0 para obter o melhor alinhamento local
    e faz traceback a partir da melhor célula.
    
    Args:
        seq1 (str): Primeira sequência.
        seq2 (str): Segunda sequência.
        submat: Objeto com método `subst(x, y)` que devolve o score de substituição.
        gap (int): Penalização de gap. Default: -1.
    
    Returns:
        Tuple[int, str, str]: (melhor_score_local, seq1_alinhada, seq2_alinhada).
    
    Raises:
        TypeError: Se alguma sequência for None (propagado de `_clean_seq`).
    
    Example:
        >>> b = Blosum62()
        >>> smith_waterman("PAWHEAE", "HEAGAWGHEE", b, gap=-1)[0] > 0
        True
    """
    s1,s2=_clean_seq(seq1),_clean_seq(seq2)
    n,m=len(s1),len(s2)
    dp=[[0]*(m+1) for _ in range(n+1)]
    bt=[[None]*(m+1) for _ in range(n+1)]
    best=(0,0,0)
    for i in range(1,n+1):
        for j in range(1,m+1):
            scores=[0,
                    dp[i-1][j-1]+submat.subst(s1[i-1],s2[j-1]),
                    dp[i-1][j]+gap,
                    dp[i][j-1]+gap]
            k=_argmax(scores)
            dp[i][j]=scores[k]
            bt[i][j]=['0','D','U','L'][k]
            if dp[i][j]>best[0]: best=(dp[i][j],i,j)
    score,i,j=best
    return score,*_sw_traceback(dp,bt,s1,s2,i,j)

# -------------------- Consenso e Alinhamento Progressivo --------------------
def consensus(alignment):
    """
    Calcula a sequência consenso de um alinhamento múltiplo.
    
    Para cada coluna, escolhe o símbolo mais frequente ignorando gaps ('-').
    Se a coluna tiver apenas gaps, o consenso coloca '-'.
    
    Args:
        alignment (List[str]): Lista de sequências alinhadas (mesmo comprimento).
    
    Returns:
        str: Sequência consenso.
    
    Example:
        >>> consensus(["A-C", "ACC", "ATC"])
        'ACC'
    """
    cols = zip(*alignment)
    cons=[]
    for col in cols:
        counts={}
        for c in col:
            if c!='-': counts[c]=counts.get(c,0)+1
        cons.append(max(counts,key=counts.get) if counts else '-')
    return ''.join(cons)

def progressive_alignment(seqs, submat, gap=-1):
    """
    Alinhamento múltiplo progressivo por consenso + alinhamento global.
    
    Estratégia:
    1) Alinha as duas primeiras sequências globalmente.
    2) Calcula consenso do alinhamento atual.
    3) Alinha o consenso (sem '-') com a próxima sequência.
    4) Insere gaps no alinhamento anterior para acompanhar o alinhamento do consenso.
    
    Args:
        seqs (List[str]): Lista de sequências (>= 2).
        submat: Objeto com método `subst(x, y)`.
        gap (int): Penalização de gap. Default: -1.
    
    Returns:
        List[str]: Lista de sequências alinhadas.
    
    Raises:
        ValueError: Se forem fornecidas menos de 2 sequências.
    
    Example:
        >>> b = Blosum62()
        >>> msa = progressive_alignment(["ATGC", "ATCC", "ATGG"], b, gap=-1)
        >>> len(msa)
        3
    """
    if len(seqs)<2: raise ValueError("Need at least two sequences")
    _,a1,a2=needleman_wunsch(seqs[0],seqs[1],submat,gap)
    alignment=[a1,a2]
    for s in seqs[2:]:
        cons=consensus(alignment)
        _,c_aln,s_aln=needleman_wunsch(cons.replace('-',''),s,submat,gap)
        idx=0
        for c in c_aln:
            if c=='-':
                for i in range(len(alignment)):
                    alignment[i]=alignment[i][:idx]+'-'+alignment[i][idx:]
            idx+=1
        alignment.append(s_aln)
    return alignment
