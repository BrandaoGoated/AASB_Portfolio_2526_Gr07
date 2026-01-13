def dna(seq):
    """
    Verifica se uma sequência é DNA válido (apenas A, C, G, T).

A validação é feita após normalização (upper + strip). Se a sequência
for vazia, é considerada inválida.

Args:
    seq (str): Sequência a validar.

Returns:
    str: "DNA Válido" se todos os caracteres forem A/C/G/T, caso contrário
    "DNA Inválido".

Example:
    >>> dna("ATGC")
    'DNA Válido'
    >>> dna("ATX")
    'DNA Inválido'

    """
    seq = seq.upper().strip()
    if not seq:
        return "DNA Inválido"
    if all(c in "ACGT" for c in seq):
        return "DNA Válido"
    return "DNA Inválido"


def rna(seq):
    """
Verifica se uma sequência é RNA válido (apenas A, C, G, U).

A validação é feita após normalização (upper + strip). Se a sequência
for vazia, é considerada inválida.

Args:
    seq (str): Sequência a validar.

Returns:
    str: "RNA Válido" se todos os caracteres forem A/C/G/U, caso contrário
    "RNA Inválido".

Example:
    >>> rna("AUGC")
    'RNA Válido'
    >>> rna("")
    'RNA Inválido'
""" 

    seq = seq.upper().strip()
    if not seq:
        return "RNA Inválido"
    if all(c in "ACGU" for c in seq):
        return "RNA Válido"
    return "RNA Inválido"


def proteina(seq):
    """
Verifica se uma sequência é proteína válida (20 aminoácidos padrão).

Considera apenas os 20 aminoácidos padrão: ACDEFGHIKLMNPQRSTVWY.
Se a sequência for vazia, é considerada inválida.

Args:
    seq (str): Sequência a validar.

Returns:
    str: "Proteína Válida" se todos os caracteres pertencerem ao conjunto
    permitido, caso contrário "Proteína Inválida".

Example:
    >>> proteina("ACDE")
    'Proteína Válida'
    >>> proteina("ACDZ")
    'Proteína Inválida'
"""

    seq = seq.upper().strip()
    if not seq:
        return "Proteína Inválida"
    aa = "ACDEFGHIKLMNPQRSTVWY"
    if all(c in aa for c in seq):
        return "Proteína Válida"
    return "Proteína Inválida"


def dna_para_rna(seq):
    """
Transcreve DNA para RNA (T -> U).

A função valida primeiro se a sequência é DNA; se não for, levanta
uma exceção.

Args:
    seq (str): Sequência de DNA.

Returns:
    str: Sequência transcrita para RNA.

Raises:
    ValueError: Se a sequência não for DNA válido.

Example:
    >>> dna_para_rna("ATGC")
    'AUGC'
"""
    seq = seq.upper().strip()
    if dna(seq) != "DNA Válido":
        raise ValueError("DNA Inválido")
    return seq.replace("T", "U")


def dna_reverso(seq):
    """
Calcula o reverso complementar de uma sequência de DNA.

Usa o emparelhamento A<->T e C<->G e devolve a sequência complementar
no sentido inverso (reverse complement).

Args:
    seq (str): Sequência de DNA.

Returns:
    str: Reverso complementar.

Raises:
    ValueError: Se a sequência não for DNA válido.

Example:
    >>> dna_reverso("ATGC")
    'GCAT'
"""
    seq = seq.upper().strip()
    if dna(seq) != "DNA Válido":
        raise ValueError("DNA Inválido")
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp[c] for c in seq[::-1])
