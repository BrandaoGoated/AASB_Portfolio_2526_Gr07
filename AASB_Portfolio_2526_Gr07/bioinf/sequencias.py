def dna(seq):
    """
    Verifica se uma sequência é DNA válido (apenas A, C, G, T).

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

    """
    seq = seq.upper().strip()
    if dna(seq) != "DNA Válido":
        raise ValueError("DNA Inválido")
    return seq.replace("T", "U")


def dna_reverso(seq):
    """
    Calcula o reverso complementar de uma sequência de DNA.

    """
    seq = seq.upper().strip()
    if dna(seq) != "DNA Válido":
        raise ValueError("DNA Inválido")
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp[c] for c in seq[::-1])
