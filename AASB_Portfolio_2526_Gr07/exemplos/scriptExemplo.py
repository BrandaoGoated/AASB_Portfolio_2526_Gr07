from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]   # .../AASB
sys.path.insert(0, str(ROOT))
from bioinf.sequencias import dna, dna_para_rna, dna_reverso
from bioinf.motifs import iupac_para_regex, procura_iupac
from bioinf.alignments import Blosum62, needleman_wunsch
from bioinf.blast import blast_simplificado
from bioinf.filo import calcular_matriz_distancias, upgma


def main():
    # Sequências
    print("DNA('ACTG') ->", dna("ACTG"))
    print("dna_para_rna('ACT') ->", dna_para_rna("ACT"))
    print("dna_reverso('ACT') ->", dna_reverso("ACT"))

    # Motifs
    print("iupac_para_regex('ATGN') ->", iupac_para_regex("ATGN"))
    print("procura_iupac('ATGCCATG', 'ATG') ->", procura_iupac("ATGCCATG", "ATG"))

    # Alinhamento
    bl = Blosum62()
    score, a1, a2 = needleman_wunsch("AC", "AG", bl)
    print("NW score:", score)
    print(a1)
    print(a2)

    # BLAST 
    query = "ATGCATGCA"
    banco = ["ATGCATGCA", "TTGCATGGA", "ATGCGTACA"]
    resultados = blast_simplificado(query, banco, k=3)
    print("Top BLAST hit:", resultados[0] if resultados else None)

    # Filogenia
    seqs = ["ATGC", "ATGA", "TTGC"]
    matriz = calcular_matriz_distancias(seqs)
    arvore = upgma(matriz)
    print("UPGMA Newick:", arvore)


if __name__ == "__main__":
    main()
