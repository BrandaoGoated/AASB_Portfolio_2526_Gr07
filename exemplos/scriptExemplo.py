from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]  # .../AASB
sys.path.insert(0, str(ROOT))

from bioinf.sequencias import dna, rna, proteina, dna_para_rna, dna_reverso

from bioinf.motifs import (
    iupac_para_regex, procura_iupac,
    prosite_para_regex, procura_prosite,
    digestao_dna,
    pwm, imprime_pwm,
    prob_gerar_sequencia, seq_mais_provavel,
    pssm_de_pwm, melhor_subsequencia
)

from bioinf.alignments import (
    Blosum62,
    needleman_wunsch,
    smith_waterman,
    progressive_alignment,
    consensus,
)

from bioinf.blast import (
    gerar_palavras,
    indexar_banco,
    encontrar_hits,
    estender_hit,
    blast_simplificado
)

from bioinf.filo import calcular_matriz_distancias, upgma


def print_header(title: str) -> None:
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)


def main():
    # =========================================================
    # 1) Sequências Biológicas
    # =========================================================
    print_header("1) Sequências Biológicas")

    print("dna('ACTG') ->", dna("ACTG"))
    print("dna('ATX') ->", dna("ATX"))  # inválida

    print("rna('AUGC') ->", rna("AUGC"))
    print("rna('AUT') ->", rna("AUT"))  # inválida (T não é RNA)

    print("proteina('ACDE') ->", proteina("ACDE"))
    print("proteina('ACDZ') ->", proteina("ACDZ"))  # inválida (Z fora do set do teu código)

    print("dna_para_rna('ACT') ->", dna_para_rna("ACT"))
    print("dna_reverso('ACT') ->", dna_reverso("ACT"))

    # exemplo de exceção
    try:
        print("dna_para_rna('ATX') ->", dna_para_rna("ATX"))
    except ValueError as e:
        print("dna_para_rna('ATX') levantou ValueError ->", e)

    # =========================================================
    # 2) Alinhamento de Sequências (global/local/múltiplo)
    # =========================================================
    print_header("2) Alinhamentos")

    bl = Blosum62()

    # Needleman-Wunsch (global)
    score_nw, a1, a2 = needleman_wunsch("AC", "AG", bl, gap=-1)
    print("Needleman-Wunsch (global)")
    print("Score:", score_nw)
    print(a1)
    print(a2)

    # Smith-Waterman (local)
    score_sw, l1, l2 = smith_waterman("PAWHEAE", "HEAGAWGHEE", bl, gap=-1)
    print("\nSmith-Waterman (local)")
    print("Score:", score_sw)
    print(l1)
    print(l2)

    # Alinhamento progressivo (múltiplo) + consenso
    seqs_msa = ["ATGC", "ATCC", "ATGG"]
    msa = progressive_alignment(seqs_msa, bl, gap=-1)
    print("\nProgressive alignment (MSA)")
    for i, s in enumerate(msa, start=1):
        print(f"MSA {i}:", s)
    print("Consenso:", consensus(msa))

    # =========================================================
    # 3) Motifs e Padrões
    # =========================================================
    print_header("3) Motifs e Padrões")

    # IUPAC
    print("iupac_para_regex('ATGN') ->", iupac_para_regex("ATGN"))
    print("procura_iupac('ATGCCATG', 'ATG') ->", procura_iupac("ATGCCATG", "ATG"))

    # PROSITE (subconjunto)
    prosite = "C-x(2)-C"
    print("\nPROSITE")
    print(f"prosite_para_regex('{prosite}') ->", prosite_para_regex(prosite))
    print("procura_prosite('ACCC', 'C-x(2)-C') ->", procura_prosite("ACCC", prosite))

    # Enzimas de restrição (com '^')
    seq_dig = "GAATTCGAATTC"
    sitio = "G^AATTC"
    cortes, fragmentos = digestao_dna(seq_dig, sitio)
    print("\nDigestão DNA")
    print("Sequência:", seq_dig)
    print("Sítio:", sitio)
    print("Cortes:", cortes)
    print("Fragmentos:", fragmentos)

    # PWM + probabilidade + janela mais provável
    print("\nPWM/PSSM")
    motifs = ["ATG", "AAG", "ATG", "ACG"]  # exemplos (mesmo comprimento)
    m_pwm = pwm(motifs, tipo="DNA", pseudocontagem=1)
    print("PWM (arredondada):", imprime_pwm(m_pwm, casas=3))

    kmer = "ATG"
    print(f"prob_gerar_sequencia('{kmer}', PWM) ->", prob_gerar_sequencia(kmer, m_pwm))

    seq_long = "CCCATGAAATGCCC"
    print("seq_mais_provavel(seq_long, PWM) ->", seq_mais_provavel(seq_long, m_pwm))

    # PSSM + melhor subsequência
    m_pssm = pssm_de_pwm(m_pwm, alfabeto="ACGT")
    best_kmer, best_pos, best_score = melhor_subsequencia(m_pssm, seq_long)
    print("melhor_subsequencia(PSSM, seq_long) ->", (best_kmer, best_pos, best_score))

    # =========================================================
    # 4) BLAST Simplificado
    # =========================================================
    print_header("4) BLAST Simplificado")

    query = "ATGCATGCA"
    banco = ["ATGCATGCA", "TTGCATGGA", "ATGCGTACA"]

    print("Query:", query)
    print("Banco:", banco)

    # Query map / palavras
    print("\nQuery words (k=3):", gerar_palavras(query, k=3))

    # Index + hits
    idx = indexar_banco(banco, k=3)
    hits = encontrar_hits(query, idx, k=3)
    print("Nº de k-mers no índice:", len(idx))
    print("Primeiros hits (até 10):", hits[:10])

    # Extensão de 1 hit (se existir)
    if hits:
        i_query, id_seq, i_seq = hits[0]
        score_ext, q_aln, s_aln = estender_hit(query, banco[id_seq], i_query, i_seq, k=3)
        print("\nExtensão do 1º hit:")
        print("hit =", (i_query, id_seq, i_seq))
        print("score =", score_ext)
        print("q_aln =", q_aln)
        print("s_aln =", s_aln)

    # BLAST final (ordenado por score)
    resultados = blast_simplificado(query, banco, k=3)
    print("\nTop BLAST hits (até 3):")
    for r in resultados[:3]:
        print(r)

    # =========================================================
    # 5) Análise Filogenética
    # =========================================================
    print_header("5) Filogenia")

    seqs = ["ATGC", "ATGA", "TTGC"]
    matriz = calcular_matriz_distancias(seqs)

    print("Sequências:", seqs)
    print("Matriz de distâncias:")
    for a in seqs:
        row = {b: round(matriz[a][b], 3) for b in seqs}
        print(a, "->", row)

    arvore = upgma(matriz)
    print("UPGMA (Newick):", arvore)


if __name__ == "__main__":
    main()
