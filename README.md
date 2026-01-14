# AASB_Portfolio_2526_Gr07
Aqui vai uma versão corrigida (sem “ficar tudo em bash”) e com a tree simplificada para não aparecer nessa linha gigante. O problema acontece quando a fence de Markdown não fecha bem; por isso deixo tudo com blocos curtos e bem fechados (cada bloco começa e termina em ```). [file:17]

Copia e cola isto tudo para o teu README.md:

text
# Portefólio de Algoritmos para Análise de Sequências Biológicas (AASB 2025/2026)

Biblioteca Python que implementa os principais algoritmos abordados na UC.
O projeto foi desenvolvido com foco em correção algorítmica, qualidade de código, documentação completa (Sphinx) e testes unitários (cobertura ≥ 80%), sendo facilmente importável e reutilizável por terceiros.[3]

---

## Estrutura do projeto (simplificada)

.
├── requirements.txt
├── bioinf/
│ ├── init.py
│ ├── alignments.py
│ ├── blast.py
│ ├── filo.py
│ ├── motifs.py
│ └── sequencias.py
├── tests/
│ ├── test_alignments_unittest.py
│ ├── test_blast.py
│ ├── test_filo.py
│ ├── test_motifs.py
│ └── test_sequencias.py
├── docs/ # Sphinx (source + build/html)
├── exemplos/
├── radon/ # relatórios Radon (json)
└── htmlcov/ # relatório coverage (HTML)

text

---

## Instalação

```bash
pip install -r requirements.txt

Uso dos códigos

    Nota: exemplos assumem execução a partir da raiz do projeto (package bioinf acessível).

    ​

bioinf/sequencias.py

python
from bioinf.sequencias import dna, rna, proteina, dna_para_rna, dna_reverso

print(dna("ACTG"))          # "DNA Válido"
print(dna("ANT"))           # "DNA Inválido"
print(rna("ACGU"))          # "RNA Válido"
print(proteina("ACDE"))     # "Proteína Válida" (exemplo)

print(dna_para_rna("ACT"))  # "ACU"
print(dna_reverso("ACT"))   # "AGT" (exemplo)

​
bioinf/alignments.py

python
from bioinf.alignments import (
    Blosum62,
    needleman_wunsch,
    smith_waterman,
    consensus,
    progressive_alignment,
)

bl = Blosum62()

score_g, a1_g, a2_g = needleman_wunsch("AC", "AG", bl)
print("Global:", score_g, a1_g, a2_g)

score_l, a1_l, a2_l = smith_waterman("AC", "AG", bl)
print("Local:", score_l, a1_l, a2_l)

aln = ["AC-", "A-G", "AAG"]
print("Consenso:", consensus(aln))

seqs = ["AC", "AG", "AT"]
print("Progressivo:", progressive_alignment(seqs, bl))

​
bioinf/motifs.py

python
from bioinf.motifs import (
    iupac_para_regex,
    procura_iupac,
    prosite_para_regex,
    procura_prosite,
    digestao_dna,
    pwm,
    pssm_de_pwm,
    score_kmer,
    melhor_subsequencia,
)

print(iupac_para_regex("CCWGG"))          # "CC[AT]GG"
print(procura_iupac("AAAA", "AA"))        # [0, 2]

print(prosite_para_regex("[AC]-x-V-x(4)-{ED}"))
print(procura_prosite("GGCATGG", "C-A-T"))  # [1]

cortes, frags = digestao_dna("GAATTCC", "G^AATTC")
print("Cortes:", cortes)       # [2]
print("Fragmentos:", frags)    # ["G", "AATTCC"]

mat = pwm(["ACG", "ACG", "ATG"], tipo="DNA", pseudocontagem=1)
pssm = pssm_de_pwm(mat, alfabeto="ACGT")

best, pos, score = melhor_subsequencia(pssm, "TTTACGAAA")
print(best, pos, score)

print("Score k-mer:", score_kmer(pssm, "ACG"))

​
bioinf/blast.py

python
from bioinf.blast import blast_simplificado

query = "ATGCATGCA"
banco = ["ATGCATGCA", "TTGCATGGA", "ATGCGTACA"]

resultados = blast_simplificado(query, banco, k=3)

for score, id_seq, q_aln, s_aln in resultados:
    print(f"Seq {id_seq}: score={score}")
    print(q_aln)
    print(s_aln)

​
bioinf/filo.py

python
from bioinf.filo import calcular_matriz_distancias, upgma

seqs = ["ATGC", "ATGA", "TTGC"]
matriz = calcular_matriz_distancias(seqs)
print("Matriz:", matriz)

arvore = upgma(matriz)
print("Árvore (Newick):", arvore)

​
Testes Unitários

Os testes encontram-se em tests/ e cobrem casos normais, casos limite e exceções.

​

Executar testes:

bash
pytest

Executar testes com cobertura:

bash
pytest --cov=bioinf --cov-report=term-missing

    O relatório HTML de cobertura é gerado em htmlcov/index.html.

    ​

Documentação (Sphinx)

A documentação é gerada com Sphinx e está disponível em HTML.

​

Gerar documentação:

bash
cd docs
make html

Abrir:

    docs/build/html/index.html

Qualidade do Código (Radon)

O trabalho exige classificação Radon mínima B (critério eliminatório).

​

bash
radon cc bioinf/ -a -s
radon mi bioinf/ -s

Autores

Grupo 07 —

    Nome 1 (nº)

    Nome 2 (nº)

    Nome 3 (nº)

UC: Algoritmos e Análise de Sistemas Biológicos
Ano letivo: 2025/2026
