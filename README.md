# AASB_Portfolio_2526_Gr07
text
# Portefólio de Algoritmos para Análise de Sequências Biológicas (AASB 2025/2026)

Biblioteca Python que implementa os principais algoritmos abordados na UC.
O projeto foi desenvolvido com foco em correção algorítmica, qualidade de código, documentação completa (Sphinx) e testes unitários (cobertura ≥ 80%), sendo facilmente importável e reutilizável por terceiros. [file:16]

---

## Estrutura do projeto

.
│ requirements.txt
├─ bioinf/
│ ├─ init.py
│ ├─ alignments.py
│ ├─ blast.py
│ ├─ filo.py
│ ├─ motifs.py
│ └─ sequencias.py
├─ docs/
│ ├─ Makefile
│ ├─ make.bat
│ ├─ source/
│ └─ build/html/
├─ exemplos/
│ └─ scriptExemplo.py
├─ tests/
│ ├─ test_alignments_unittest.py
│ ├─ test_blast.py
│ ├─ test_filo.py
│ ├─ test_motifs.py
│ └─ test_sequencias.py
├─ radon/
└─ htmlcov/

text

> Nota: as pastas `radon/` e `htmlcov/` contêm, respetivamente, relatórios Radon e relatórios de cobertura (coverage HTML). [file:16]

---

## Conteúdo do Projeto

### 1. Sequências Biológicas (`bioinf/sequencias.py`)
Funções principais (usadas nos testes): `dna`, `rna`, `proteina`, `dna_para_rna`, `dna_reverso`. [file:27]

### 2. Alinhamento de Sequências (`bioinf/alignments.py`)
Inclui: `Blosum62`, `_clean_seq`, `_argmax`, `simple_substitution_matrix`, `needleman_wunsch`, `smith_waterman`, `consensus`, `progressive_alignment`. [file:24]

### 3. Motifs e Padrões (`bioinf/motifs.py`)
Inclui: `iupac_para_regex`, `procura_iupac`, `prosite_para_regex`, `procura_prosite`, `digestao_dna`, `pwm`, `prob_gerar_sequencia`, `seq_mais_provavel`, `pssm_de_pwm`, `score_kmer`, `melhor_subsequencia`. [file:23]

### 4. BLAST Simplificado (`bioinf/blast.py`)
Inclui: `gerar_palavras`, `indexar_banco`, `encontrar_hits`, `estender_hit`, `blast_simplificado`. [file:25]

### 5. Análise Filogenética (`bioinf/filo.py`)
Inclui: `calcular_matriz_distancias`, `upgma` (gera árvore em formato Newick). [file:26]

---

## Instalação

```bash
pip install -r requirements.txt

Uso dos códigos

    Nota: os exemplos abaixo assumem que estás a executar a partir da raiz do projeto (ou que o package bioinf está acessível no PYTHONPATH). [file:16]

bioinf/sequencias.py

python
from bioinf.sequencias import dna, rna, proteina, dna_para_rna, dna_reverso

# Validação
print(dna("ACTG"))          # "DNA Válido"
print(dna("ANT"))           # "DNA Inválido"
print(rna("ACGU"))          # "RNA Válido"
print(proteina("ACDE"))     # "Proteína Válida" (exemplo)

# Transformações
print(dna_para_rna("ACT"))  # "ACU"
print(dna_reverso("ACT"))   # "AGT" (exemplo)

[file:27]
bioinf/alignments.py

python
from bioinf.alignments import (
    Blosum62,
    simple_substitution_matrix,
    needleman_wunsch,
    smith_waterman,
    consensus,
    progressive_alignment,
)

bl = Blosum62()

# Alinhamento global (Needleman–Wunsch)
score_g, a1_g, a2_g = needleman_wunsch("AC", "AG", bl)
print("Global:", score_g, a1_g, a2_g)

# Alinhamento local (Smith–Waterman)
score_l, a1_l, a2_l = smith_waterman("AC", "AG", bl)
print("Local:", score_l, a1_l, a2_l)

# Consenso
aln = ["AC-", "A-G", "AAG"]
print("Consenso:", consensus(aln))

# Alinhamento progressivo
seqs = ["AC", "AG", "AT"]
print("Progressivo:", progressive_alignment(seqs, bl))

[file:24]
bioinf/motifs.py

python
from bioinf.motifs import (
    iupac_para_regex,
    procura_iupac,
    prosite_para_regex,
    procura_prosite,
    digestao_dna,
    pwm,
    prob_gerar_sequencia,
    seq_mais_provavel,
    pssm_de_pwm,
    score_kmer,
    melhor_subsequencia,
)

# IUPAC
print(iupac_para_regex("CCWGG"))       # "CC[AT]GG"
print(procura_iupac("AAAA", "AA"))     #[13]

# PROSITE
print(prosite_para_regex("[AC]-x-V-x(4)-{ED}"))
print(procura_prosite("GGCATGG", "C-A-T"))  #[13]

# Digestão por enzimas de restrição
cortes, frags = digestao_dna("GAATTCC", "G^AATTC")
print("Cortes:", cortes)       #[14]
print("Fragmentos:", frags)    # ["G", "AATTCC"]

# PWM / PSSM
mat = pwm(["ACG", "ACG", "ATG"], tipo="DNA", pseudocontagem=1)
pssm = pssm_de_pwm(mat, alfabeto="ACGT")
best, pos, score = melhor_subsequencia(pssm, "TTTACGAAA")
print(best, pos, score)

print("Score k-mer:", score_kmer(pssm, "ACG"))

[file:23]
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

[file:25]
bioinf/filo.py

python
from bioinf.filo import calcular_matriz_distancias, upgma

seqs = ["ATGC", "ATGA", "TTGC"]
matriz = calcular_matriz_distancias(seqs)
print("Matriz:", matriz)

arvore = upgma(matriz)
print("Árvore (Newick):", arvore)

[file:26]
Testes Unitários

Os testes encontram-se em tests/ e cobrem casos normais, casos limite e exceções. [file:16][file:23][file:24][file:25][file:26][file:27]
Executar testes

bash
pytest

[file:16]
Executar testes com cobertura

bash
pytest --cov=bioinf --cov-report=term-missing

[file:16]

    O relatório HTML de cobertura é gerado em htmlcov/index.html. [file:16]

Documentação (Sphinx)

A documentação é gerada com Sphinx na pasta docs/ e inclui HTML gerado. [file:16]
Gerar documentação HTML

bash
cd docs
make html

[file:16]

Abrir:

    docs/build/html/index.html [file:16]

Qualidade do Código (Radon)

O trabalho exige classificação Radon mínima B (critério eliminatório). [file:16]

bash
radon cc bioinf/ -a -s
radon mi bioinf/ -s

[file:16]
Autores

Grupo 07 —

    Nome 1 (nº)

    Nome 2 (nº)

    Nome 3 (nº)

UC: Algoritmos e Análise de Sistemas Biológicos
Ano letivo: 2025/2026
