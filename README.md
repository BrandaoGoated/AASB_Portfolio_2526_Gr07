# AASB_Portfolio_2526_Gr07
# Portefólio de Algoritmos para Análise de Sequências Biológicas (AASB 2025/2026)

Biblioteca Python que implementa os principais algoritmos abordados na UC.
O projeto foi desenvolvido com foco em correção algorítmica, qualidade de código, documentação completa (Sphinx) e testes unitários (cobertura ≥ 80%), sendo facilmente importável e reutilizável por terceiros.

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

> Nota: as pastas `radon/` e `htmlcov/` contêm, respetivamente, relatórios Radon e relatórios de cobertura (coverage HTML).

---

## Conteúdo do Projeto

### 1. Sequências Biológicas (`bioinf/sequencias.py`)
Funções principais: `dna`, `rna`, `proteina`, `dna_para_rna`, `dna_reverso`.

### 2. Alinhamento de Sequências (`bioinf/alignments.py`)
Inclui `Blosum62`, `simple_substitution_matrix`, `needleman_wunsch`, `smith_waterman`, `consensus`, `progressive_alignment`.

### 3. Motifs e Padrões (`bioinf/motifs.py`)
Inclui `iupac_para_regex`, `procura_iupac`, `prosite_para_regex`, `procura_prosite`, `digestao_dna`, `pwm`, `prob_gerar_sequencia`, `seq_mais_provavel`, `pssm_de_pwm`, `score_kmer`, `melhor_subsequencia`.

### 4. BLAST Simplificado (`bioinf/blast.py`)
Inclui `gerar_palavras`, `indexar_banco`, `encontrar_hits`, `estender_hit`, `blast_simplificado`.

### 5. Análise Filogenética (`bioinf/filo.py`)
Inclui `calcular_matriz_distancias` e `upgma`.

---

## Instalação

```bash
pip install -r requirements.txt

Uso dos códigos
bioinf/sequencias.py

python
from bioinf.sequencias import dna, rna, proteina, dna_para_rna, dna_reverso

print(dna("ACTG"))
print(rna("ACGU"))
print(proteina("ACDEFGHIKLMNPQRSTVWY"))

print(dna_para_rna("ACT"))
print(dna_reverso("ACT"))

bioinf/alignments.py

python
from bioinf.alignments import Blosum62, needleman_wunsch, smith_waterman, consensus, progressive_alignment

bl = Blosum62()

score_g, a1_g, a2_g = needleman_wunsch("AC", "AG", bl)
print("Global:", score_g, a1_g, a2_g)

score_l, a1_l, a2_l = smith_waterman("AC", "AG", bl)
print("Local:", score_l, a1_l, a2_l)

aln = ["AC-", "A-G", "AAG"]
print("Consenso:", consensus(aln))

seqs = ["AC", "AG", "AT"]
print("Progressivo:", progressive_alignment(seqs, bl))

bioinf/motifs.py

python
from bioinf.motifs import (
    iupac_para_regex, procura_iupac,
    prosite_para_regex, procura_prosite,
    digestao_dna,
    pwm, prob_gerar_sequencia, seq_mais_provavel,
    pssm_de_pwm, score_kmer, melhor_subsequencia
)

print(iupac_para_regex("CCWGG"))
print(procura_iupac("AAAA", "AA"))

print(prosite_para_regex("[AC]-x-V-x(4)-{ED}"))
print(procura_prosite("GGCATGG", "C-A-T"))

cortes, frags = digestao_dna("GAATTCC", "G^AATTC")
print(cortes, frags)

mat = pwm(["ACG", "ACG", "ATG"], tipo="DNA", pseudocontagem=1)
pssm = pssm_de_pwm(mat, alfabeto="ACGT")
best, pos, score = melhor_subsequencia(pssm, "TTTACGAAA")
print(best, pos, score)

print(score_kmer(pssm, "ACG"))

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

bioinf/filo.py

python
from bioinf.filo import calcular_matriz_distancias, upgma

seqs = ["ATGC", "ATGA", "TTGC"]
matriz = calcular_matriz_distancias(seqs)
print(matriz)

print(upgma(matriz))

Testes Unitários
Executar testes

bash
pytest

Executar testes com cobertura

bash
pytest --cov=bioinf --cov-report=term-missing

    O relatório HTML de cobertura é gerado em htmlcov/index.html.

Documentação (Sphinx)
Gerar documentação HTML

bash
cd docs
make html

Abrir:

    docs/build/html/index.html

Qualidade do Código (Radon)

bash
radon cc bioinf/ -a -s
radon mi bioinf/ -s

Autores

Grupo 07 —
David Brandão PG 59418 - Contruibuição: Código do sequencias.py e de motifs.py e o relacionado aos mesmos.

UC: Algoritmos e Análise de Sistemas Biológicos
Ano letivo: 2025/2026
