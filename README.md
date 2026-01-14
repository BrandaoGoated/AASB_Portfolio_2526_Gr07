# AASB_Portfolio_2526_Gr07

# Portefólio de Algoritmos para análise de sequências biológicas (AASB 2025/2026)

Biblioteca Python que implementa os principais algoritmos abordados na UC. O projeto foi desenvolvido com foco em correção algorítmica, qualidade de código, documentação completa (Sphinx) e testes unitários (cobertura ≥ 80%), sendo facilmente importável e reutilizável por terceiros. [file:16]

---

## Conteúdo do Projeto

Este portefólio inclui implementações dos seguintes tópicos: [file:16]

### 1. Sequências Biológicas
- Validação de sequências (DNA, RNA e proteínas). [file:22]
- Transcrição DNA → RNA. [file:22]
- Reverso-complemento de DNA. [file:22]

### 2. Alinhamento de Sequências
- Matriz de substituição BLOSUM62. [file:18]
- Needleman–Wunsch (alinhamento global): score + reconstrução. [file:18]
- Smith–Waterman (alinhamento local): score + reconstrução. [file:18]
- Alinhamento múltiplo progressivo + consenso. [file:18]

### 3. Motifs e Padrões
- Procura de padrões com ambiguidades (IUPAC). [file:23]
- Conversão PROSITE → expressões regulares. [file:23]
- Digestão por enzimas de restrição (regex/cortes/fragmentação). [file:23]
- PWM e PSSM: construção, probabilidade e melhor subsequência. [file:23]

### 4. BLAST Simplificado
- Geração de k-mers, indexação do banco, hits, extensão e melhor alinhamento. [file:25]

### 5. Análise Filogenética
- Construção de matriz de distâncias por p-distância. [file:21]
- UPGMA: construção de árvore em formato Newick. [file:21]

---

# Uso dos códigos

## `sequencias.py`

```python
from bioinf.sequencias import dna, rna, proteina, dna_para_rna, dna_reverso

# Validação
print(dna("ACTG"))                     # "DNA Válido"
print(rna("ACGU"))                     # "RNA Válido"
print(proteina("ACDEFGHIKLMNPQRSTVWY"))# "Proteína Válida"

# Transcrição DNA -> RNA
print(dna_para_rna("ATGC"))            # "AUGC"

# Reverso-complemento
print(dna_reverso("ATGC"))             # "GCAT"

