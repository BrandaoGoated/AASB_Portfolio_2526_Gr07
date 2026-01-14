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
```

## `alignments.py`

```python
from bioinf.alignments import (
    Blosum62,
    simple_substitution_matrix,
    needleman_wunsch,
    smith_waterman,
    consensus,
    progressive_alignment,
)

bl = Blosum62()

# Needleman–Wunsch (global)
score_g, a1_g, a2_g = needleman_wunsch("AC", "AG", bl, gap=-1)
print("Global:", score_g, a1_g, a2_g)

# Smith–Waterman (local)
score_l, a1_l, a2_l = smith_waterman("AC", "AG", bl, gap=-1)
print("Local:", score_l, a1_l, a2_l)

# Consenso de um alinhamento (exemplo)
aln = ["AC-", "A-G", "AAG"]
print("Consenso:", consensus(aln))

# Alinhamento múltiplo progressivo (exemplo)
seqs = ["AC", "AG", "AT"]
aln_prog = progressive_alignment(seqs, bl, gap=-1)
print("Progressivo:", aln_prog)

# Matriz de substituição simples (opcional)
mat = simple_substitution_matrix(match=2, mismatch=-1, alphabet="ACGT")
print(mat[("A", "A")], mat[("A", "C")])
```
## `blast.py`
```python
from bioinf.blast import (
    gerar_palavras,
    indexar_banco,
    encontrar_hits,
    estender_hit,
    blast_simplificado,
)

query = "ATGCATGCA"
banco = ["ATGCATGCA", "TTGCATGGA", "ATGCGTACA"]
k = 3

# k-mers (palavras)
print(gerar_palavras(query, k))

# Indexação do banco
idx = indexar_banco(banco, k)
print("Index keys:", list(idx.keys())[:5])

# Hits
hits = encontrar_hits(query, idx, k)
print("Hits (exemplo):", hits[:5])

# Extensão de um hit (exemplo)
score, aq, as_ = estender_hit(query, banco, i_query=3, i_seq=3, k=k, match=2, mismatch=-1)
print(score)
print(aq)
print(as_)

# BLAST simplificado (ranking de resultados)
resultados = blast_simplificado(query, banco, k=3)
for score, id_seq, q_aln, s_aln in resultados:
    print(f"Seq {id_seq}: score={score}")
    print(q_aln)
    print(s_aln)
```
## filo.py
``` python
from bioinf.filo import calcular_matriz_distancias, upgma

seqs = ["ATGC", "ATGA", "TTGC"]

# Matriz de distâncias (p-distância)
matriz = calcular_matriz_distancias(seqs)
print(matriz)

# UPGMA -> devolve Newick (string terminada em ';')
newick = upgma(matriz)
print(newick)

```
## motifs.py
``` python
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
print(iupac_para_regex("CCWGG"))      # "CC[AT]GG"
print(procura_iupac("AAAA", "AA"))    #  (sem overlap)[1]

# PROSITE
print(prosite_para_regex("[AC]-x-V-x(4)-{ED}"))
print(procura_prosite("GGCATGG", "C-A-T"))  #[1]

# Digestão por enzima de restrição
cortes, frags = digestao_dna("GAATTCC", "G^AATTC")
print(cortes, frags)

# PWM / PSSM
mat = pwm(["ACG", "ACG", "ATG"], tipo="DNA", pseudocontagem=1)
print(prob_gerar_sequencia("ACG", mat))

subs = seq_mais_provavel("TTTAAATTT", pwm(["AAA", "AAA"], tipo="DNA", pseudocontagem=1))
print(subs)

pssm = pssm_de_pwm(mat, alfabeto="ACGT")
best, pos, score = melhor_subsequencia(pssm, "TTTACGAAA")
print(best, pos, score)

print(score_kmer(pssm, "ACG"))

```
## Testes Unitários
Os testes foram desenvolvidos com unittest / pytest e cobrem casos normais, casos limite e exceções.
```
pytest
```
Executar testes com cobertura
```
pytest --cov=bioinf --cov-report=term-missing
```

## Documentação (Sphinx)
Gerar documentação:
```
cd docs
make html
```
Abrir:
docs/build/html/index.html

## Qualidade do Código (Radon)
```
radon cc bioinf/ -a -s
radon mi bioinf/ -s
```
## Autores
Grupo 07 —
    David Brandão - PG59418
    Nome 2 (nº)
    Nome 3 (nº)

UC: Algoritmos e Análise de Sistemas Biológicos
Ano letivo: 2025/2026
