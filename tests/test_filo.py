# tests/test_filo_unittest.py

import unittest
from bioinf.filo import calcular_matriz_distancias, upgma

class TestFilogenia(unittest.TestCase):
    # ---------- calcular_matriz_distancias ----------
    def test_matriz_distancias_basico(self):
        sequencias = ["ATGC", "ATGA", "TTGC"]
        matriz = calcular_matriz_distancias(sequencias)
        # Verifica tamanho da matriz
        self.assertEqual(len(matriz), 3)
        # Distância entre seq1 e seq2 (diferença em 1 posição)
        self.assertEqual(matriz["ATGC"]["ATGA"], 0.25)
        # Distância entre seq1 e seq3 (diferença em 1 posição)
        self.assertEqual(matriz["ATGC"]["TTGC"], 0.25)
        # Distância entre seq2 e seq3 (diferença em 2 posições)
        self.assertEqual(matriz["ATGA"]["TTGC"], 0.5)

    def test_matriz_distancias_tamanhos_diferentes(self):
        sequencias = ["ATGC", "ATG"]
        with self.assertRaises(ValueError):
            calcular_matriz_distancias(sequencias)

    def test_matriz_distancias_uma_sequencia(self):
        sequencias = ["ATGC"]
        matriz = calcular_matriz_distancias(sequencias)
        self.assertEqual(matriz["ATGC"]["ATGC"], 0.0)

    # ---------- upgma ----------
    def test_upgma_basico(self):
        sequencias = ["ATGC", "ATGA", "TTGC"]
        matriz = calcular_matriz_distancias(sequencias)
        arvore = upgma(matriz)
        # Deve terminar com ponto e vírgula (formato Newick)
        self.assertTrue(arvore.endswith(";"))
        # Deve conter os nomes das sequências
        for seq in sequencias:
            self.assertIn(seq, arvore)
        # Deve conter parênteses indicando agrupamento
        self.assertIn("(", arvore)
        self.assertIn(")", arvore)

    def test_upgma_duas_sequencias(self):
        sequencias = ["ATGC", "ATGA"]
        matriz = calcular_matriz_distancias(sequencias)
        arvore = upgma(matriz)
        # Formato esperado: (ATGC:x,ATGA:y);
        self.assertTrue(arvore.startswith("("))
        self.assertTrue(arvore.endswith(");"))

    def test_upgma_uma_sequencia(self):
        sequencias = ["ATGC"]
        matriz = calcular_matriz_distancias(sequencias)
        arvore = upgma(matriz)
        # Deve apenas devolver a sequência com ponto e vírgula
        self.assertEqual(arvore, "ATGC;")

if __name__ == "__main__":
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
