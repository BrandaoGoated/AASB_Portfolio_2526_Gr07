import unittest

from bioinf.sequencias import dna, rna, proteina, dna_para_rna, dna_reverso


class TestSequencias(unittest.TestCase):

    # -------------------------
    # Validação
    # -------------------------
    def test_dna_valido(self):
        self.assertEqual(dna("ACTG"), "DNA Válido")

    def test_dna_invalido_caracter(self):
        self.assertEqual(dna("ANT"), "DNA Inválido")

    def test_dna_invalido_vazio(self):
        self.assertEqual(dna(""), "DNA Inválido")

    def test_rna_valido(self):
        self.assertEqual(rna("ACGU"), "RNA Válido")

    def test_rna_invalido_vazio(self):
        self.assertEqual(rna(""), "RNA Inválido")

    def test_proteina_valida(self):
        self.assertEqual(proteina("ACDEFGHIKLMNPQRSTVWY"), "Proteína Válida")

    def test_proteina_invalida_caracter(self):
        self.assertEqual(proteina("ACD*"), "Proteína Inválida")

    def test_case_insensitive_e_espacos(self):
        self.assertEqual(dna("  actg "), "DNA Válido")
        self.assertEqual(rna(" acgu "), "RNA Válido")

    # -------------------------
    # Transformações
    # -------------------------
    def test_dna_para_rna_ok(self):
        self.assertEqual(dna_para_rna("ACT"), "ACU")

    def test_dna_para_rna_raises(self):
        # tem U, logo não é DNA
        with self.assertRaises(ValueError):
            dna_para_rna("ACU")

    def test_dna_reverso_ok(self):
        self.assertEqual(dna_reverso("ACT"), "AGT")

    def test_dna_reverso_raises(self):
        with self.assertRaises(ValueError):
            dna_reverso("ACU")


if __name__ == "__main__":
    unittest.main()
