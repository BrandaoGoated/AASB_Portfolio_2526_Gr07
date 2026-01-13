import unittest
from bioinf.alignments import (
    Blosum62, _clean_seq, _argmax, simple_substitution_matrix,
    needleman_wunsch, smith_waterman, consensus, progressive_alignment
)

# -------------------- Classe de testes --------------------

class TestAlignments(unittest.TestCase):

    # -------------------- Utilitários --------------------
    def test_clean_seq_normal(self):
        self.assertEqual(_clean_seq("acgt"), "ACGT")
        self.assertEqual(_clean_seq(" aCg "), "ACG")

    def test_clean_seq_none(self):
        with self.assertRaises(TypeError):
            _clean_seq(None)

    def test_argmax_basico(self):
        self.assertEqual(_argmax([1, 3, 2]), 1)
        self.assertEqual(_argmax([5, 5, 2]), 0)

    # -------------------- Blosum62 --------------------
    def test_blosum62_subst(self):
        bl = Blosum62()
        self.assertEqual(bl.subst('A','A'), 4)
        self.assertEqual(bl.subst('A','G'), 0)
        self.assertEqual(bl.subst('W','Y'), 2)

    # -------------------- Matrizes simples --------------------
    def test_simple_substitution_matrix(self):
        mat = simple_substitution_matrix(match=2, mismatch=-1)
        self.assertEqual(mat[('A','A')], 2)
        self.assertEqual(mat[('A','C')], -1)
        mat2 = simple_substitution_matrix(match=1, mismatch=0, alphabet="AB")
        self.assertEqual(mat2[('A','B')], 0)

    # -------------------- Needleman-Wunsch --------------------
    def test_needleman_wunsch_basico(self):
        bl = Blosum62()
        score, a1, a2 = needleman_wunsch("AC", "AG", bl)
        self.assertIsNotNone(score)
        self.assertEqual(len(a1), len(a2))
        # Testa alinhamento com sequências vazias
        score2, a3, a4 = needleman_wunsch("", "", bl)
        self.assertEqual(score2, 0)
        self.assertEqual(a3, "")
        self.assertEqual(a4, "")

    # -------------------- Smith-Waterman --------------------
    def test_smith_waterman_basico(self):
        bl = Blosum62()
        score, a1, a2 = smith_waterman("AC", "AG", bl)
        self.assertGreaterEqual(score, 0)
        self.assertEqual(len(a1), len(a2))
        score2, a3, a4 = smith_waterman("AAA", "CCC", bl)
        self.assertGreaterEqual(score2, 0)

    # -------------------- Consenso --------------------
    def test_consensus_basico(self):
        aln = ["AC-", "A-G", "AAG"]
        cons = consensus(aln)
        self.assertEqual(cons, "ACG")

    # -------------------- Alinhamento progressivo --------------------
    def test_progressive_alignment_basico(self):
        bl = Blosum62()
        seqs = ["AC", "AG", "AT"]
        aln = progressive_alignment(seqs, bl)
        self.assertEqual(len(aln), 3)
        for seq in aln:
            self.assertEqual(len(seq), len(aln[0]))
        # Testa erro com menos de duas sequências
        with self.assertRaises(ValueError):
            progressive_alignment(["A"], bl)

# -------------------- Executar os testes --------------------
if __name__ == "__main__":
    unittest.main()