import unittest

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


class TestMotifs(unittest.TestCase):

    # -------------------------
    # IUPAC
    # -------------------------
    def test_iupac_para_regex(self):
        self.assertEqual(iupac_para_regex("CCWGG"), "CC[AT]GG")

    def test_iupac_para_regex_invalid(self):
        with self.assertRaises(ValueError):
            iupac_para_regex("CCZGG")

    def test_procura_iupac_sem_overlap(self):
        # "AA" em "AAAA" sem overlap -> [0, 2]
        self.assertEqual(procura_iupac("AAAA", "AA"), [0, 2])

    # -------------------------
    # PROSITE
    # -------------------------
    def test_prosite_conversion(self):
        prosite = "[AC]-x-V-x(4)-{ED}"
        regex = prosite_para_regex(prosite)
        self.assertEqual(regex, "[AC].V.{4}[^ED]")

    def test_procura_prosite(self):
        # "C-A-T" vira "CAT"
        seq = "GGCATGG"
        padrao = "C-A-T"
        pos = procura_prosite(seq, padrao)
        self.assertEqual(pos, [2])

    # -------------------------
    # Digestão (restrição)
    # -------------------------
    def test_digestao_dna_ecori(self):
        cortes, frags = digestao_dna("GAATTCC", "G^AATTC")
        self.assertEqual(cortes, [1])
        self.assertEqual(frags, ["G", "AATTCC"])

    def test_digestao_dna_requires_caret(self):
        with self.assertRaises(ValueError):
            digestao_dna("GAATTCC", "GAATTC")

    def test_digestao_dna_none_raises(self):
        with self.assertRaises(TypeError):
            digestao_dna(None, "G^AATTC")

    # -------------------------
    # PWM
    # -------------------------
    def test_pwm_laplace(self):
        # Seqs ["A","A"], DNA (ACGT), pseudocontagem=1:
        # denom = n + a*pseudo = 2 + 4 = 6
        # A = (2+1)/6 = 3/6 ; C = (0+1)/6 = 1/6
        mat = pwm(["A", "A"], tipo="DNA", pseudocontagem=1)
        self.assertAlmostEqual(mat[0]["A"], 3/6)
        self.assertAlmostEqual(mat[0]["C"], 1/6)

    def test_prob_gerar_sequencia_len_raises(self):
        mat = pwm(["AA", "AT"], tipo="DNA", pseudocontagem=1)
        with self.assertRaises(ValueError):
            prob_gerar_sequencia("A", mat)

    def test_seq_mais_provavel(self):
        mat = pwm(["AAA", "AAA"], tipo="DNA", pseudocontagem=1)
        subs = seq_mais_provavel("TTTAAATTT", mat)
        self.assertEqual(subs, ["AAA"])  # devolve lista

    # -------------------------
    # PSSM
    # -------------------------
    def test_pssm_e_melhor_subsequencia(self):
        mat = pwm(["ACG", "ACG", "ATG"], tipo="DNA", pseudocontagem=1)
        pssm = pssm_de_pwm(mat, alfabeto="ACGT")
        best, pos, score = melhor_subsequencia(pssm, "TTTACGAAA")
        self.assertEqual(best, "ACG")
        self.assertEqual(pos, 3)
        self.assertNotEqual(score, float("-inf"))

    def test_score_kmer_len_raises(self):
        mat = pwm(["ACG", "ACG", "ATG"], tipo="DNA", pseudocontagem=1)
        pssm = pssm_de_pwm(mat, alfabeto="ACGT")
        with self.assertRaises(ValueError):
            score_kmer(pssm, "AC")

    def test_melhor_subsequencia_short_seq_raises(self):
        mat = pwm(["ACG", "ACG", "ATG"], tipo="DNA", pseudocontagem=1)
        pssm = pssm_de_pwm(mat, alfabeto="ACGT")
        with self.assertRaises(ValueError):
            melhor_subsequencia(pssm, "AC")


if __name__ == "__main__":
    unittest.main()
