
# tests/test_blast_unittest.py
import unittest
from bioinf.blast import (
    gerar_palavras,
    indexar_banco,
    encontrar_hits,
    estender_hit,
    blast_simplificado,
)

class TestBLASTLike(unittest.TestCase):
    # ---------- gerar_palavras ----------
    def test_gerar_palavras_basico(self):
        consulta = "ATAT"
        k = 2
        esperado = ["AT", "TA", "AT"]
        self.assertEqual(gerar_palavras(consulta, k), esperado)

    def test_gerar_palavras_k_igual_tamanho(self):
        consulta = "ATGC"
        k = 4
        self.assertEqual(gerar_palavras(consulta, k), ["ATGC"])

    def test_gerar_palavras_sem_palavras(self):
        consulta = "AT"
        k = 3
        # O teu código devolve lista vazia quando k > len(seq)
        self.assertEqual(gerar_palavras(consulta, k), [])

    # ---------- indexar_banco ----------
    def test_indexar_banco_basico(self):
        banco = ["ATGC", "TTGC"]
        k = 2
        idx = indexar_banco(banco, k)
        # Chaves esperadas
        self.assertIn("AT", idx)
        self.assertIn("TG", idx)
        self.assertIn("GC", idx)
        # Posições em cada sequência
        self.assertIn((0, 0), idx["AT"])  # "AT" na seq 0, pos 0
        self.assertIn((0, 1), idx["TG"])  # "TG" na seq 0, pos 1
        self.assertIn((0, 2), idx["GC"])  # "GC" na seq 0, pos 2
        self.assertIn((1, 1), idx["TG"])  # "TG" na seq 1, pos 1
        self.assertIn((1, 2), idx["GC"])  # "GC" na seq 1, pos 2

    def test_indexar_banco_vazio(self):
        banco = []
        idx = indexar_banco(banco, k=3)
        self.assertEqual(idx, {})

    # ---------- encontrar_hits ----------
    def test_encontrar_hits_basico(self):
        consulta = "ATAT"
        banco = ["GATAT", "ATGAT"]
        k = 2
        idx = indexar_banco(banco, k)
        hits = encontrar_hits(consulta, idx, k)
        # hits: lista de tuples (i_query, id_seq, pos)
        # Espera-se que "AT" (i=0 e i=2) ocorra em ambas as sequências
        self.assertTrue(any(h[0] == 0 for h in hits))  # i_query=0
        self.assertTrue(any(h[0] == 2 for h in hits))  # i_query=2
        self.assertTrue(any(h[1] == 0 for h in hits))  # id_seq=0
        self.assertTrue(any(h[1] == 1 for h in hits))  # id_seq=1

    def test_encontrar_hits_sem_correspondencias(self):
        consulta = "AAAA"
        banco = ["TTTTTT", "CCCCCC"]
        idx = indexar_banco(banco, k=2)
        hits = encontrar_hits(consulta, idx, k=2)
        self.assertEqual(hits, [])

    # ---------- estender_hit ----------
    def test_estender_hit_perfeito_seed_centro(self):
        consulta = "ATGCATGCA"
        sequencia = "ATGCATGCA"
        k = 3
        # escolher o seed "CAT" ao centro (i=3)
        score, aq, as_ = estender_hit(consulta, sequencia, i_query=3, i_seq=3, k=k, match=2, mismatch=-1)
        # Para sequências iguais, extensão cobre tudo e só há matches:
        # seed: 3*2 = 6; à esquerda: 3*2 = 6; à direita: 3*2 = 6; total 18
        self.assertEqual(score, 18)
        self.assertEqual(aq, consulta)
        self.assertEqual(as_, sequencia)

    def test_estender_hit_com_mismatches(self):
        consulta = "ATGCATGCA"
        sequencia = "ATGCGTGCA"  # mismatch no meio (CAT vs GTG)
        k = 3
        # usar seed "ATG" no início (i=0)
        score, aq, as_ = estender_hit(consulta, sequencia, i_query=0, i_seq=0, k=k, match=2, mismatch=-1)
        self.assertEqual(aq[:3], "ATG")
        self.assertEqual(as_[:3], "ATG")
        self.assertEqual(len(aq), len(as_))  # alinhamentos devem ter o mesmo comprimento

    # ---------- blast_simplificado ----------
    def test_blast_simplificado_ordenacao_e_melhor_hit(self):
        query = "ATGCATGCA"
        banco = ["ATGCATGCA", "TTGCATGGA", "ATGCGTACA"]
        resultados = blast_simplificado(query, banco, k=3)
        # Melhor deve ser contra a sequência id=0 (igual à query)
        score0, id0, q0, s0 = resultados[0]
        self.assertEqual(id0, 0)
        self.assertEqual(q0, query)
        self.assertEqual(s0, banco[0])
        # Scores desc
        scores = [r[0] for r in resultados]
        self.assertEqual(scores, sorted(scores, reverse=True))

    def test_blast_simplificado_sem_hits(self):
        query = "AAAAAA"
        banco = ["TTTTTT", "CCCCCC"]
        resultados = blast_simplificado(query, banco, k=3)
        self.assertEqual(resultados, [])

    def test_blast_simplificado_banco_vazio(self):
        query = "ATGCATGCA"
        resultados = blast_simplificado(query, [], k=3)
        self.assertEqual(resultados, [])


if __name__ == "__main__":
    unittest.main()