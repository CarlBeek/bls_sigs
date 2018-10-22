import unittest
from random import uniform

from fields import Fq, Fq2, Fq6, Fq12
from ec import EC, TwistedEC
from params import q, r, g1_x, g1_y, g2_x, g2_y, FROB_FQ2, FROB_FQ6_C1, FROB_FQ6_C2, FROB_FQ12_C1
from BLS import sign, verify, key_gen, compress, decompress
from paring import twist, untwist, paring


class TestFields(unittest.TestCase):
    def setUp(self):
        self.a1 = Fq(3, 11)
        self.b1 = Fq(4, 11)
        self.c1 = Fq(7, 11)
        self.d1 = Fq(12, 11)
        self.e1 = Fq(2, 11)
        self.f1 = Fq(9, 11)
        self.g1 = Fq(6, 11)

        self.a2 = Fq2(self.a1, self.b1)
        self.b2 = Fq2(self.b1, self.c1)
        self.c2 = Fq2(self.c1, self.d1)
        self.d2 = Fq2(self.d1, self.e1)
        self.e2 = Fq2(self.e1, self.f1)
        self.f2 = Fq2(self.f1, self.g1)

        self.a6 = Fq6(self.a2, self.b2, self.c2)
        self.b6 = Fq6(self.b2, self.c2, self.d2)
        self.c6 = Fq6(self.d2, self.e2, self.f2)

        self.a12 = Fq12(self.a6, self.b6)
        self.b12 = Fq12(self.b6, self.c6)

    def test_Fq(self):
        self.assertEqual(self.a1 + self.a1, self.a1 * 2)
        self.assertTrue((self.a1 - self.a1).is_zero())
        self.assertEqual(self.a1 + self.a1 + self.a1 + self.b1 + self.b1,  3 * self.a1 + 2 * self.b1)
        self.assertTrue((self.a1 * self.b1 / self.b1 / self.a1).is_one())
        self.assertEqual(self.a1.sqrt().square(), self.a1)
        self.assertEqual(self.a1 ** 3, self.a1 * self.a1 * self.a1)

    def test_Fq2(self):
        self.assertEqual(self.a2 + self.a2, self.a2 * 2)
        self.assertTrue((self.a2 - self.a2).is_zero())
        self.assertEqual(self.a2 + self.a2 + self.a2 + self.b2 + self.b2, 3 * self.a2 + 2 * self.b2)
        self.assertTrue((self.a2 * self.b2 / self.b2 / self.a2).is_one())
        self.assertEqual(self.a2.sqrt().square(), self.a2)
        self.assertEqual(self.a2 ** 3, self.a2 * self.a2 * self.a2)

    def test_Fq6(self):
        self.assertEqual(self.a6 + self.a6, self.a6 * 2)
        self.assertTrue((self.a6 - self.a6).is_zero())
        self.assertEqual(self.a6 + self.a6 + self.a6 + self.b6 + self.b6, 3 * self.a6 + 2 * self.b6)
        self.assertTrue((self.a6 * self.b6 / self.b6 / self.a6).is_one())
        self.assertEqual(self.a6 ** 3, self.a6 * self.a6 * self.a6)

    def test_Fq12(self):
        self.assertEqual(self.a12 + self.a12, self.a12 * 2)
        self.assertTrue((self.a12 - self.a12).is_zero())
        self.assertEqual(self.a12 + self.a12 + self.a12 + self.b12 + self.b12, 3 * self.a12 + 2 * self.b12)
        self.assertTrue((self.a12 * self.b12 / self.b12 / self.a12).is_one())
        self.assertEqual(self.a12 ** 3, self.a12 * self.a12 * self.a12)


class TestParams(unittest.TestCase):
    def test_frobenius_coefficients(self):
        all_one_poly = Fq2(Fq.one(q), Fq.one(q))
        self.assertEqual((Fq(1, q), Fq(-1, q)), FROB_FQ2)
        self.assertEqual(tuple(all_one_poly ** (((q ** i) - 1) // 3) for i in range(0, 6)), FROB_FQ6_C1)
        self.assertEqual(tuple(all_one_poly ** (((2 * q ** i) - 2) // 3) for i in range(0, 6)), FROB_FQ6_C2)
        self.assertEqual(tuple(all_one_poly ** (((q ** i) - 1) // 6) for i in range(0, 12)), FROB_FQ12_C1)


class TestEC(unittest.TestCase):
    def setUp(self):
        self.a1 = Fq(48, 199)
        self.b1 = Fq(50, 199)
        self.c1 = Fq(62, 199)
        self.A1 = EC.get_point_from_x(self.a1)
        self.B1 = EC.get_point_from_x(self.b1)

        self.a2 = Fq2(self.b1, self.a1)
        self.b2 = Fq2(self.b1, self.c1)
        self.A2 = TwistedEC.get_point_from_x(self.a2)
        self.B2 = TwistedEC.get_point_from_x(self.b2)

        self.g1 = EC.from_affine(g1_x, g1_y)
        self.g2 = TwistedEC.from_affine(g2_x, g2_y)

    def test_EC_over_Fq(self):
        self.assertTrue(self.A1.is_on_curve())
        self.assertEqual(self.A1 + self.A1, self.A1.double())
        self.assertTrue((self.A1 + self.A1).is_on_curve())
        self.assertTrue(self.A1.double().is_on_curve())
        self.assertTrue((self.A1 * 100).is_on_curve())
        self.assertTrue((self.A1 - self.A1).is_infinity())
        self.assertEqual((self.A1 + EC.infinity(self.a1)), self.A1)
        self.assertEqual(self.A1 + self.A1 + self.A1 + self.B1 + self.B1, 3 * self.A1 + 2 * self.B1)

    def test_twisted_EC_over_Fq2(self):
        self.assertTrue(self.A2.is_on_curve())
        self.assertEqual(self.A2 + self.A2, self.A2.double())
        self.assertTrue((self.A2 + self.A2).is_on_curve())
        self.assertTrue(self.A2.double().is_on_curve())
        self.assertTrue((self.A2*100).is_on_curve())
        self.assertTrue((self.A2 - self.A2).is_infinity())
        self.assertEqual((self.A2 + EC.infinity(self.a2)), self.A2)
        self.assertEqual(self.A2 + self.A2 + self.A2 + self.B2 + self.B2, 3 * self.A2 + 2 * self.B2)

    def test_g1_on_curve(self):
        self.assertTrue(self.g1.is_on_curve())
        self.assertTrue((self.g1 * 5).is_on_curve())

    def test_g2_on_curve(self):
        self.assertTrue(self.g2.is_on_curve())
        self.assertTrue((self.g2 * 5).is_on_curve())


class TestParing(unittest.TestCase):
    def setUp(self):
        self.g1 = EC.from_affine(g1_x, g1_y)
        self.g2 = TwistedEC.from_affine(g2_x, g2_y)

        self.sk_0 = int(uniform(1, r))
        self.sk_1 = int(uniform(1, r))
        self.pk_0 = key_gen(self.sk_0)
        self.pk_1 = key_gen(self.sk_1)
        self.msg_0 = self.g2 * int(uniform(1, r))
        self.msg_1 = self.g2 * int(uniform(1, r))
        self.sigma_0_0 = sign(self.msg_0, self.sk_0)
        self.sigma_0_1 = sign(self.msg_0, self.sk_1)
        self.sigma_1_0 = sign(self.msg_1, self.sk_0)
        self.sigma_1_1 = sign(self.msg_1, self.sk_1)

    def test_twisting(self):
        self.assertEqual(self.g2, twist(untwist(self.g2)))

    def test_pk_on_curve(self):
        self.assertTrue(self.pk_0.is_on_curve())
        self.assertTrue(self.pk_1.is_on_curve())

    def test_msg_on_curve(self):
        self.assertTrue(self.msg_0.is_on_curve())
        self.assertTrue(self.msg_1.is_on_curve())

    def test_signature_is_on_curve(self):
        self.assertTrue(self.msg_0.is_on_curve())
        self.assertTrue(self.msg_1.is_on_curve())

    def test_paring_func(self):
        self.assertEqual(paring(self.g1 * self.sk_0, self.g2, True), paring(self.g1, self.g2 * self.sk_0, True))

    def test_invalid_sig(self):
        self.assertFalse(verify(self.msg_0, self.sigma_0_0, self.pk_1))
        self.assertFalse(verify(self.msg_0, self.sigma_1_1, self.pk_0))
        self.assertFalse(verify(self.msg_0, self.sigma_1_1, self.pk_1))
        self.assertFalse(verify(self.msg_1, self.sigma_0_0, self.pk_0))
        self.assertFalse(verify(self.msg_1, self.sigma_0_0, self.pk_1))
        self.assertFalse(verify(self.msg_1, self.sigma_1_1, self.pk_0))

    def test_valid_sig(self):
        self.assertTrue(verify(self.msg_0, self.sigma_0_0, self.pk_0))
        self.assertTrue(verify(self.msg_1, self.sigma_1_1, self.pk_1))

    def test_sigs_aggregate(self):
        self.assertTrue(verify(self.msg_0, self.sigma_0_0 + self.sigma_0_1, self.pk_0 + self.pk_1))

    def test_point_compression(self):
        self.assertEqual(self.sigma_0_0, decompress(compress(self.sigma_0_0)))
        self.assertEqual(self.sigma_1_1, decompress(compress(self.sigma_1_1)))


if __name__ == '__main__':
    unittest.main()
