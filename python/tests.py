import unittest

from fields import Fq, Fq2, Fq6, Fq12
from ec import EC
import params


class TestFields(unittest.TestCase):
    def setUp(self):
        self.a1 = Fq(3, 11)
        self.b1 = Fq(4, 11)
        self.c1 = Fq(7, 11)
        self.d1 = Fq(12, 11)
        self.e1 = Fq(2, 11)
        self.f1 = Fq(9, 11)
        self.g1 = Fq(8, 11)

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
        

class TestEC(unittest.TestCase):
    def setUp(self):
        self.a1 = Fq(48, 199)
        self.b1 = Fq(50, 199)
        self.c1 = Fq(62, 199)
        self.A1 = EC.get_point_from_x(self.a1)
        self.B1 = EC.get_point_from_x(self.b1)

        self.a2 = Fq2(self.b1, self.a1)
        self.b2 = Fq2(self.b1, self.c1)
        self.A2 = EC.get_point_from_x(self.a2)
        self.B2 = EC.get_point_from_x(self.b2)

        self.g1 = EC.from_affine(Fq(params.g1_x, params.q), Fq(params.g1_y, params.q))
        self.g2 = EC.from_affine(Fq2(Fq(params.g2_x0, params.q), Fq(params.g2_x1, params.q)),
                                 Fq2(Fq(params.g2_y0, params.q), Fq(params.g2_y1, params.q)))

    def test_EC_over_Fq(self):
        self.assertTrue(self.A1.is_on_curve())
        self.assertTrue((self.A1 - self.A1).is_infinity())
        self.assertEqual((self.A1 + EC.infinity(self.a1)), self.A1)
        self.assertEqual(self.A1 + self.A1 + self.A1 + self.B1 + self.B1, 3 * self.A1 + 2 * self.B1)

    def test_EC_over_Fq2(self):
        self.assertTrue(self.A2.is_on_curve())
        self.assertTrue((self.A2 - self.A2).is_infinity())
        self.assertEqual((self.A2 + EC.infinity(self.a2)), self.A2)
        self.assertEqual(self.A2 + self.A2 + self.A2 + self.B2 + self.B2, 3 * self.A2 + 2 * self.B2)

    def test_g1_on_curve(self):
        self.assertTrue(self.g1.is_on_curve())

    def test_g2_on_curve(self):
        self.assertTrue(self.g2.is_on_curve())


if __name__ == '__main__':
    unittest.main()
