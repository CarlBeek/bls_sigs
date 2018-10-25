class Fq(int):
    def __new__(cls, x: int, q: int):
        x = x % q
        ret = super().__new__(cls, x)
        ret.q = q
        return ret

    def __str__(self):
        return super().__str__()

    def __iter__(self):
        return iter([self])

    def __len__(self):
        return 1

    def __add__(self, other):
        return Fq(super().__add__(other), self.q)

    def __radd__(self, other):
        return self.__add__(other)

    def __neg__(self):
        return Fq(super().__neg__(), self.q)

    def __sub__(self, other):
        return Fq(super().__sub__(other), self.q)

    def __rsub__(self, other):
        return Fq(super().__sub__(other).__neg__(), self.q)

    def __mul__(self, other):
        return Fq(super().__mul__(other), self.q)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        other = Fq(other, self.q)
        return Fq(self * other.inverse(), self.q)

    def __rdiv__(self, other):
        return Fq(other * self.inverse(), self.q)

    def __pow__(self, power):
        # Basic square and multiply algorithm
        power = bin(int(power))
        power = power[2:]  # Removes '0b' from number
        ret = self.__class__.one(self.q)
        for i in power:
            ret = ret.square()
            if i == '1':
                ret *= self
        return ret

    def inverse(self):
        t = 0
        new_t = 1
        r = self.q
        new_r = self

        while new_r:
            q = r // new_r
            t, new_t = new_t, t - q * new_t
            r, new_r = new_r, r - q * new_r
        return Fq(t, self.q)

    def sqrt(self):
        # Simplified Tonelli-Shanks for q==3 mod 4
        # https://eprint.iacr.org/2012/685.pdf  Algorithm 2
        assert self.q % 4 == 3
        a1 = self ** ((self.q - 3) // 4)  # Todo: this (maybe) can be sped up with frobrenius endos
        a0 = a1.square() * self
        if a0 == Fq(-1, a0.q):
            raise ArithmeticError('The square root does not exist.')
        return a1*self

    def square(self):
        return self * self

    def is_zero(self):
        return self == 0

    def is_one(self):
        return self == 1

    @property
    def degree(self):
        return 1

    @classmethod
    def zero(cls, q):
        return cls(0, q)

    @classmethod
    def one(cls, q):
        return cls(1, q)

    @classmethod
    def to_cls(cls, obj, q):
        if isinstance(obj, int):
            return cls(obj, q)
        raise NotImplementedError


class ExtensionField(tuple):
    def __new__(cls, elems, q):
        elems = (Fq(elem, q) for elem in elems)
        ret = super().__new__(cls, elems)
        ret.q = q
        ret.degree = len(ret)
        return ret

    def __neg__(self):
        return self.__class__((-a for a in self), self.q)

    def __add__(self, other):
        other = self.__class__.to_cls(other, self.degree, self.q)
        c = (c_self + c_other for c_self, c_other in zip(self, other))
        return self.__class__(c, self.q)

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self + other.__neg__()

    def __rsub__(self, other):
        return self.__neg__() + other

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        return self * other.inverse()

    def __rdiv__(self, other):
        return other * self.inverse()

    def __pow__(self, power):
        # Basic square and multiply algorithm
        power = bin(int(power))
        power = power[2:]  # Removes '0b' from number
        ret = self.__class__.one(self.degree, self.q)
        for i in power:
            ret = ret.square()
            if i == '1':
                ret *= self
        return ret

    def __str__(self):
        return 'Fq%i' % self.degree + super().__str__()

    def __gt__(self, other):
        for value in zip(reversed(self), reversed(other)):
            if value[0] == value[1]:
                continue
            else:
                return value[0] > value[1]
        return False

    def __lt__(self, other):
        return other > self

    def __eq__(self, other):
        other = self.to_cls(other, self.degree, self.q)
        return super().__eq__(other)

    def square(self):
        return self * self

    @classmethod
    def zero(cls, degree, q):
        return cls.__new__((0,) * degree, q)

    @classmethod
    def one(cls, degree, q):
        return cls((1,) + (0,) * (degree-1), q)

    @classmethod
    def to_cls(cls, obj, degree, q):
        if isinstance(obj, cls):
            return obj
        elif isinstance(obj, (Fq, int)):
            obj = (obj,)
        return cls((obj + (0,)*(degree - len(obj))), q)



class Fq2(ExtensionField):
    def __mul__(self, other):
        other = self.to_cls(other, self.degree, self.q)
        aa = self[0] * other[0]
        bb = self[1] * other[1]
        o = other[0] + other[1]
        c1 = self[1] + self[0]
        c1 *= o
        c1 -= aa
        c1 -= bb
        c0 = aa - bb
        return self.__class__((c0, c1), self.q)

    def mul_by_nonresidue(self):
        return Fq2((self[0] - self[1], self[0] + self[1]), self.q)

    def inverse(self):
        t1 = self[1] * self[1]
        t0 = self[0] * self[0]
        t0 += t1
        t0 = t0.inverse()
        a = self[0]*t0
        b = self[1]*t0
        b = -b
        return self.__class__((a, b), self.q)

    def sqrt(self):
        # Modified Tonelli-Shanks for q==3 mod 4
        # https://eprint.iacr.org/2012/685.pdf  Algorithm 9
        assert self.q % 4 == 3 # q%4 == 3 This can ultimately be removed for BLS12-381
        a1 = self ** ((self.q - 3) // 4) # Todo: this can be sped up with frobrenius endos
        alpha = a1.square() * self
        a0 = alpha**(self.q+1)
        one = Fq2.one(2, self.q)

        if a0 == - one:
            raise ArithmeticError('The square root does not exist.')
        x0 = a1 * self
        if alpha == - one:
            # This returns i*x0 (i=sqrt(-1))
            return self.__class__(0, Fq(-1, self.q).sqrt(), self.q)*x0
        b = (alpha + one)**((self.q - 1)//2) # Todo: this can be sped up with frobrenius endos
        return b * x0

    def frobenius_endo(self, power):
        from params import FROB_FQ2
        return self.__class__(self[0], self[1] * FROB_FQ2[power % 2])


class Fq6(ExtensionField):
    def __mul__(self, other):
        other = self.to_cls(other, self.degree, self.q)
        aa, bb, cc = self.c0, self.c1, self.c2

        aa *= other.c0
        bb *= other.c1
        cc *= other.c2

        # t1
        t1 = other.c1
        t1 += other.c2
        tmp = self.c1 + self.c2
        t1 *= tmp
        t1 -= bb
        t1 -= cc
        t1 = t1.mul_by_nonresidue()
        t1 += aa

        # t3
        t3 = other.c0
        t3 += other.c2
        tmp = self.c0 + self.c2
        t3 *= tmp
        t3 -= aa
        t3 += bb
        t3 -= cc

        # t2
        t2 = other.c0
        t2 += other.c1
        tmp = self.c0 + self.c1
        t2 *= tmp
        t2 -= aa
        t2 -= bb
        t2 += cc.mul_by_nonresidue()

        return self.__class__((*t1, *t2, *t3), self.q)

    def inverse(self):
        c0 = self.c2
        c0 = c0.mul_by_nonresidue()
        c0 *= self.c1
        c0 = -c0
        c0 += self.c0*self.c0

        c1 = self.c2
        c1 *= c1
        c1 = c1.mul_by_nonresidue()
        c1 -= self.c0*self.c1

        c2 = self.c1
        c2 = c2.square()
        c2 -= self.c0*self.c2

        tmp1 = self.c2*c1
        tmp2 = self.c1*c2
        tmp1 += tmp2
        tmp1 = tmp1.mul_by_nonresidue()
        tmp2 = self.c0*c0
        tmp1 += tmp2

        tmp1 = tmp1.inverse()
        if tmp1 == 0:
            raise ArithmeticError
        c0 *= tmp1
        c1 *= tmp1
        c2 *= tmp1
        return self.__class__((*c0, *c1, *c2), self.q)

    def mul_by_nonresidue(self):
        c1, c2, c0 = self.c0, self.c1, self.c2
        c0 = c0.mul_by_nonresidue()
        return self.__class__((*c0, *c1, *c2), self.q)

    def frobenius_endo(self, power):
        from params import FROB_FQ6_C1, FROB_FQ6_C2
        c0 = self.c0.frobenius_endo(power)
        c1 = self.c1.frobenius_endo(power)
        c2 = self.c2.frobenius_endo(power)

        c1 *= FROB_FQ6_C1[power % 6]
        c2 *= FROB_FQ6_C2[power % 6]
        return self.__class__((*c0, *c1, *c2), self.q)

    @property
    def c0(self):
        return Fq2(self[0:2], self.q)

    @property
    def c1(self):
        return Fq2(self[2:4], self.q)

    @property
    def c2(self):
        return Fq2(self[4:6], self.q)



class Fq12(ExtensionField):
    def __mul__(self, other):
        other = self.to_cls(other, self.degree, self.q)
        aa = self.c0 * other.c0
        bb = self.c1 * other.c1
        o = other.c0 + other.c1
        c1 = self.c1 + self.c0
        c1 *= o
        c1 -= aa
        c1 -= bb
        c0 = bb.mul_by_nonresidue()
        c0 += aa
        return self.__class__((*c0, *c1), self.q)

    def inverse(self):
        c0 = self.c0 * self.c0
        c1 = self.c1 * self.c1
        c1 = c1.mul_by_nonresidue()
        c0 -= c1

        t = c0.inverse()
        t0 = t * self.c0
        t1 = t * self.c1
        t1 = -t1
        return self.__class__((*t0, *t1), self.q)

    def frobenius_endo(self, power):
        from params import FROB_FQ12_C1
        c0 = self.c0.frobenius_endo(power)
        c1 = self.c1.frobenius_endo(power)

        c1c0 = c1.c0 * FROB_FQ12_C1[power % 12]
        c1c1 = c1.c1 * FROB_FQ12_C1[power % 12]
        c1c2 = c1.c2 * FROB_FQ12_C1[power % 12]

        return self.__class__((*c0, *c1c0, *c1c1, *c1c2), self.q)

    @property
    def c0(self):
        return Fq6(self[0:6], self.q)

    @property
    def c1(self):
        return Fq6(self[6:12], self.q)
