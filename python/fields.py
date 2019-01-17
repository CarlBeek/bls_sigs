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
        other = Fq.to_cls(other, self.q)
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
        a1 = self ** ((self.q - 3) // 4)  # Todo: this can be spead up with frobrenius endos
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

    @classmethod
    def zero(cls, q):
        return cls(0, q)

    @classmethod
    def one(cls, q):
        return cls(1, q)

    @classmethod
    def all_one_poly(cls, q):
        return cls.one(q)

    @classmethod
    def to_cls(cls, obj, q):
        if isinstance(obj, cls):
            return obj
        elif isinstance(obj, int):
            return cls(obj, q)
        raise NotImplementedError


class Fq2(tuple):
    def __new__(cls, c0: Fq, c1: Fq):
        return super().__new__(cls, (c0, c1))

    def __add__(self, other):
        other = self.to_cls(other, self.q)
        return self.__class__(self.c0 + other.c0, self.c1 + other.c1)

    def __radd__(self, other):
        return self.__add__(other)

    def __neg__(self):
        return self.__class__(-self.c0, -self.c1)

    def __sub__(self, other):
        other = self.to_cls(other, self.q)
        return self + - other

    def __rsub__(self, other):
        other = self.to_cls(other, self.q)
        return other + - self

    def __mul__(self, other):
        other = self.to_cls(other, self.q)
        aa = self.c0 * other.c0
        bb = self.c1 * other.c1
        o = other.c0 + other.c1
        c1 = self.c1 + self.c0
        c1 *= o
        c1 -= aa
        c1 -= bb
        c0 = aa - bb
        return self.__class__(c0, c1)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        other = self.to_cls(other, self.q)
        return self * other.inverse()

    def __rdiv__(self, other):
        other = self.to_cls(other, self.q)
        return other * self.inverse()

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

    def __gt__(self, other):
        other = self.to_cls(other, self.q)
        if self.c1 > other.c1:
            return True
        elif self.c1 < other.c1:
            return False
        elif self.c0 > other.c0:
            return True
        return False

    def __lt__(self, other):
        other = self.to_cls(other, self.q)
        return not (self > other and self == other)

    def __str__(self):
        return 'Fq2(' + str(self.c0) + ' + ' + str(self.c1) + ' * u)'

    def mul_by_nonresidue(self):
        return Fq2(self.c0 - self.c1, self.c0 + self.c1)

    def inverse(self):
        t1 = self.c1 * self.c1
        t0 = self.c0 * self.c0
        t0 += t1
        t0 = t0.inverse()
        a = self.c0*t0
        b = self.c1*t0
        b = -b
        return self.__class__(a, b)

    def square(self):
        return self * self

    def sqrt(self):
        # Modified Tonelli-Shanks for q==3 mod 4
        # https://eprint.iacr.org/2012/685.pdf  Algorithm 9
        assert self.q % 4 == 3 # q%4 == 3 This can ultimately be removed for BLS12-381
        a1 = self ** ((self.q - 3) // 4) # Todo: this can be spead up with frobrenius endos
        alpha = a1.square() * self
        a0 = alpha**(self.q+1)

        if a0 == -Fq2.one(self.q):
            raise ArithmeticError('The square root does not exist.')
        x0 = a1 * self
        if alpha == -Fq2.one(self.q):
            # This returns i*x0 (i=sqrt(-1))
            return self.__class__(Fq(0, self.q), Fq(-1, self.q).sqrt())*x0
        b = (alpha + 1)**((self.q - 1)//2) # Todo: this can be spead up with frobrenius endos
        return b * x0

    def frobenius_endo(self, power):
        from params import FROB_FQ2
        return self.__class__(self.c0, self.c1 * FROB_FQ2[power % 2])

    def is_zero(self):
        return self.c0.is_zero() and self.c1.is_zero()

    def is_one(self):
        return self.c0.is_one() and self.c1.is_zero()

    @classmethod
    def zero(cls, q):
        return cls(Fq.zero(q), Fq.zero(q))

    @classmethod
    def one(cls, q):
        return cls.to_cls(Fq.one(q), q)

    @classmethod
    def all_one_poly(cls, q):
        return cls(Fq.one(q), Fq.one(q))

    @classmethod
    def to_cls(cls, obj, q):
        if isinstance(obj, cls):
            return obj
        elif isinstance(obj, (int, Fq)):
            return cls(Fq.to_cls(obj, q), Fq.zero(q))
        raise NotImplementedError

    @property
    def q(self):
        return self.c0.q

    @property
    def c0(self):
        return self[0]
    
    @property
    def c1(self):
        return self[1]


class Fq6(tuple):
    def __new__(cls, c0, c1, c2):
        return super().__new__(cls, (c0, c1, c2))

    def __neg__(self):
        c = (-c_self for c_self in self)
        return self.__class__(*c)

    def __add__(self, other):
        other = self.to_cls(other, self.q)
        c = (c_self + c_other for c_self, c_other in zip(self, other))
        return self.__class__(*c)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        other = self.to_cls(other, self.q)
        return self + - other

    def __rsub__(self, other):
        return other + - self

    def __mul__(self, other):
        other = self.to_cls(other, self.q)
        aa, bb, cc = self

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

        return self.__class__(t1, t2, t3)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        other = self.to_cls(other, self.q)
        return self * other.inverse()

    def __rdiv__(self, other):
        return other * self.inverse()

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

    def __gt__(self, other):
        other = self.to_cls(other, self.q)
        if self.c2 > other.c2:
            return True
        elif self.c2 < other.c2:
            return False
        elif self.c1 > other.c1:
            return True
        elif self.c1 < other.c1:
            return False
        elif self.c0 > other.c0:
            return True
        return False

    def __lt__(self, other):
        other = self.to_cls(other, self.q)
        return not (self > other and self == other)

    def __str__(self):
        return 'Fq6(' + str(self.c0) + ' + ' + str(self.c1) + ' * v + ' + str(self.c2) + ' * v^2)'

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
        if tmp1.is_zero():
            raise ArithmeticError
        c0 *= tmp1
        c1 *= tmp1
        c2 *= tmp1
        return self.__class__(c0, c1, c2)

    def square(self):
        return self * self

    def mul_by_nonresidue(self):
        c1, c2, c0 = self
        c0 = c0.mul_by_nonresidue()
        return self.__class__(c0, c1, c2)

    def frobenius_endo(self, power):
        from params import FROB_FQ6_C1, FROB_FQ6_C2
        c0 = self.c0.frobenius_endo(power)
        c1 = self.c1.frobenius_endo(power)
        c2 = self.c2.frobenius_endo(power)

        c1 *= FROB_FQ6_C1[power % 6]
        c2 *= FROB_FQ6_C2[power % 6]
        return self.__class__(c0, c1, c2)

    def is_zero(self):
        return self.c0.is_zero() and self.c1.is_zero() and self.c2.is_zero()

    def is_one(self):
        return self.c0.is_one() and self.c1.is_zero() and self.c2.is_zero()

    @classmethod
    def zero(cls, q):
        return cls(Fq2.zero(q), Fq2.zero(q), Fq2.zero(q))

    @classmethod
    def one(cls, q):
        return cls.to_cls(Fq.one(q), q)

    @classmethod
    def to_cls(cls, obj, q):
        if isinstance(obj, cls):
            return obj
        elif isinstance(obj, (int, Fq, Fq2)):
            return cls(Fq2.to_cls(obj, q), Fq2.zero(q), Fq2.zero(q))
        raise NotImplementedError

    @property
    def q(self):
        return self[0].q

    @property
    def c0(self):
        return self[0]

    @property
    def c1(self):
        return self[1]

    @property
    def c2(self):
        return self[2]


class Fq12(tuple):
    def __new__(cls, c0, c1):
        return super().__new__(cls, (c0, c1))

    def __neg__(self):
        c = (-c_self for c_self in self)
        return self.__class__(*c)

    def __add__(self, other):
        other = self.to_cls(other, self.q)
        c = (c_self + c_other for c_self, c_other in zip(self, other))
        return self.__class__(*c)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self + -other

    def __rsub__(self, other):
        return other + -self

    def __mul__(self, other):
        other = self.to_cls(other, self.q)
        aa = self.c0 * other.c0
        bb = self.c1 * other.c1
        o = other.c0 + other.c1
        c1 = self.c1 + self.c0
        c1 *= o
        c1 -= aa
        c1 -= bb
        c0 = bb.mul_by_nonresidue()
        c0 += aa
        return self.__class__(c0, c1)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        other = self.to_cls(other, self.q)
        return self * other.inverse()

    def __rdiv__(self, other):
        return other * self.inverse()

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

    def __gt__(self, other):
        other = self.to_cls(other, self.q)
        if self.c1 > other.c1:
            return True
        elif self.c1 < other.c1:
            return False
        elif self.c0 > other.c0:
            return True
        return False

    def __lt__(self, other):
        other = self.to_cls(other, self.q)
        return not (self > other and self == other)

    def __str__(self):
        return 'Fq12(' + str(self.c0) + ' + ' + str(self.c1) + ' * w)'

    def inverse(self):
        c0 = self.c0 * self.c0
        c1 = self.c1 * self.c1
        c1 = c1.mul_by_nonresidue()
        c0 -= c1

        t = c0.inverse()
        t0 = t * self.c0
        t1 = t * self.c1
        t1 = -t1
        return self.__class__(t0, t1)

    def square(self):
        return self * self

    def frobenius_endo(self, power):
        from params import FROB_FQ12_C1
        c0 = self.c0.frobenius_endo(power)
        c1 = self.c1.frobenius_endo(power)

        c1c0 = c1.c0 * FROB_FQ12_C1[power % 12]
        c1c1 = c1.c1 * FROB_FQ12_C1[power % 12]
        c1c2 = c1.c2 * FROB_FQ12_C1[power % 12]

        return self.__class__(c0, Fq6(c1c0, c1c1, c1c2))

    def is_zero(self):
        return self.c0.is_zero() and self.c1.is_zero()

    def is_one(self):
        return self.c0.is_one() and self.c1.is_zero()

    @classmethod
    def zero(cls, q):
        return cls(Fq6.zero(q), Fq6.zero(q))

    @classmethod
    def one(cls, q):
        return cls.to_cls(Fq.one(q), q)

    @classmethod
    def to_cls(cls, obj, q):
        if isinstance(obj, cls):
            return obj
        elif isinstance(obj, (int, Fq, Fq2, Fq6)):
            return cls(Fq6.to_cls(obj, q), Fq6.zero(q))
        raise NotImplementedError

    @property
    def q(self):
        return self[0].q

    @property
    def c0(self):
        return self[0]

    @property
    def c1(self):
        return self[1]


def shitty_print(elem):
    if isinstance(elem, Fq):
        return elem
    else:
        return [shitty_print(s) for s in elem]

if __name__ == '__main__':
    q = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
    # a = Fq6(Fq2(Fq(123, q), Fq(456, q).inverse()), Fq2(Fq(987345897, q), Fq(7235986723, q)), Fq2(Fq(82635739, q), Fq(37854362965, q))).inverse() + (Fq6(Fq2(Fq(358927054890, q), Fq(1640628, q)).inverse(), Fq2.zero(q), Fq2(Fq(0, q), Fq(1, q)))*Fq6(Fq2(Fq(0, q), Fq(1, q)), Fq2.zero(q), Fq2.zero(q)))
    # b = (((a**97843289).inverse() + 100) * a + 10983420934987321879312890312879786321897123879629873698471).inverse()

    a = Fq6(Fq2.all_one_poly(q), Fq2.all_one_poly(q), Fq2.all_one_poly(q))
    b = Fq6(Fq2(Fq.zero(q), Fq(1, q)), Fq2.zero(q), Fq2.zero(q))
    print(a)
    print(b)
    print(a*b)

    a2 = Fq2(Fq(8953249078321897032109123578, q), Fq(18796543789076543, q))
    print(a2.frobenius_endo(1))
    # print(Fq(2, q).inverse())

    b = Fq6(Fq2(Fq(1391480139073249803319359180048344339509516026819802451214654929250728703234810770821865915244430291176113392586742, q),
                Fq(1612317712542663890081113311376019869718629909620604387996387280234976046649155112835436167973344680809449978882112, q)),
            Fq2(Fq(3220862940368626293296014966215844512046684494314082373471564105700750447097473715456639923414620241478699983825925, q),
                Fq(2369452030556451996763222640647413327139789520135959582970792194712267126750568809654671258200824218662919681321530, q)),
            Fq2(Fq(3190587454610077223176343295169157502250630032730423831188786790767781681364466089654706471019908794689950716541218, q),
                Fq(3398916547685317465766376584328149757439489843737047752629685847650026248993628549908476031583435785800225436201570, q)))

    print(b**(q**3))
    print(shitty_print(b**(q**2)))

    c = Fq12(Fq6(Fq2(Fq(297337868834170793670139122744378726479201867292723741905403927492789473149162169485284037341440980065776726951721, q),
                     Fq(622274856095598354639655456495881934563293727989726050622905348990239346286607670888369985633840340011707251042506, q)),
                 Fq2(Fq(947253790885885905586255606557418403407380972814444400312716434855027996393348706322495893345981307600323372846736, q),
                     Fq(1461719620741406082999455235802950617408296494933304225929075519545122340579389886772661651525168876967652679327665,q)),
                 Fq2(Fq(2122739563431538881275873484801695771009172888776844144742191441698926455264153592718300047006552948158807911143091, q),
                     Fq(1508606287621993082785371759790924577047145322726831795904184556180235695248233616406433991455567053134294259100296, q))),
             Fq6(Fq2(Fq(297337868834170793670139122744378726479201867292723741905403927492789473149162169485284037341440980065776726951721, q),
                     Fq(622274856095598354639655456495881934563293727989726050622905348990239346286607670888369985633840340011707251042506, q)),
                 Fq2(Fq(947253790885885905586255606557418403407380972814444400312716434855027996393348706322495893345981307600323372846736, q),
                     Fq(1461719620741406082999455235802950617408296494933304225929075519545122340579389886772661651525168876967652679327665, q)),
                 Fq2(Fq(2122739563431538881275873484801695771009172888776844144742191441698926455264153592718300047006552948158807911143091, q),
                     Fq(1508606287621993082785371759790924577047145322726831795904184556180235695248233616406433991455567053134294259100296, q))))

    print([shitty_print(c**(q**i)) for i in range(12)])
