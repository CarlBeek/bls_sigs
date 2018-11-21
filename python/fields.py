q = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
# Fq2_poly_coeffs = (1, 0)
# Fq12_poly_coeffs =


def fq_inv(x, q):
    # Extended euclidean alg
    t = 0
    new_t = 1
    r = q
    new_r = x
    while new_r:
        q = r // new_r
        t, new_t = new_t, t - q * new_t
        r, new_r = new_r, r - q * new_r
    return t

class Fq(object):
    def __init__(self, n):
        self.n = n

    def __add__(self, other):
        return Fq((self.n + other.n) % q)

    def __radd__(self, other):
        return self.__add__(other)

    def __neg__(self):
        return Fq((- self.n) % q)

    def __sub__(self, other):
        return Fq((self.n + other.n.__neg__()) % q)

    def __rsub__(self, other):
        return Fq((other.n + self.n.__neg__()) % q)

    def __mul__(self, other):
        return Fq((self.n * other.n) % q)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        return Fq(self * other.inverse())

    def __rdiv__(self, other):
        return Fq(other * self.inverse())

    def __pow__(self, power):
        # Basic square and multiply algorithm
        power = bin(int(power))
        power = power[2:]  # Removes '0b' from number
        ret = 1
        for i in power:
            ret = ret.square()
            if i == '1':
                ret *= self
        return ret

    def __eq__(self, other):
        return self.n == other.n

    def __gt__(self, other):
        return self.n > other.n

    def __lt__(self, other):
        return self.n < other.n

    def inverse(self):
        return Fq(fq_inv(self.n, q))

    def sqrt(self):
        # Simplified Tonelli-Shanks for q==3 mod 4
        # https://eprint.iacr.org/2012/685.pdf  Algorithm 2
        assert q % 4 == 3 # Todo: this can be removed for BLS 12 381
        a1 = self.n ** ((q - 3) // 4)  # Todo: this can be spead up with frobrenius endos
        a0 = a1 * a1 * self.n
        if a0 == Fq(-1 % q):
            raise ArithmeticError('The square root does not exist.')
        return a1*self.n


def poly_long_div(num, denum, q):
    res_order = len(denum) - len(num)
    for i in range(-1, res_order-2, -1):
        mul = num[i] * fq_inv(denum[-1], q)
        for j in range(-2, -len(denum)-1, -1):
            num[i+j+1] -= denum[j]*mul
    return num[res_order-1:]


class Fqx(object):
    def __init__(self, elems):
        elems = tuple(elems) if not isinstance(elems, tuple) else elems # Needed to call len(elems) later (if generator)
        self.elems = elems
        self.deg = len(elems)
        self.mod_poly = [0] # This is intentionally useless as mod the 0 poly is a fixed point

    def __add__(self, other):
        return Fqx(((s + o) % q for s, o in zip(self.elems, other.elems)))

    def __neg__(self):
        return Fqx(((-s) % q for s in self.elems))

    def __sub__(self, other):
        return Fqx(((s - o) % q for s, o in zip(self.elems, other.elems)))

    def __rsub__(self, other):
        return Fqx(((o - s) % q for s, o in zip(self.elems, other.elems)))

    def __mul__(self, other):
        ret = [0]*(self.deg + other.deg -2)
        # Multiply all elememnts of the same order
        for i, elem_self in enumerate(self.elems):
            for j, elem_other in enumerate(other.elems):
                ret[i + j] += elem_self * elem_other
        # Todo: divide by polynomial
        ret = poly_long_div(ret, self.mod_poly)
        return self.__class__((n % q for n in ret))

    def __rmul__(self, other):
        self * other

    def __div__(self, other):
        return self * other.inverse()

    def __truediv__(self, other):
        return self.__div__(other)

    def __rdiv__(self, other):
        return self.inverse() * other

    def __str__(self):
        return 'Fq%i' % self.deg + self.elems.__str__()

    def __gt__(self, other):
        for value in zip(reversed(self.elems), reversed(other.elems)):
            if value[0] == value[1]:
                continue
            else:
                return value[0] > value[1]
        return False

    def __lt__(self, other):
        return other > self

    def __eq__(self, other):
        return self.elems == other.elems

    def __pow__(self, power):
        # Basic square and multiply algorithm
        power = bin(int(power))
        power = power[2:]  # Removes '0b' from number
        ret = 1
        for i in power:
            ret = ret.square()
            if i == '1':
                ret *= self
        return ret

    def int_mul(self, n):
        return self.__class__(((x*n) % q for x in self.elems))

    def inverse(self):
        pass


class Fq2(Fqx):
    def __init__(self, elems):
        self.elems = elems
        self.deg = 2
        self.mod_poly = [0, 1]

    def sqrt(self):
        # Modified Tonelli-Shanks for q==3 mod 4
        # https://eprint.iacr.org/2012/685.pdf  Algorithm 9
        assert q % 4 == 3 # Todo: This can ultimately be removed for BLS12-381
        a1 = self ** ((q - 3) // 4)  # Todo: this can be spead up with frobrenius endos
        alpha = a1.square() * self
        a0 = alpha**(q+1)
        one = Fq2(1, 0)
        if a0 == - one:
            raise ArithmeticError('The square root does not exist.')
        x0 = a1 * self
        if alpha == - one:
            # This returns i*x0 (i=sqrt(-1))
            return Fq2((0, Fq(-1).sqrt().n))*x0
        b = (alpha + one)**((q - 1)//2)  # Todo: this can be spead up with frobrenius endos
        return b * x0

    def frobenius_endo(self, power):
        from params import FROB_FQ2
        return Fq2((self.elems[0], self.elems[1] * FROB_FQ2[power % 2]))


class Fq6(Fqx):
    def frobenius_endo(self, power):
        from params import FROB_FQ6_C1, FROB_FQ6_C2
        c0 = self.c0.frobenius_endo(power)
        c1 = self.c1.frobenius_endo(power)
        c2 = self.c2.frobenius_endo(power)

        c1 *= FROB_FQ6_C1[power % 6]
        c2 *= FROB_FQ6_C2[power % 6]
        return self.__class__((*c0, *c1, *c2))

    @property
    def c0(self):
        return Fq2((self.elems[0:2]))

    @property
    def c1(self):
        return Fq2((self.elems[2:4]))

    @property
    def c2(self):
        return Fq2((self.elems[4:6]))


class Fq12(Fqx):
    def frobenius_endo(self, power):
        from params import FROB_FQ12_C1
        c0 = self.c0.frobenius_endo(power)
        c1 = self.c1.frobenius_endo(power)

        c1c0 = c1.c0 * FROB_FQ12_C1[power % 12]
        c1c1 = c1.c1 * FROB_FQ12_C1[power % 12]
        c1c2 = c1.c2 * FROB_FQ12_C1[power % 12]

        return self.__class__(c0, Fq6(c1c0, c1c1, c1c2))

    @property
    def c0(self):
        return Fq6((self.elems[0:6]))

    @property
    def c1(self):
        return Fq6((self.elems[6:12]))



def costly_func():
    a = Fq2((2, 2))
    return(a + a)

if __name__ == '__main__':
    import timeit
    # print(timeit.timeit(costly_func))
    a = [1, 1, 1, 1]
    b = [1, 1]
    print(poly_long_div(a, b, q))

