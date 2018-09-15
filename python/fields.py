class Fq(int):
    def __new__(cls, x: int, q: int):
        x = x % q
        ret = super().__new__(cls, x)
        ret.q = q
        return ret

    def __str__(self):
        return super().__str__()

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
        return Fq(self.inverse() * other, self.q)

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

    def square(self):
        return self * self


class Fq2(tuple):
    def __new__(cls, c0: Fq, c1: Fq):
        ret = super().__new__(cls, (c0, c1))
        ret.q = c0.q
        return ret

    def __add__(self, other):
        return Fq2(self[0] + other[0], self[1] + other[1])

    def __radd__(self, other):
        return self.__add__(other)

    def __neg__(self):
        return Fq2(self[0].__neg__(), self[1].__neg__())

    def __sub__(self, other):
        return Fq2(self[0]-other[0], self[1] - other[1])

    def __rsub__(self, other):
        return Fq2(other[0] - self[0], other[1] - self[1])

    def __mul__(self, other):
        aa = self[0] * other[0]
        bb = self[0] * other[1]
        o = other[0] + other[1]
        c1 = self[1] + self[0]
        c1 *= o
        c1 -= aa
        c1 -= bb
        c0 = aa - bb
        return Fq2(c0, c1)

    def __str__(self):
        return str(self[1]) + '*u + ' + str(self[0])

    def mul_by_nonresidue(self):
        return Fq2(self[0] - self[1], self[0] + self[1])

    def inverse(self):
        raise NotImplementedError


class Fq6(tuple):
    def __new__(cls, c0, c1, c2):
        ret = super().__new__(cls, (c0, c1, c2))
        ret.q = c0.q
        return ret

    def __neg__(self):
        c0 = -self[0]
        c1 = -self[1]
        c2 = -self[2]
        return Fq6(c0, c1, c2)

    def __add__(self, other):
        c0 = self[0] + other[0]
        c1 = self[1] + other[1]
        c2 = self[2] + other[2]
        return Fq6(c0, c1, c2)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        c0 = self[0] - other[0]
        c1 = self[1] - other[1]
        c2 = self[2] - other[2]
        return Fq6(c0, c1, c2)

    def __rsub__(self, other):
        c0 = other[0] - self[0]
        c1 = other[1] - self[1]
        c2 = other[2] - self[2]
        return Fq6(c0, c1, c2)

    def __str__(self):
        a = [element for tupl in self for element in tupl]
        c = 5
        a = list(map(lambda x: str(x) + 'u^' + str(c), a[::-1]))
        return str(a)


if __name__ == '__main__':
    a = Fq2(Fq(3, 17), Fq(4, 17))
    b = Fq2(Fq(5, 17), Fq(11, 17))
    c = Fq2(Fq(8, 17), Fq(14, 17))

    d = Fq6(a, b, c)
    print(d)
