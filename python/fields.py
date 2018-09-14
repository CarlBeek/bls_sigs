class Fq(int):
    def __new__(cls, x: int, q: int):
        x = x % q
        ret = super().__new__(cls, x)
        ret.q = q
        return ret

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

    def inverse(self):
        t1 = self[0]


if __name__ == '__main__':
    a = Fq2(Fq(3, 17), Fq(4, 17))
    b = Fq2(Fq(5, 17), Fq(11, 17))
    print(a*b)