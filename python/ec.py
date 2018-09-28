import params
from fields import Fq, Fq2


class g1:
    def __init__(self, X, Y, Z):
        self.X = X
        self.Y = Y
        self.Z = Z

        # This method sets X and Y = 0 in the case of infinity
        self.is_infinity()

    def __str__(self):
        if self.is_infinity():
            return 'Infinity'
        return 'X = ' + str(self.X) + ', Y = ' + str(self.Y) + ', Z = ' + str(self.Z)

    def is_infinity(self):
        if self.Z.is_zero():
            self.X = self.basefield.zero()
            self.Y = self.basefield.zero()
            return True
        return False

    def is_on_curve(self):
        if self.is_infinity():
            return True
        return self.Y**2 == self.X**3 + params.b*(self.Z**6)

    def as_affine(self):
        return self.x, self.y, self.is_infinity()

    @property
    def x(self):
        return self.X * (self.Z**2).inverse()

    @property
    def y(self):
        return self.Y * (self.Z**3).inverse()

    @property
    def basefield(self):
        return self.X.zero()

    @classmethod
    def from_affine(cls, x, y, infinity=False):
        if infinity:
            cls(x, y, x.zero())
        return cls(x, y, x.one())

    @classmethod
    def get_point_from_x(cls, X, greatest):
        X3b = X ** 3 + params.b
        Y = X3b.sqrt()
        if greatest and Y < - Y:
            Y = - Y
        return cls(X, Y, X.one())


if __name__ == '__main__':
    A = g1.get_point_from_x(Fq(3, 11), True)
    print(A.as_affine())
