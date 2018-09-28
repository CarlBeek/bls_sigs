import params
from fields import Fq, Fq2


class g1Affine:
    def __init__(self, x, y, infinity=False):
        self.basefield = x.zero()
        self.x = x
        self.y = y

        self.infinity = infinity
        # This method sets x and y = 0 in the case of infinity
        self.is_infinity()

    def __str__(self):
        if self.is_infinity():
            return 'Infinity'
        return 'x = ' + str(self.x) + '\ny = ' + str(self.y)

    def is_infinity(self):
        if self.infinity:
            self.x = self.basefield.zero()
            self.y = self.basefield.zero()
        return self.infinity

    @classmethod
    def get_point_from_x(cls, x, greatest=False):
        x3b = x**3 + params.b
        y = x3b.sqrt()
        if greatest and y < - y:
            y = - y
        return cls(x, y)

    @classmethod
    def from_projective(cls, proj):
        if proj.infinity():
            return cls()
        x = proj.X/(proj.Z**2)
        y = proj.y/(proj.Z**3)
        return cls(x, y)


class g1Projective:
    def __init__(self, X, Y, Z):
        self.basefield = X.zero()
        self.X = X
        self.Y = Y
        self.Z = Z

        # This method sets X and Y = 0 in the case of infinity
        self.is_infinity()

    def __str__(self):
        if self.is_infinity():
            return 'Infinity'
        return 'X = ' + str(self.X) + '\nY = ' + str(self.Y) + '\nZ = ' + str(self.Z)

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

    @classmethod
    def from_affine(cls, affine):
        if affine.is_infinity():
            cls(affine.x, affine.y, affine.basefield.zero())
        return cls(affine.x, affine.y, affine.basefield.one())


if __name__ == '__main__':
    a = g1Affine.get_point_from_x(Fq(3, 11), True)
    print(g1Projective.from_affine(a))
