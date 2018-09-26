import params
from fields import Fq, Fq2

class g1Affine:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.infinity = False # Todo: Check this rather than assuming not inf

    def __str__(self):
        if self.infinity:
            return 'Infinity'
        return 'x = ' + str(self.x) + '\ny = ' + str(self.y)

    @classmethod
    def get_point_from_x(cls, x, greatest=False):
        x3b = x**3 + params.b
        y = x3b.sqrt()
        if greatest and y < - y:
            y = - y
        return g1Affine(x, y)


if __name__ == '__main__':
    a = g1Affine.get_point_from_x(Fq(3, 11))
    print(a)
