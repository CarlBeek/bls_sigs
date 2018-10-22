from paring import paring
from params import g1_x, g1_y, g2_x, g2_y, q, q_bits
from ec import EC, TwistedEC
from fields import Fq, Fq2

g1 = EC.from_affine(g1_x, g1_y)
g2 = TwistedEC.from_affine(g2_x, g2_y)


def sign(message, private_key):
    return message * private_key


def verify(message, signature, public_key):
    if not signature.is_on_curve():
        return False
    return paring(public_key, message, exp=True) == paring(g1, signature, exp=True)


def key_gen(private_key):
    return g1 * private_key


def compress(point):
    res = 0
    for i, x in enumerate(point.x):
        res += x << (i*q_bits)
    res += (point.y > - point.y) << (2 * q_bits)  # Encodes the parity bit for the Y co-ord
    return res


def decompress(point):
    greatest = point >> 2*q_bits
    x_c0 = Fq(point & (2 ** q_bits - 1), q)
    x_c1 = Fq((point >> q_bits) & (2 ** q_bits - 1), q)
    x = Fq2(x_c0, x_c1)
    return TwistedEC.get_point_from_x(x, greatest)
