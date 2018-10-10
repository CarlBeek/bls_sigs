from paring import paring
from params import g1_x, g1_y, g2_x, g2_y
from ec import EC, TwistedEC

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
