import params
from fields import Fq, Fq2, Fq6, Fq12
from ec import EC, TwistedEC


def paring(P, Q, exp=False):
    R = millers_alg(P, Q)
    if exp:
        R = final_exp(R)
    return R


def untwist(P):
    q = P.X.q
    root = Fq6(Fq2.zero(q), Fq2.one(q), Fq2.zero(q))
    zero = Fq6.zero(q)
    omega2 = Fq12(root, zero)
    omega3 = Fq12(zero, root)
    return EC.from_affine(omega2.inverse() * P.x, omega3.inverse() * P.y)


def twist(P):
    q = P.X.q
    root = Fq6(Fq2.zero(q), Fq2.one(q), Fq2.zero(q))
    zero = Fq6.zero(q)
    omega2 = Fq12(root, zero)
    omega3 = Fq12(zero, root)
    c0 = omega2 * P.x
    c1 = omega3 * P.y
    return TwistedEC.from_affine(c0.c0.c0, c1.c0.c0)


def line(_R, _Q, _P):
    if _R == _Q:
        _R = untwist(_R)
        grad = (3 * _R.x.square()) / (2 * _R.y)
        offset = _R.y - grad * _R.x
    else:
        _R = untwist(_R)
        _Q = untwist(_Q)
        if _R == - _Q:
            print(_P.x - _R.x)
            return _P.x - _R.x
        grad = (_Q.y - _R.y) / (_Q.x - _R.x)
        offset = (_Q.y * _R.x - _R.y * _Q.x) / (_R.x - _Q.x)
    return (- grad * _P.x) + _P.y - offset


def millers_alg(P, Q):
    R = Q
    f = Fq12.one(params.q)
    for r_i in bin(params.BLS_x )[3:]:
        l_RR = line(R, R, P)
        f = f.square() * l_RR
        R *= 2
        if r_i == '1':
            l_RQ = line(R, Q, P)
            f *= l_RQ
            R += Q
    # if params.BLS_negative:
    #     f = -f
    return f


def final_exp(r):
    return r ** ((params.q ** 12 - 1) // params.r)
    # f2 = r.inverse()
    # r = -r
    # r *= f2
    # r = r.frobenius_endo(2)
    # r *= f2
    #
    # def exp_by_x(_f, _x):
    #     _f = _f ** _x
    #     if params.BLS_negative:
    #         _f = -_f
    #     return _f
    #
    # x = params.BLS_x
    # y0 = r
    # y0 = y0.square()
    # y1 = exp_by_x(y0, x)
    # y2 = exp_by_x(y1, x >> 1)
    # y3 = -r
    # y1 *= y3
    # y1 = -y1
    # y1 *= y2
    # y2 = exp_by_x(y1, x)
    # y3 = exp_by_x(y2, x)
    # y3 *= -y1
    # y1 = y1.frobenius_endo(3)
    # y2 = y2.frobenius_endo(2)
    # y1 *= y2
    # y2 = exp_by_x(y3, x)
    # y2 *= y0
    # y2 *= r
    # y1 *= y2
    # y2 = y3.frobenius_endo(1)
    # y1 *= y2
    # return y1

