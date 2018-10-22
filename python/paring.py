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
    # Todo: Speed up by saving omega**2 and omega**3 as constatnt
    omega2 = Fq12(root, zero)
    omega3 = Fq12(zero, root)
    return EC.from_affine(omega2.inverse() * Fq12.to_cls(P.x, q), omega3.inverse() * Fq12.to_cls(P.y, q))


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
        grad = (_R.x.to_cls(3, _R.x.q) * _R.x.square()) / (_R.x.to_cls(2, _R.x.q) * _R.y)
        offset = _R.y - grad * _R.x
    else:
        _R = untwist(_R)
        _Q = untwist(_Q)
        if _R == - _Q:
            return _P.x - _R.x
        grad = (_Q.y - _R.y) / (_Q.x - _R.x)
        offset = (_Q.y * _R.x - _R.y * _Q.x) / (_R.x - _Q.x)
    # print(type(grad))
    # print(type(_P.x))
    _Px = grad.to_cls(_P.x, _P.x.q)
    _Py = grad.to_cls(_P.y, _P.y.q)
    return (- grad * _Px) + _Py - offset


def millers_alg(P, Q):
    R = Q
    f = Fq12.one(params.q)
    for r_i in bin(params.BLS_x)[3:]:
        l_RR = line(R, R, P)
        f = f.square() * l_RR
        R *= 2
        if r_i == '1':
            l_RQ = line(R, Q, P)
            f *= l_RQ
            R += Q
    if params.BLS_negative:
        f = -f
    return f


def final_exp(r):
    # Calculate exponentiation via "easy" and "hard" parts
    # https://eprint.iacr.org/2008/490.pdf

    # Easy part:
    m0 = r.frobenius_endo(6)
    m0 *= r.inverse()
    m = m0.frobenius_endo(2)
    m *= m0

    def exp_by_x(_f, _x):
        _f = _f ** _x
        if params.BLS_negative:
            _f = -_f
        return _f

    return m ** ((params.q**4 - params.q**2 + 1)//params.r)
    # Hard part:
    # Setup needed vars
    m_x = exp_by_x(m, params.BLS_x)
    m_x_2 = exp_by_x(m_x, params.BLS_x)
    m_x_3 = exp_by_x(m_x_2, params.BLS_x)
    m_p = m.frobenius_endo(1)
    m_p_2 = m.frobenius_endo(2)
    m_p_3 = m.frobenius_endo(3)
    m_x_p = m_x.frobenius_endo(1)
    m_x_2_p = m_x_2.frobenius_endo(1)
    m_x_3_p = m_x_3.frobenius_endo(1)
    m_x_2_p_2 = m_x_2.frobenius_endo(2)

    y0 = m_p * m_p_2 * m_p_3
    y1 = - m
    y2 = m_x_2_p_2
    y3 = - m_x_p
    y4 = - (m_x * m_x_2_p)
    y5 = - m_x_2
    y6 = - (m_x_3 * m_x_3_p)

    t0 = y6.square()
    t0 *= y4
    t0 *= y5
    t1 = y3 * y5
    t1 *= t0
    t0 *= y2
    t1 = t1.square()
    t1 *= t0
    t1 = t1.square()
    t0 = t1 * y1
    t1 *= y0
    t0 = t0.square()
    t0 *= t1
    return t0
