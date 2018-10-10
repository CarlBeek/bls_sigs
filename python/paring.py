import params
from fields import Fq, Fq2, Fq6, Fq12
from ec import EC


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


def millers_alg(P, Q):
    def line(P_from, P_to, P_eval):
        if P_from == P_to:
            P_from = untwist(P_from)
            grad = (3 * P_from.x.square())/(2 * P_from.y)
            offset = P_from.y - grad * P_from.x
        else:
            P_from = untwist(P_from)
            P_to = untwist(P_to)
            grad = (P_to.y - P_from.y)/(P_to.x - P_from.x)
            offset = (P_to.y * P_from.x - P_from.y * P_to.x) * (P_from.x - P_to.y)
        return (- grad * P_eval.x) + P_eval.y + offset

    R = Q
    f = Fq12.one(params.q)

    for r_i in params.BLS_x[1:]:
        l_RR = line(R, R, P)
        R *= 2
        f = f.square() * l_RR

        if r_i == '1':
            l_RQ = line(R, Q, P)
            R += Q
            f *= l_RQ

    if params.BLS_negative:
        f = -f
    return f


def final_exp(R):
    #Todo: make this even remotely efficient...
    return R**((params.q**12-1)//params.r)
