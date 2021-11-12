from sympy import symbols, poly
from xnrt import xnrt
from utils import s, q, bs

x, y = symbols('x, y')
a, x1, n = symbols('a, x1, n')


def xn(avar, pol, x1var, nvar):
    p = poly(pol, x, y)
    b = bs(p)
    qsn = q(nvar)[s(nvar)]

    if nvar == 1:
        return x1var
    elif nvar > 1:
        return sum(b[r][t] * xnrt(nvar, r, t, avar) for r, t in p.monoms()) \
           + ((2 * avar) ** qsn + (2 * avar - 1) * (nvar * avar ** qsn - (2 * avar) ** qsn)) * x1var


def xnvar(avar, pol, x1var):
    return lambda nvar: xn(avar, pol, x1var, nvar)
