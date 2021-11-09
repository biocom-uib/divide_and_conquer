from sympy import symbols, poly
from xnrt import xnrt
from utils import s, q

x, y = symbols('x, y')
a, x1, n = symbols('a, x1, n')


def bs(pol):     # de momento, asumo que la entrada *es* un polinomio de sympy
    degx = pol.degree(x)
    degy = pol.degree(y)
    rts = pol.monoms()

    bs = [[0 for _ in range(degy+1)] for _ in range(degx+1)]
    for r, t in rts:
        bs[r-1][t-1] = pol.coeff_monomial(x**r * y**t)
    return bs


def xn(avar, pol, x1var, nvar):   # q[i-1] // AQUÍ LA ENTRADA NO ES POL DE SYMPY
    p = poly(pol, x, y)
    b = bs(p)
    qsn = q(nvar)[s(nvar)]

    return sum(b[r][t] * xnrt(nvar, r, t, avar) for r, t in p.monoms()) \
           + ((2 * avar) ** qsn + (2 * avar - 1) * (nvar * avar ** qsn - (2 * avar) ** qsn)) * x1var
