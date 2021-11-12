from sympy import log, symbols


x, y = symbols('x, y')
a, x1, n = symbols('a, x1, n')


def bs(pol):
    degx = pol.degree(x)
    degy = pol.degree(y)
    rts = pol.monoms()

    bs = [[0 for _ in range(degy+1)] for _ in range(degx+1)]
    for r, t in rts:
        bs[r][t] = pol.coeff_monomial(x**r * y**t)
    return bs


def M(n, i):
    if i == 0:
        return 0
    else:
        return sum(2**qj for qj in q(n)[i:])


def delta(p):
    if p:
        return 1
    else:
        return 0


def q(n):
    bn = format(n, 'b')[::-1]
    lg = len(bn)
    return [0] + [i for i in range(lg) if bn[i] == '1']


def s(n):
    return int(sum(int(i) for i in format(n, 'b')))


def lg(a, t):
    return int(log(a, 2)) - t
