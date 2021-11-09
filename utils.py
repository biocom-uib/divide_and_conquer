from sympy import log


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


def q(n):       # siguiendo la idea de Cesc de que q0 = 0
    bn = format(n, 'b')[::-1]
    lg = len(bn)
    return [0] + [i for i in range(lg) if bn[i] == '1']


def s(n):
    return int(sum(int(i) for i in format(n, 'b')))


def lg(a, t):
    return int(log(a, 2)) - t

