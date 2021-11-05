from sympy import symbols, poly, bernoulli, binomial, log

import sys
sys.setrecursionlimit(2000)

x, y = symbols('x, y')


def T(d, n, x): # I use the notations from the article // alerta la X
    if d == 0 and x == 1:
        return n
    elif d == 0:
        return (x**n - 1) / (x - 1)
    elif d == 1 and x == 1:
        return (n * (n-1)) / 2
    elif d == 1:
        return (n * x**n * (x - 1) - x * (x**n - 1)) / (x - 1)**2
    else:
        return sum(k**d * x**k for k in range(n))


def M(n, i):
    return sum(2**qj for qj in q(n)[i-1:])


def delta(p):
    if p:
        return 1
    else:
        return 0


def q(n):
    bn = format(n, 'b')[::-1]
    lg = len(bn)
    return [i for i in range(lg) if bn[i] == '1']


def s(n):
    return int(sum(int(i) for i in format(n, 'b')))


def alpha(n, d, m, a):  # q[i-1]
    return 1 / (2*(m+1)) * sum(sum(binomial(m+1, j) * bernoulli(j) * 2**j * M(n, i)**(m+1-j)
                                   * (T(d, q(n)[i-1], a*2**(j-m)) - T(d, q(n)[i-2], a*2**(j-m)))
                                   for j in range(m+1)) for i in range(1, s(n) + 1)) \
           + sum(q(n)[i]**d * (a / 2**(-m))**q(n)[i-1] * (n - M(n, i)) * M(n, i+1)**m
                 for i in range(1, s(n))) - T(d, q(n)[s(n) - 1], 2 * a) * delta(m == 0)


def bs(pol):     # de momento, asumo que la entrada *es* un polinomio de sympy
    degx = pol.degree(x)
    degy = pol.degree(y)
    rts = pol.monoms()
    bs = [[0 for _ in range(degy)] for _ in range(degx)]
    for r, t in rts:
        bs[r-1][t-1] = pol.coeff_monomial(x**r * y**t)
    return bs


def lg(a, t):
    return int(log(a, 2)) - t


def xnrt(a, n, r, t): # faltan los CASES
    qsn = q(n)[s(n)-1]
    xnrt1 = sum(sum(binomial(r, i-t-1) * binomial(i, k) * bernoulli(i-k) / (i*(2**(i-1)-a))
                   for i in range(k, r+t+1) if i != t + lg(a, t) + 1) * n**k for k in range(1, r+t+1))
    xnrt2 = 1 / (a - 1) * delta(r > 0 and t == 0 and a != 1)
    xnrt3 = (1 - sum(binomial(r, l) / (2**(t+l) - a) for l in range(r) if l != lg(a, t))) \
           * (T(0, q(n)[s(n)-1], 2*a) + n*a**(q(n)[s(n)-1]) - (2*a)**(q(n)[s(n)-1]))
    xnrt4 = sum(2**(-i) * binomial(r+t, i) - 2**(-i+1) * binomial(r, i-t)
                 - sum(binomial(r, l) * binomial(t+l, i) / (2**(t+l) - a) for l in range(i-t+1, r)) * alpha(n, 0, i, a)
                 for i in range(r+t))
    xnrt5 = delta(r > 0 and lg(a, t) in range(r)) / a * binomial(r, lg(a, t)) * (T(1,q(n)[s(n)-1], 2*a)
                                                                                + (n*a**(q(n)[s(n)-1]))*q(n)[s(n)-1]
                                                                                + sum(binomial(t+lg(a, t), i)
                                                                                      * alpha(n, 1, i, a)
                                                                                      for i in range(t+lg(a, t))))
    return xnrt1 + xnrt2 + xnrt3 + xnrt4 + xnrt5


def xn0t(a, n, t):
    qsn = q(n)[s(n)-1]
    if a == 1/2:
        return sum(2**(-i) * binomial(t, i) * alpha(n, 0, i, a) for i in range(t)) + qsn \
               + n * 2**(-qsn - 1)
    else:
        return sum(2**(-i) * binomial(t, i) * alpha(n, 0, i, a) for i in range(t)) + ((2*a)**qsn - 1) / (2*a - 1) \
               + n * a**qsn - (2*a)**qsn


def xnr0(a, n, r):
    qsn = q(n)[s(n)-1]
    if a == 1/2:
        sum1 = 2 * sum((sum(binomial(r, j-1) * binomial(j, k) * bernoulli(j-k) / (j * (2**j - 1))
                            for j in range(k, r+1))) * n**k for k in range(1, r+1))
        sum2 = (1 - 2 * sum(binomial(r,l) / (2**(l+1) - 1) for l in range(r))) * (qsn + n * 2**(-qsn)
                                                                                  - 1)
        sum3 = sum(2**(-i) * binomial(r,i) + 2 * sum(binomial(r,l) * binomial(l,i) / (2**(l+1) - 1) for l in range(i+1,r))
                   * alpha(n, 0, i, a) for i in range(r)) + 2
        return sum1 + sum2 - sum3
    elif a == 1:
        sum1 = sum(sum(binomial(r, j-1) * binomial(j, k) * bernoulli(j-k) / (j * (2**(j-1)-1)) for j in range(k, r+1))
                   * n**k for k in range(2, r+1))
        sum2 = (qsn + 1 + sum(binomial(r, j) * (bernoulli(j) - 1) / (2**j - 1) for j in range(1, r))) * n
        sum3 = 1 + sum(binomial(r, j) / (2**j - 1) for j in range(1, r)) - 2**(qsn + 1)
        sum4 = sum((2**(-i) * binomial(r, i) + sum(binomial(r, l) * binomial(l, i) / (2**l - 1) for l in range(i+1, r)))
                   * alpha(n, 0, i, a) for i in range(r))
        return sum1 + sum2 + sum3 - sum4
    elif any(a == 2**ell for ell in range(1, r)):
        ell = int(log(a, 2))
        sum1 = (1 - sum(binomial(r, l) / (2**l - a) for l in range(r) if l != ell)) * (((2*a)**qsn - 1) / (2*a - 1)
                                                                                       + n * a**qsn - (2*a)**qsn)
        sum2 = sum(sum(binomial(r, i-1) * binomial(i, k) * bernoulli(i-k) / (i*(2**(i-1) - a))
                       for i in range(k, r+1) if i != ell+1) * n**k for k in range(1, r+1)) + 1 / (a-1)
        sum3 = sum((2**(-i) * binomial(r, i) + sum(binomial(r, l) * binomial(l, i) / (2**l - a)
                                                   for l in range(i+1, r) if l != ell)) * alpha(n, 0, i, a) for i in range(r))
        sum4 = binomial(r, ell) / a * (((qsn - 1) * (2*a)**(qsn + 1) - qsn*(2*a)**qsn + 2*a) / (2*a - 1)**2
                                       + n*qsn*a**qsn + qsn*(2*a)**qsn)
        sum5 = binomial(r, ell) / a * sum(binomial(ell, i) * alpha(n, 1, i, a) for i in range(ell))
        return sum1 + sum2 - sum3 + sum4 + sum5
    else:
        sum1 = sum((sum(binomial(r, j-1) * binomial(j, k) * bernoulli(j-k) for j in range(k, r+1))) * n**k
                   for k in range(1, r+1))
        sum2 = (1 - sum(binomial(r, l) / (2**l - a) for l in range(r))) * (((2*a)**qsn - 1) / (2*a - 1) + n*a**qsn
                                                                           - (2*a)**qsn) + 1 / (a - 1)
        sum3 = sum((2**(-i) * binomial(r, i) - 2**(-i+1) * binomial(r, i) - sum(binomial(r, l) * binomial(l, i) / (2**l - a)
                                                                                for l in range(i+1, r))) * alpha(n, 0, i, a)
                   for i in range(r))
        return sum1 + sum2 + sum3


def xn(a, pol, x1, n):   # q[i-1] // AQUÍ LA ENTRADA NO ES POL DE SYMPY
    p = poly(pol)
    b = bs(p)
    qsn = q(n)[s(n)-1]
    return sum(b[r][t] * xnrt(a, n, r, t) for r, t in p.monoms()) \
           + ((2*a) ** q(n)[s(n) - 1] + (2*a - 1) * (n*a**qsn - (2*a) ** q(n)[s(n) - 1])) * x1