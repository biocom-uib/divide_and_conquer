from utils import *
from sympy import binomial, bernoulli


def T(d, n, x):
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


def alpha(n, d, m, a):
    qsn = q(n)[s(n)]
    if m == 0:
        sum1 = 1 / 2 * sum(M(n, i) * (T(d, q(n)[i], a) - T(d, q(n)[i-1], a)) for i in range(1, s(n) + 1))
        sum2 = sum(q(n)[i]**d * a**q(n)[i] * (n - M(n, i)) for i in range(1, s(n))) - T(d, qsn, 2 * a)

    else:
        sum1 = 1 / (2*(m+1)) * sum(sum(binomial(m+1, j) * bernoulli(j) * 2**j * M(n, i)**(m+1-j)
                                       * (T(d, q(n)[i], a*2**(j-m)) - T(d, q(n)[i-1], a*2**(j-m)))
                                       for j in range(m+1)) for i in range(1, s(n) + 1))
        sum2 = sum(q(n)[i]**d * (a * 2**(-m))**q(n)[i] * (n - M(n, i)) * M(n, i+1)**m
                     for i in range(1, s(n)))
    return sum1 + sum2


def xnrt(n, r, t, a):
    if r == 0:
        return xn0t(n, t, a)
    elif t == 0:
        return xnr0(n, r, a)
    else:
        qsn = q(n)[s(n)]

        if a == 1/2:
            sum1 = 2 * sum(sum(binomial(r, j-t-1) * binomial(j, k) * bernoulli(j-k) / (j*(2**j - 1))
                               for j in range(k, r+t+1)) * n**k for k in range(1, r+t+1))
            sum2 = (1 - 2 * sum(binomial(r, l) / (2**(t+l+1) - 1) for l in range(r))) * (qsn + n*2**(-qsn) - 1)
            sum3 = sum((2**(-i) * binomial(r+t, i) - 2**(-i+1) * binomial(r, i-t) -
                       2 * sum(binomial(r, l) * binomial(t+l, i) / (2**(t+l+1) - 1) for l in range(i-t+1, r)))
                       * alpha(n, 0, i, a) for i in range(r+t))
            return sum1 + sum2 + sum3

        elif a == 1:
            sum1 = sum(sum(binomial(r, j-t-1) * binomial(j, k) * bernoulli(j-k) / (j * (2**(j-1) - 1))
                               for j in range(k, r+t+1)) * n**k for k in range(2, r+t+1))
            sum2 = sum(binomial(r, j-t-1) * bernoulli(j-1) / (2**(j-1) - 1)
                       for j in range(2, r+t+1)) * n
            sum3 = (1 - sum(binomial(r, l) / (2**(t+l) - 1) for l in range(r))) * (n - 1)
            sum4 = sum((2**(-i) * binomial(r+t, i) - 2**(-i+1) * binomial(r, i-t)
                       - sum(binomial(r, l) * binomial(t+l, i) / (2**(t+l) - 1) for l in range(i-t+1, r))) * alpha(n, 0, i, 1)
                       for i in range(r+t))
            return sum1 + sum2 + sum3 + sum4

        elif any(a == 2**(t+ell) for ell in range(r)):
            ell = int(log(a, 2)) - t
            sum1 = sum((sum(binomial(r, j-t-1) * binomial(j, k) / (j * (2**(j-1) - a)) for j in range(k, r+t+1)
                            if j != t + ell + 1)) * n**k for k in range(1, r+t+1))
            sum2 = binomial(r, ell) / a * (((qsn - 1) * (2*a)**(qsn + 1) - qsn * (2*a)**qsn + 2*a) / (2*a - 1)**2
                                           + n*qsn*a**qsn - qsn*(2*a)**qsn)
            sum3 = (1 - sum(binomial(r, l) / (2**(t+l) - a) for l in range(r) if l != ell)) \
                   * (((2*a)**qsn - 1) / (2*a - 1) + n*a**qsn - (2*a)**qsn)
            sum4 = sum((2**(-i) * binomial(r+t, i) - 2**(-i+1) * binomial(r, i-t) - sum(binomial(r, l)
                                                                                        * binomial(t+l, i) / (2**(t+l) - a)
                        for l in range(i-t+1, r) if l != ell)) * alpha(n, 0, i, a) for i in range(r+t-1))
            sum5 = binomial(r, ell) / a * sum(binomial(t + ell, i) * alpha(n, 1, i, a) for i in range(t + ell))
            return sum1 + sum2 + sum3 + sum4 + sum5

        else:
            sum1 = sum(sum(binomial(r, j-t-1) * binomial(j, k) * bernoulli(j-k) / (j * (2**(j-1) - a))
                               for j in range(k, r+t+1)) * n**k for k in range(1, r+t+1))
            sum2 = (1 - sum(binomial(r, l) / (2**(t+l) - a) for l in range(r))) * (((2*a)**qsn - 1) / (2*a - 1)
                                                                                   + n*a**qsn - (2*a)**qsn)
            sum3 = sum((2**(-i) * binomial(r+t, i) - 2**(-i+1) * binomial(r, i-t)
                       - sum(binomial(r, l) * binomial(t+l, i) / (2**(t+l) - a) for l in range(i-t+1, r))) * alpha(n, 0, i, a)
                       for i in range(r+t))
            return sum1 + sum2 + sum3


def xn0t(n, t, a):
    qsn = q(n)[s(n)]

    if a == 1/2:
        return sum(2**(-i) * binomial(t, i) * alpha(n, 0, i, a) for i in range(t)) + qsn \
               + n * 2**(-qsn) - 1

    else:
        return sum(2**(-i) * binomial(t, i) * alpha(n, 0, i, a) for i in range(t)) + ((2*a)**qsn - 1) / (2*a - 1) \
               + n * a**qsn - (2*a)**qsn


def xnr0(n, r, a):
    qsn = q(n)[s(n)]

    if a == 1/2:
        sum1 = 2 * sum((sum(binomial(r, j-1) * binomial(j, k) * bernoulli(j-k) / (j * (2**j - 1))
                            for j in range(k, r+1))) * n**k for k in range(1, r+1))
        sum2 = (1 - 2 * sum(binomial(r, l) / (2**(l+1) - 1) for l in range(r))) * (qsn + n * 2**(-qsn) - 1)
        sum3 = sum((2**(-i) * binomial(r, i) + 2 * sum(binomial(r, l) * binomial(l, i) / (2**(l+1) - 1) for l in range(i+1, r)))
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
        sum4 = binomial(r, ell) / a * ((qsn * (2*a)**(qsn+1) - qsn*(2*a) - (2*a)**(qsn+1) + 2*a) / (2*a - 1)**2
                                       + qsn * a**qsn * (n - 2**qsn))
        sum5 = binomial(r, ell) / a * sum(binomial(ell, i) * alpha(n, 1, i, a) for i in range(ell))
        return sum1 + sum2 - sum3 + sum4 + sum5

    else:
        sum1 = sum((sum(binomial(r, j-1) * binomial(j, k) * bernoulli(j-k) / (j * (2**(j-1) - a)) for j in range(k, r+1)))
                   * n**k for k in range(1, r+1))
        sum2 = (1 - sum(binomial(r, l) / (2**l - a) for l in range(r))) * (((2*a)**qsn - 1) / (2*a - 1) + n*a**qsn
                                                                           - (2*a)**qsn) + 1 / (a - 1)
        sum3 = sum((2**(-i) * binomial(r, i) - 2**(-i+1) * binomial(r, i) - sum(binomial(r, l) * binomial(l, i) / (2**l - a)
                                                                                for l in range(i+1, r))) * alpha(n, 0, i, a)
                   for i in range(r))
        return sum1 + sum2 + sum3
