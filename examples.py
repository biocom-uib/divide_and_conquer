"""
All the examples presented here are those of the article, but they
are not computed in the most efficient way possible -- indeed, this
appendix is only meant to test the validity of the algorithms.
In some cases problems arise due to polynomial simplifications and
errors in the evaluation of the equality between `int` and `float`
instances.

For some methods, better algorithms already exist in the Python
repository https://github.com/biocom-uib/biotrees
"""


from sympy import symbols, binomial

a, c, x1 = symbols("a, c, x1")


def ex1(x1, n):
    res = n % 2
    if n == 1:
        return x1
    else:
        return a*ex1(x1, (n+res)//2) + a*ex1(x1, (n-res)//2) + c


def ex2(x1, n):
    res = n % 2
    if n == 1:
        return x1
    else:
        return ex2(x1, (n + res) // 2) + ex2(x1, (n - res) // 2) + n


def ex3(x1, n):
    res = n % 2
    if n == 1:
        return x1
    else:
        return ex3(x1, (n + res) // 2) + ex3(x1, (n - res) // 2) + (n-res) // 2


def ex4(x1, n):
    res = n % 2
    if n == 1:
        return x1
    else:
        return ex4(x1, (n + res) // 2) + ex4(x1, (n - res) // 2) + res


def ex5(x1, n):
    res = n % 2
    if n == 1:
        return x1
    else:
        return ex5(x1, (n + res) // 2) + ex5(x1, (n - res) // 2) + n**2 - res


def ex6(x1, n):
    res = n % 2
    if n == 1:
        return x1
    else:
        return - ex6(x1, (n+res)//2) - ex6(x1, (n-res)//2) + (n - res) // 2


def ex7(x1, n):
    res = n % 2
    if n == 1:
        return x1
    else:
        return ex7(x1, (n + res) // 2) + ex7(x1, (n - res) // 2) + binomial((n+res)//2, 2) + binomial((n-res)//2, 2)


def ex8(x1, n):
    res = n % 2
    if n == 1:
        return x1
    else:
        return ex8(x1, (n + res) // 2) + ex8(x1, (n - res) // 2) + binomial((n+res)//2, 2)*binomial((n-res)//2, 2)


def ex9(x1, n):
    res = n % 2
    if n == 1:
        return x1
    else:
        return 1/2*(ex9(x1, (n + res) // 2) + ex9(x1, (n - res) // 2) + res)


def ex10(i, x1, n):
    res = n % 2
    if n == 0:
        return
    elif n == 1:
        return x1
    else:
        g1 = (n-res)//2 * res
        g3 = 4 * (n-res)//2 + 2 * (n+res)//2 - 6 - g1
        if i == 1:
            return 2*ex10(1, x1, (n + res) // 2) + 2*ex10(1, x1, (n - res) // 2) + g1
        elif i == 3:
            return 2*ex10(3, x1, (n + res) // 2) + 2*ex10(3, x1, (n - res) // 2) + g3
        else:
            return "i not in {1,3}"
