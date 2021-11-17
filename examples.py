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


def ex1(n, x1var=x1):
    res = n % 2
    if n == 1:
        return x1var
    else:
        return a * ex1((n + res) // 2, x1var) + a * ex1((n - res) // 2, x1var) + c


def ex2(n, x1var=0):
    res = n % 2
    if n == 1:
        return x1var
    else:
        return ex2((n + res) // 2, x1var) + ex2((n - res) // 2, x1var) + n


def ex3(n, x1var=0):
    res = n % 2
    if n == 1:
        return x1var
    else:
        return ex3((n + res) // 2, x1var) + ex3((n - res) // 2, x1var) + (n-res) // 2


def ex4(n, x1var=1):
    res = n % 2
    if n == 1:
        return x1var
    else:
        return ex4((n + res) // 2, x1var) + ex4((n - res) // 2, x1var) + res


def ex5(n, x1var=1):
    res = n % 2
    if n == 1:
        return x1var
    else:
        return ex5((n + res) // 2, x1var) + ex5((n - res) // 2, x1var) + n ** 2 - res


def ex6(n, x1var=0):
    res = n % 2
    if n == 1:
        return x1var
    else:
        return - ex6((n + res) // 2, x1var) - ex6((n - res) // 2, x1var) + (n - res) // 2


def ex7(n, x1var=0):
    res = n % 2
    if n == 1:
        return x1var
    else:
        return ex7((n + res) // 2, x1var) + ex7((n - res) // 2, x1var) + binomial((n + res) // 2, 2) + binomial((n - res) // 2, 2)


def ex8(n, x1var=0):
    res = n % 2
    if n == 1:
        return x1var
    else:
        return ex8((n + res) // 2, x1var) + ex8((n - res) // 2, x1var) + binomial((n + res) // 2, 2) * binomial((n - res) // 2, 2)


def ex9(n, x1var=0):
    res = n % 2
    if n == 1:
        return x1var
    else:
        return 1/2*(ex9((n + res) // 2, x1var) + ex9((n - res) // 2, x1var) + res)


def ex10(i, n, x1var=0):
    res = n % 2
    if n == 0:
        return
    elif n == 1:
        return x1var
    else:
        g1 = (n-res)//2 * res
        g3 = 4 * (n-res)//2 + 2 * (n+res)//2 - 6 - g1
        if i == 1:
            return 2 * ex10(1, (n + res) // 2, x1var) + 2 * ex10(1, (n - res) // 2, x1var) + g1
        elif i == 3:
            return 2 * ex10(3, (n + res) // 2, x1var) + 2 * ex10(3, (n - res) // 2, x1var) + g3
        else:
            return "i not in {1,3}"
