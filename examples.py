"""
All the examples presented here are those of the article, but they
are not computed in the most efficient way possible -- indeed, this
appendix is only meant to test the validity of the algorithms.
In some cases problems arise due to polynomial simplifications and
errors in the evaluation of the equality between `int` and `float`
instances.

In all instances, `x1` is assumed to be the same as in the article,
so in order to invoke the method it is as easy as writing, say,
`ex1(3)` without having to specify the initial condition.

For some methods, better algorithms already exist in the Python
repository https://github.com/biocom-uib/biotrees
"""


from sympy import symbols, binomial

a, c, x1 = symbols("a, c, x1")


def ex1(n, x1var=x1):
    """
    A divide-and-conquer equation with a constant independent term.
    Compare it with Xn = xnvar(a, c, x1), and then in the case a=1
    with `Yn = xnvar(1, c, x1)`, as in `ex1(4) == Yn(4)`. Notice that
    since the result in both     cases is a polynomial, they might
    not be considered to be the same at first sight and might need
    to be expanded.
    :param n: `int` instance
    :param x1var: `SymPy` expression, `int` instance
    :return: `SymPy` expression, `int` instance
    """
    res = n % 2
    if n == 1:
        return x1var
    else:
        return a * ex1((n + res) // 2, x1var) + a * ex1((n - res) // 2, x1var) + c


def ex2(n, x1var=0):
    """
    Minimum total Sackin index Sn of a rooted bifurcating tree with
    `n` leaves.
    Compare it with `Sn = xnvar(1, "x+y", 0)`, as in `ex2(4) == Sn(4)`.
    :param n: `int` instance
    :param x1var: `SymPy` expression, `int` instance
    :return: `SymPy` expression, `int` instance
    """
    res = n % 2
    if n == 1:
        return x1var
    else:
        return ex2((n + res) // 2, x1var) + ex2((n - res) // 2, x1var) + n


def ex3(n, x1var=0):
    """
    Minimum Colless index Sn of a rooted bifurcating tree with `n` leaves.
    Compare it with `Cn = xnvar(1, "x-y", 0)`, as in `ex3(4) == Cn(4)`.
    :param n: `int` instance
    :param x1var: `SymPy` expression, `int` instance
    :return: `SymPy` expression, `int` instance
    """
    res = n % 2
    if n == 1:
        return x1var
    else:
        return ex3((n + res) // 2, x1var) + ex3((n - res) // 2, x1var) + (n-res) // 2


def ex4(n, x1var=1):
    """
    The parabola `n**2`.
    Compare it with `Xn = xnvar(2, "y-x", 1)`, as in `ex4(4) == Xn(4)`.
    :param n: `int` instance
    :param x1var: `SymPy` expression, `int` instance
    :return: `SymPy` expression, `int` instance
    """
    res = n % 2
    if n == 1:
        return x1var
    else:
        return 2*ex4((n + res) // 2, x1var) + 2*ex4((n - res) // 2, x1var) - res


def ex5(n, x1var=1):
    """
    The sequence of Lebesgue constants of the Walsh system.
    Compare it with `Ln = xnvar(1/2, "1/2*(x-y)", 1)`, as in `ex5(4) == Ln(4)`.
    :param n: `int` instance
    :param x1var: `SymPy` expression, `int` instance
    :return: `SymPy` expression, `int` instance
    """
    res = n % 2
    if n == 1:
        return x1var
    else:
        return 1/2 * ex5((n + res) // 2, x1var) + 1/2 * ex5((n - res) // 2, x1var) + 1/2 * res


def ex6(n, x1var=0):
    """
    Sequence A005536 in the OEIS.
    Compare it with `Xn = xnvar(-1, "y", 0)`, as in `ex6(4) == Xn(4)`.
    :param n: `int` instance
    :param x1var: `SymPy` expression, `int` instance
    :return: `SymPy` expression, `int` instance
    """
    res = n % 2
    if n == 1:
        return x1var
    else:
        return - ex6((n + res) // 2, x1var) - ex6((n - res) // 2, x1var) + (n - res) // 2


def ex7(n, x1var=0):
    """
    Sequence A087733 in the OEIS.
    Compare it with `Xn = xnvar(-1, "x*y", 0)`, as in `ex7(4) == Xn(4)`.
    :param n: `int` instance
    :param x1var: `SymPy` expression, `int` instance
    :return: `SymPy` expression, `int` instance
    """
    res = n % 2
    if n == 1:
        return x1var
    else:
        return - ex7((n + res) // 2, x1var) - ex7((n - res) // 2, x1var) + (n - res) * (n + res) // 4


def ex8(n, x1var=0):
    """
    Minimum total cophenetic index of a rooted bifurcating tree with `n` leaves.
    Compare it with `Fn = xnvar(1, "1/2*(x**2 + y**2 - x - y)", 0)`, as in `ex8(4) == Fn(4)`.
    :param n: `int` instance
    :param x1var: `SymPy` expression, `int` instance
    :return: `SymPy` expression, `int` instance
    """
    res = n % 2
    if n == 1:
        return x1var
    else:
        return ex8((n + res) // 2, x1var) + ex8((n - res) // 2, x1var) + binomial((n + res) // 2, 2) + binomial((n - res) // 2, 2)


def ex9(n, x1var=0):
    """
    Minimum rooted quartet index index of a rooted bifurcating tree with `n` leaves.
    Compare it with `QIn = xnvar(1, "1/4*(x*y - x**2 * y - x * y**2  + x**2 * y**2)", 0)`,
    as in `ex9(4) == QIn(4)`.
    :param n: `int` instance
    :param x1var: `SymPy` expression, `int` instance
    :return: `SymPy` expression, `int` instance
    """
    res = n % 2
    if n == 1:
        return x1var
    else:
        return ex9((n + res) // 2, x1var) + ex9((n - res) // 2, x1var) + binomial((n + res) // 2, 2) * binomial((n - res) // 2, 2)


def ex10(i, n, x1var=0):
    """
    Sequences A006581 and A006583 in the OEIS, depending whether `i == 1` or `i == 3`.
    They depend on a paramete `i` that can be 1 or 3. Compare them with
    * if i == 1, `An = xnvar(2, "x*y - y**2", 0)`, as in `ex10(1, 4) == An(4)`;
    * if i == 3, `Bn = xnvar(2, "4*y + 2*x - 6 - (x*y - y**2)", 0)`, as in `ex10(3, 4) == Bn(4)`.
    :param i: `int` instance
    :param n: `int` instance
    :param x1var: `SymPy` expression, `int` instance
    :return: `SymPy` expression, `int` instance
    """
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
