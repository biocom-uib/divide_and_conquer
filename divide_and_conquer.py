"""
This Python repository is a companion to the article *Explicit solution of a divide-and-conquer
dividing by a half recurrence with polynomial independent term*, by Tomás M. Coronado, Arnau Mir
and Francesc Rosselló.

The polynomial variable deserves a little explanation. In the paper, we deal with polynomials
over the variables ceiling(n / 2) and floor(n / 2). In order to ease the task for the reader,
we have set x to represent ceiling(n / 2) and y to represent floor(n / 2). So that, for example,
if we want to write n we need to write the polynomial "x + y" as a string.
"""


from sympy import symbols, poly
from xnrt import xnrt
from utils import s, q, bs

x, y = symbols('x, y')
a, x1, n = symbols('a, x1, n')
"""
The user might choose to use variables in order to keep their computations abstract. In that case, the
method `.subs()` from the SymPy library must be used. 
"""


def xn(avar, pol, x1var, nvar):
    """
    Gives the value, computed by means of the formulas in the aforementioned article,
    of a specific value of a specific recurrence x_n, given:
    :param avar: `float` or `int` instance
    :param pol: `string` instance
    :param x1var: `float` or `int` instance
    :param nvar: `int` instance
    :return: `float` or `int` or SymPy object instance
    """
    p = poly(pol, x, y)
    b = bs(p)
    qsn = q(nvar)[s(nvar)]

    if nvar == 1:
        return x1var
    elif nvar > 1:
        return sum(b[r][t] * xnrt(nvar, r, t, avar) for r, t in p.monoms()) \
           + ((2 * avar) ** qsn + (2 * avar - 1) * (nvar * avar ** qsn - (2 * avar) ** qsn)) * x1var


def xnvar(avar, pol, x1var):
    """
    Gives a lambda expression that, given an `int`, returns the value of a recurrence in that integer, given:
    :param avar: `float` or `int` instance
    :param pol: `string` instance
    :param x1var: `float` or `int` instance
    :return: `lambda` expression from `int` instance to `float` or `int` or SymPy object instance
    """
    return lambda nvar: xn(avar, pol, x1var, nvar)
