# divide_and_conquer

This Python repository is a companion to the article *Explicit solution of a divide-and-conquer dividing by a half recurrence with polynomial independent term*, by Tomás M. Coronado, Arnau Mir and Francesc Rosselló.

The main file of the repository is `divide_and_conquer`; the files `utils` and `xnrt` are auxiliar to it, while `examples` gives the recurrent definition of the examples found in the article.

In order to use the program in `divide_and_conquer`, the user needs to have a working version of **SymPy** (https://www.sympy.org/en/index.html). 

In the main file we find two programs, `xn` and `xnvar`:
* `xn` one gives the value, computed by means of the formulas in the aforementioned article, of a specific value of a specific recurrence *x_n*,
* `xnvar` allows us to define a variable, say `Xn`, as `Xn = xnvar(1, "x+y", 0)`, so that we can then compute `Xn(1)`, `Xn(2)`, `Xn(3)` and so on...

The polynomial variable deserves a little explanation. In the paper, we deal with polynomials over the variables *ceiling(n / 2)* and *floor(n / 2)*. In order to ease the task for the reader, we have set `x` to represent *ceiling(n / 2)* and `y` to represent *floor(n / 2)*. So that, for example, if we want to write *n* we need to write the polynomial `"x + y"` as a string.
