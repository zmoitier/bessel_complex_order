""" Helper for printing """

from sympy import Poly


def num2str(num, prec=18):
    """num2str"""
    if num == 0:
        return 0

    return f"{num.evalf(prec):.18e}"


def eval_coeffs(expr, var):
    """Eval coefficient of series."""
    if expr.is_number:
        print(num2str(expr))

    elif isinstance(expr, Poly):
        for c in expr.all_coeffs():
            print(num2str(c))

    else:
        for k in reversed(range(expr.getn())):
            print(num2str(expr.coeff(var, k)))
