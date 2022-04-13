""" Helper for printing """


def num2str(num, prec=18):
    """num2str"""
    if num == 0:
        return 0

    return f"{num.evalf(prec):.18e}"


def eval_coeffs(expr, var):
    """Eval coefficient of series."""
    if expr.is_number:
        print([num2str(expr)])
    else:
        print([num2str(expr.coeff(var, k)) for k in range(expr.getn())])
