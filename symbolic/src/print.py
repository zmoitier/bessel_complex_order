""" Helper for printing """


def num2str(num, prec=18):
    """num2str"""
    if num == 0:
        return 0

    return f"{num.evalf(prec):.18e}"
