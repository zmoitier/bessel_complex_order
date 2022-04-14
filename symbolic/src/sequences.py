""" Define sequences """

from sympy import Poly, Rational


def seq_lambda_mu(j_max):
    """Sequences u and v."""
    u = [Rational(1)]
    v = [Rational(1)]
    for j in range(1, j_max + 1):
        u.append(Rational((6 * j - 5) * (6 * j - 1) * u[j - 1], 48 * j))
        v.append(Rational((6 * j + 1) * u[j], 1 - 6 * j))
    return (u, v)


def seq_poly_UV(j_max, X):
    """Polinomial U and V."""
    U = [Poly(1, X, domain="QQ")]
    V = [Poly(1, X, domain="QQ")]

    D = Poly(X**2 * (1 - X**2) / 2, X)
    I = Poly((1 - 5 * X**2) / 8, X)

    D0 = Poly(X * (1 - X**2) / 2, X)
    D1 = Poly(X**2 * (1 - X**2), X)

    for j in range(1, j_max + 1):
        tmp = D * U[j - 1].diff() + (I * U[j - 1]).integrate()
        U.append(tmp)
        V.append(tmp - D0 * U[j - 1] - D1 * U[j - 1].diff())
    return (U, V)
