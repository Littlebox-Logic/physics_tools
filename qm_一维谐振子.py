#!/usr/bin/env python3
# Logic Sun 2026.09.19

from sympy import *
from sympy.physics.quantum.constants import hbar

x, t = symbols("x t", real = True)
m, n, omega = symbols("m n omega", positive = True)

psi_0 = (m * omega / (pi * hbar)) ** Rational(1, 4) \
        * exp(- m * omega * x ** 2 / (2 * hbar))

V  = lambda x   : m * (omega * x) ** 2 / 2
p  = lambda psi : - I * hbar * diff(psi, x)
H  = lambda psi : simplify(- hbar ** 2 / (2 * m) * diff(psi, x, 2) + V(x) * psi)

ad = lambda psi : 1 / sqrt(2 * hbar * m * omega) * (- I * p(psi) + m * omega * x * psi)
a  = lambda psi : 1 / sqrt(2 * hbar * m * omega) * (+ I * p(psi) + m * omega * x * psi)

norm = lambda psi : simplify(psi / sqrt(integrate(Abs(psi) ** 2, (x, -oo, +oo))))

x_aver  = lambda psi : simplify(integrate(Abs(psi) ** 2 * x, (x, -oo, +oo)))
x2_aver = lambda psi : simplify(integrate((x * Abs(psi)) ** 2, (x, -oo, +oo)))
p_aver  = lambda psi : simplify(integrate(conjugate(psi) * p(psi), (x, -oo, +oo)))
p2_aver = lambda psi : simplify(integrate(conjugate(psi) * p(p(psi)), (x, -oo, +oo)))

sigma_x = lambda psi : sqrt(x2_aver(psi) - x_aver(psi) ** 2)
sigma_p = lambda psi : sqrt(p2_aver(psi) - p_aver(psi) ** 2)

def state(n : int):
    psi = psi_0
    while (n := n - 1) + 1:
        psi = ad(psi)
    return norm(psi)

def time(psi, *elem):
    Psi = 0
    psi_n, last_n = psi_0, 0
    for n in sorted(elem):
        while (last_n < n):
            psi_n = ad(psi_n)
            last_n += 1
        Psi += integrate(conjugate(norm(psi_n)) * psi, (x, -oo, +oo)) * psi_n * exp(- I * omega * t * (n + Rational(1, 2)))
    return Psi

show = lambda expr : preview(simplify(expr), dvioptions = ['-D', '180'], euler = False)

if __name__ == "__main__":
    psi = psi_0
    formula = "\\begin{aligned}"
    for n in range(3):
        psi = norm(ad(psi)) if n else psi_0
        formula += f"\\psi_{n}(x) &= {latex(psi)}\\\\" \
                   f"E_{n} &= {latex(H(psi) / psi)}\\\\" \
                   f"\\langle\\psi_{n}|x|\\psi_{n}\\rangle &= {latex(x_aver(psi))}\\\\" \
                   f"\\langle\\psi_{n}|x^2|\\psi_{n}\\rangle &= {latex(x2_aver(psi))}\\\\" \
                   f"\\langle\\psi_{n}|\\hat{{p}}|\\psi_{n}\\rangle &= {latex(p_aver(psi))}\\\\" \
                   f"\\langle\\psi_{n}|\\hat{{p}}^2|\\psi_{n}\\rangle &= {latex(p2_aver(psi))}\\\\" \
                   f"\\sigma_{{x{n}}}\\cdot\\sigma_{{p{n}}} &= {latex(sigma_x(psi) * sigma_p(psi))}\\\\"
    formula += "\\end{aligned}"

    print(formula)
    preview(f"$${formula}$$", dvioptions = ['-D', '180'], euler = False)
