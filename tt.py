from sympy import *
from sympy.physics.quantum.constants import hbar
A = symbols('A:4')
B = symbols('B:4')
k0, k1 = symbols('k:2')
V_0, a, b, m, E = symbols('V_0, a, b, m, E', positive = True)
x = symbols('x', real = True)
psi1 = exp(k0*x*I)+A[0]*exp(-k0*x*I)
psi2 = B[0] * exp(k1 * x) + B[1] * exp(-k1 * x)
psi3 = A[1] * exp(k0 * x * I) + A[2] * exp(-k0 * x *I)
psi4 = B[2] * exp(k1 * x) + B[3] * exp(-k1 * x)
psi5 = A[3] * exp(k0 * x * I)
eqs = (
    Eq(psi1.subs(x, 0), psi2.subs(x, 0)),
    Eq(psi1.diff(x).subs(x, 0), psi2.diff(x).subs(x, 0)),
    Eq(psi2.subs(x, a), psi3.subs(x, a)),
    Eq(psi2.diff(x).subs(x, a), psi3.diff(x).subs(x, a)),
    Eq(psi3.subs(x, a+b), psi4.subs(x, a+b)),
    Eq(psi3.diff(x).subs(x, a+b), psi4.diff(x).subs(x, a+b)),
    Eq(psi4.subs(x, 2*a+b), psi5.subs(x, 2*a+b)),
    Eq(psi4.diff(x).subs(x, 2*a+b), psi5.diff(x).subs(x, 2*a+b))
)
vars = A+B
# sol = solve(eqs, A+B)
