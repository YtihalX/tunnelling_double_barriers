using Symbolics
@variables A[1:4]
@variables B[1:4]
@variables k[1:3]
@variables x, V_0, a, b, m, E, hbar
psi1 = exp(sqrt(2*m*E)/hbar * x * im) + A[1] * exp(-sqrt(2*m*E)/hbar * x * im)
# psi1 = substitute(psi1, Dict(sqrt(2*m*E)/hbar => sqrt(2*m*E)/hbar), fold = false)

psi2 = B[1] * exp(sqrt(2*m*(V_0-E))/hbar * x * im) + B[2] * exp(-sqrt(2*m*(V_0-E))/hbar * x * im)
# psi2 = substitute(psi2, Dict(sqrt(2*m*(V_0-E))/hbar => sqrt(2*m*(V_0-E))/hbar), fold = false)

psi3 = A[2] * exp(sqrt(2*m*E)/hbar * x * im) + A[3] * exp(-sqrt(2*m*E)/hbar * x *im)
# psi3 = substitute(psi3, Dict(sqrt(2*m*E)/hbar => sqrt(2*m*E)/hbar), fold = false)

psi4 = B[3] * exp(sqrt(2*m*(V_0-E))/hbar * x * im) + B[4] * exp(-sqrt(2*m*(V_0-E))/hbar * x * im)
# psi4 = substitute(psi4, Dict(sqrt(2*m*(V_0-E))/hbar => sqrt(2*m*(V_0-E))/hbar), fold = false)

psi5 = A[4] * exp(sqrt(2*m*E)/hbar * x * im)
# psi5 = substitute(psi5, Dict(sqrt(2*m*E)/hbar => sqrt(2*m*E)/hbar), fold = false)
df = Differential(x)
sys_of_wave = [
    substitute(psi1, Dict(x => 0), fold = false) ~ substitute(psi2, Dict(x => 0), fold = false),
    substitute(df(psi1), Dict(x => 0), fold = false) ~ substitute(df(psi2), Dict(x => 0), fold = false),
    substitute(psi2, Dict(x => a), fold = false) ~ substitute(psi3, Dict(x => a), fold = false),
    substitute(df(psi2), Dict(x => a), fold = false) ~ substitute(df(psi3), Dict(x => a), fold = false),
    substitute(psi3, Dict(x => a + b), fold = false) ~ substitute(psi4, Dict(x => a + b), fold = false),
    substitute(df(psi3), Dict(x => a + b), fold = false) ~ substitute(df(psi4), Dict(x => a + b), fold = false),
    substitute(psi4, Dict(x => 2*a + b), fold = false) ~ substitute(psi5, Dict(x => 2*a + b), fold = false)
    substitute(df(psi4), Dict(x => 2*a + b), fold = false) ~ substitute(df(psi5), Dict(x => 2*a + b), fold = false)
]
Symbolics.solve_for([sys_of_wave], [A, B])
