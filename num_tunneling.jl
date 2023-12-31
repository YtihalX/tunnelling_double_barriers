using LinearAlgebra
using Plots
using Interact

const a, b = 1, 1

function potential(x, V_0)
    if x<=2a || x>=2a+b && x <= 3a + b || x >= 3a + 2b
        return 0
    else
        return V_0
    end
end

function schrodinger_1d(N, L, V_0)
    dx = L / (N + 1)
    x = LinRange(0, L, N)
    diag_main = 2*ones(N)
    diag_off = -1*ones(N - 1)
    T = (1 / dx^2) * SymTridiagonal(diag_main, diag_off)
    V = Diagonal(potential.(x, V_0))
    return T + V
end

plt = plot()

# @manipulate for V_0 in 10:0.01:100
#     N = 1000
#     L = 4
#     H = schrodinger_1d(N, L, V_0)
#     energes, psi = eigen(H)
#     x = [i*L/(N+1) for i in 1:N]
#     plot(x, psi[:, 1:4])
# end

N = 5000
L = 8
H = schrodinger_1d(N, L, 100)
E = 80
energies, psi = eigen(H)
x = LinRange(L/(N+1), L, N)
plot(x, psi[:, 1])
