# Define a struct to encapsulate a single-mode quantum optical system
struct SingleModeSystem
    N::Int                          # Dimension of Hilbert space
    A::Matrix{Float64}             # Annihilation operator
    Ad::Matrix{Float64}            # Creation operator
    psi0::Vector{Float64}          # Vacuum state |0⟩
    rho0::Matrix{Float64}          # Vacuum density matrix |0⟩⟨0|
end

# Constructor function to create the system
function SingleModeSystem(N::Int)
    A = zeros(N, N)
    for k in 1:(N-1)
        A[k, k+1] = sqrt(k)        # Standard annihilation operator in Fock basis
    end
    Ad = A'                        # Creation operator is Hermitian adjoint of A
    psi0 = zeros(N)
    psi0[1] = 1.0                  # |0⟩ state
    rho0 = psi0 * psi0'            # Density matrix ρ = |0⟩⟨0|
    return SingleModeSystem(N, A, Ad, psi0, rho0)
end


# Generate a Fock state density matrix |n⟩⟨n|
function fockn(sys::SingleModeSystem, n::Int)
    return sys.Ad^n * sys.rho0 * sys.A^n / factorial(n)
end

# Displacement operator D(α)
function D(sys::SingleModeSystem, alpha)
    return exp(alpha * sys.Ad .- conj(alpha) * sys.A)
end

# Squeezing operator S(α)
function Sq(sys::SingleModeSystem, alpha)
    return exp(0.5 * (alpha * sys.Ad * sys.Ad .- conj(alpha) * sys.A * sys.A))
end

# Kerr evolution operator exp(-iλ a†² a²)
function Kerr(sys::SingleModeSystem, lambda)
    return exp(-1im * lambda * sys.Ad * sys.Ad * sys.A * sys.A)
end

# Create a coherent state ρ = D(α) ρ₀ D†(α)
function coherent_state(sys::SingleModeSystem, alpha)
    D1 = D(sys, alpha)
    return D1 * sys.rho0 * D1'
end
