module AxialWoodsSaxonWithSphericalBasis

using Plots
using LinearAlgebra
using Parameters
using Base.Threads
using Memoize

include("./basis.jl")
include("./angular_momentum.jl")
include("./single_particle_states.jl")
include("./three_body.jl")
include("./wigner_matrix.jl")

export PhysicalParam

@with_kw struct PhysicalParam{T} @deftype Float64
    ħc = 197.
    mc² = 938.
    M = ħc^2/2mc²
    
    # particle number
    Z::Int64 = 8; @assert iseven(Z) === true
    N::Int64 = Z; @assert iseven(N) === true
    A::Int64 = Z + N; @assert A === Z + N
    
    # parameters of Woods-Saxon potential
    V₀ = -42.86 # [MeV]
    r₀ = 1.27 # [fm]
    R₀ = r₀*A^(1/3) # [fm]
    a = 0.67 # [fm]
    κ = 0.44
    
    # model space
    Emax = 30 # [MeV]
    lmax::Int64 = 10 
    Λmax::Int64 = 2lmax+1; @assert isodd(Λmax)
    
    # radial mesh
    Nr::Int64 = 100
    Δr = 0.2
    rs::T = range(Δr, Nr*Δr, length=Nr)
end







end # module
