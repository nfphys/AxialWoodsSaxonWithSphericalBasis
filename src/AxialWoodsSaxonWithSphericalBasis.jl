module AxialWoodsSaxonWithSphericalBasis

using Plots
using LinearAlgebra
using SparseArrays
using Parameters
using Base.Threads
using KrylovKit
using ProgressMeter
using Roots

include("./basis.jl")
include("./angular_momentum.jl")
include("./single_particle_states.jl")
include("./legendre.jl")
include("./wigner_matrix.jl")

include("./three_body_Hamiltonian.jl")
include("./three_body_density.jl")
include("./three_body.jl")
include("./greens_function.jl")

export PhysicalParam, calc_MEs_ang

@with_kw struct PhysicalParam{T} @deftype Float64
    ħc = 197.
    mc² = 938.
    
    # particle number of the core nucleus
    Z::Int64 = 6; @assert iseven(Z) 
    N::Int64 = 14; @assert iseven(N) 
    A::Int64 = Z + N; @assert A === Z + N

    M = ħc^2/2mc²*(1 + 1/A)
    
    # parameters of Woods-Saxon potential
    V₀ = -38.76 # [MeV]
    V₁ = -25.63/1.25^2
    r₀ = 1.25 # [fm]
    R₀ = r₀*A^(1/3) # [fm]
    a = 0.65 # [fm]
    V_gaus = 4.66
    μ_gaus = 0.09
    
    # model space
    Emax = 5 # [MeV]
    lmax::Int64 = 5
    Λmax::Int64 = 2lmax+1; @assert isodd(Λmax)

    # radial mesh
    Nr::Int64 = 60
    Δr = 0.5
    rs::T = range(Δr, Nr*Δr, length=Nr)
    ir_matching::Int64 = floor(Int, R₀/Δr)
    Jmax::Int64 = 6

    # parameters of neutron-neutron interaction 
    a_nn = -15.0 # [fm]
    v₀_nn = 2π^2*ħc^2/mc² * 2a_nn/(π - 2a_nn*sqrt(mc²*Emax/ħc^2))
    v_rho = -v₀_nn 
    R_rho = r₀*A^(1/3)
    a_rho = 0.67 

    # angular matrix elements 
    MEs_ang::Array{Float64, 4} = calc_MEs_ang(lmax)
    MEs_ang2::Array{Float64, 4} = calc_MEs_ang2(lmax)
end


function calc_MEs_ang(lmax)
    dim = 0
    for l in 0:lmax, j in 2l+1: -2: max(2l-1,0), m in j : -2: -j
        dim += 1
    end
    #@show dim

    MEs_ang = zeros(Float64, dim, dim, dim, dim)

    prog = Progress(dim^4, 1, "Calculating angular matrix elements...")

    n₄_ljm = 0
    for l₄ in 0:lmax, j₄ in 2l₄+1: -2: max(2l₄-1,0), m₄ in j₄: -2: -j₄
        n₄_ljm += 1

        n₃_ljm = 0
        for l₃ in 0:lmax, j₃ in 2l₃+1: -2: max(2l₃-1,0), m₃ in j₃: -2: -j₃
            n₃_ljm += 1

            n₂_ljm = 0
            for l₂ in 0:lmax, j₂ in 2l₂+1: -2: max(2l₂-1,0), m₂ in j₂ : -2: -j₂
                n₂_ljm += 1

                n₁_ljm = 0
                for l₁ in 0:lmax, j₁ in 2l₁+1: -2: max(2l₁-1,0), m₁ in j₁ : -2: -j₁
                    n₁_ljm += 1

                    if m₁+m₂ ≠ m₃+m₄ 
                        next!(prog)
                        continue
                    end
                    M = div(m₁+m₂, 2)

                    # angular matrix element (L=0)
                    ME_ang = 0.0
                    Jmin = max(div(abs(j₁-j₂),2), div(abs(j₃-j₄),2), abs(M))
                    Jmax = min(div(j₁+j₂,2), div(j₃+j₄,2))
                    for J in Jmin:Jmax 
                        if isodd(l₁+l₂+J) || isodd(l₃+l₄+J)
                            continue 
                        end
                        ME_ang += (-1)^(l₁+l₃) * sqrt((j₁+1)*(j₃+1))/4π *
                        clebsch(j₁, 1, 2J, 0, j₂, 1) * clebsch(j₃, 1, 2J, 0, j₄, 1) * 
                        clebsch(j₁, m₁, j₂, m₂, 2J, 2M) * clebsch(j₃, m₃, j₄, m₄, 2J, 2M)
                    end 

                    MEs_ang[n₁_ljm, n₂_ljm, n₃_ljm, n₄_ljm] = ME_ang

                    next!(prog)
                end
            end
        end
    end

    return MEs_ang
end

function calc_MEs_ang2(lmax)
    dim = 0
    for l in 0:lmax, j in 2l+1: -2: max(2l-1,0), m in j : -2: -j
        dim += 1
    end
    #@show dim

    MEs_ang2 = zeros(Float64, dim, dim, dim, dim)

    prog = Progress(dim^4, 1, "Calculating angular matrix elements...")

    n₄_ljm = 0
    for l₄ in 0:lmax, j₄ in 2l₄+1: -2: max(2l₄-1,0), m₄ in j₄: -2: -j₄
        n₄_ljm += 1

        n₃_ljm = 0
        for l₃ in 0:lmax, j₃ in 2l₃+1: -2: max(2l₃-1,0), m₃ in j₃: -2: -j₃
            n₃_ljm += 1

            n₂_ljm = 0
            for l₂ in 0:lmax, j₂ in 2l₂+1: -2: max(2l₂-1,0), m₂ in j₂ : -2: -j₂
                n₂_ljm += 1

                n₁_ljm = 0
                for l₁ in 0:lmax, j₁ in 2l₁+1: -2: max(2l₁-1,0), m₁ in j₁ : -2: -j₁
                    n₁_ljm += 1

                    if m₁+m₂ ≠ m₃+m₄ 
                        next!(prog)
                        continue
                    end
                    M = div(m₁+m₂, 2)

                    # angular matrix element (L=0)
                    ME_ang = 0.0
                    for J₃₄ in max(div(abs(j₃-j₄),2), abs(M)) : div(j₃+j₄,2)
                        if isodd(l₃+l₄+J₃₄) 
                            continue 
                        end
                        for J₁₂ in max(div(abs(j₁-j₂),2), abs(M)) : div(j₁+j₂,2)
                            if isodd(l₁+l₂+J₁₂)
                                continue 
                            end
                            temp = (-1)^M * sqrt((2J₁₂+1)*(2J₃₄+1)/20π) * 
                            clebsch(2J₁₂,   0, 2J₃₄,  0, 2*2, 0) * 
                            clebsch(2J₁₂, -2M, 2J₃₄, 2M, 2*2, 0)
                            if isapprox(temp, 0.0; rtol=1e-10)
                                continue 
                            end

                            ME_ang += (-1)^(l₁+l₃) * sqrt((j₁+1)*(j₃+1))/4π * temp * 
                            clebsch(j₁, 1, 2J₁₂, 0, j₂, 1) * clebsch(j₃, 1, 2J₃₄, 0, j₄, 1) * 
                            clebsch(j₁, m₁, j₂, m₂, 2J₁₂, 2M) * clebsch(j₃, m₃, j₄, m₄, 2J₃₄, 2M)
                        end
                    end

                    MEs_ang2[n₁_ljm, n₂_ljm, n₃_ljm, n₄_ljm] = ME_ang

                    next!(prog)
                end
            end
        end
    end

    return MEs_ang2
end






end # module
