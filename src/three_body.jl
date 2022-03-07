
export test_calc_BE1_strength

function calc_dipole_matrix_element(param, spstates, n₁, n₂, M)
    @unpack Z, A, Nr, rs, Δr, Emax, lmax = param
    @unpack nstates, ψs, spEs, qnums, occ = spstates 

    i₁ = cld(n₁,2)
    Λ₁ = qnums[i₁].Λ
    Π₁ = qnums[i₁].Π
    if iseven(n₁)
        Λ₁ = -Λ₁
    end

    i₂ = cld(n₂,2)
    Λ₂ = qnums[i₂].Λ
    Π₂ = qnums[i₂].Π
    if iseven(n₂)
        Λ₂ = -Λ₂
    end

    ME_dipole = 0.0
    for l₂ in 0:lmax, j₂ in 2l₂+1: -2: max(2l₂-1,0)
        if j₂ < abs(Λ₂) || (-1)^l₂ ≠ Π₂
            continue 
        end
        n₂_lj = calc_n_lj(l₂, j₂)

        for l₁ in 0:lmax, j₁ in 2l₁+1: -2: max(2l₁-1,0)
            if j₁ < abs(Λ₁) || (-1)^l₁ ≠ Π₁
                continue 
            end
            n₁_lj = calc_n_lj(l₁, j₁)

            # radial matrix element 
            ME_rad = 0.0
            for ir in 1:Nr 
                ME_rad += ψs[ir, n₁_lj, i₁] * rs[ir] * ψs[ir, n₂_lj, i₂]
            end
            ME_rad *= Δr

            ME_ang = calc_angular_matrix_element(l₁, j₁, Λ₁, 1, M, l₂, j₂, Λ₂)
            
            ME_dipole += -(Z/A) * ME_rad * ME_ang 
        end
    end
    
    return ME_dipole
end



function calc_BE1_strength(param, spstates, coeff_gs, coeff_excited, M)
    @assert M === 0 || M === -1 || M === +1
    
    @unpack Z, A, Nr, rs, Δr, Emax, lmax = param
    @unpack nstates, ψs, spEs, qnums, occ = spstates 

    Λ_gs = 0
    Π_gs = 1
    dim_gs = length(coeff_gs)

    Λ_excited = 2M 
    Π_excited = -1
    dim_excited = length(coeff_excited)

    #prog = Progress(dim_gs*dim_excited, 1, "Calculating B(E1) strength...")

    ME = 0.0 # matrix element

    n₃₄ = 0 # ground states
    for n₄ in 1:2nstates
        i₄ = cld(n₄, 2)
        if occ[i₄] == 1.0
            continue 
        end

        spE₄ = spEs[i₄]
        Λ₄ = qnums[i₄].Λ
        Π₄ = qnums[i₄].Π
        if iseven(n₄)
            Λ₄ = -Λ₄
        end

        for n₃ in 1:n₄-1 # n₃ < n₄
            i₃ = cld(n₃, 2)
            if occ[i₃] == 1.0
                continue 
            end

            spE₃ = spEs[i₃]
            Λ₃ = qnums[i₃].Λ
            Π₃ = qnums[i₃].Π
            if iseven(n₃)
                Λ₃ = -Λ₃
            end

            if Λ_gs ≠ Λ₃ + Λ₄ || Π_gs ≠ Π₃*Π₄ || spE₃ + spE₄ > Emax 
                continue 
            end
            n₃₄ += 1

            n₁₂ = 0 # excited states 
            for n₂ in 1:2nstates 
                i₂ = cld(n₂, 2)
                if occ[i₂] == 1.0
                    continue 
                end 

                spE₂ = spEs[i₂]
                Λ₂ = qnums[i₂].Λ
                Π₂ = qnums[i₂].Π
                if iseven(n₂)
                    Λ₂ = -Λ₂
                end

                for n₁ in 1:n₂-1 # n₁ < n₂
                    i₁ = cld(n₁, 2)
                    if occ[i₁] == 1.0
                        continue
                    end

                    spE₁ = spEs[i₁]
                    Λ₁ = qnums[i₁].Λ
                    Π₁ = qnums[i₁].Π
                    if iseven(n₁)
                        Λ₁ = -Λ₁
                    end

                    if Λ_excited ≠ Λ₁ + Λ₂ || Π_excited ≠ Π₁*Π₂ || spE₁ + spE₂ > Emax
                        continue 
                    end
                    n₁₂ += 1

                    if n₂ === n₄ 
                        ME += coeff_excited[n₁₂] * coeff_gs[n₃₄] * 
                        calc_dipole_matrix_element(param, spstates, n₁, n₃, M)
                    end

                    if n₁ === n₃ 
                        ME += coeff_excited[n₁₂] * coeff_gs[n₃₄] *
                        calc_dipole_matrix_element(param, spstates, n₂, n₄, M)
                    end

                    if n₂ === n₃
                        ME -= coeff_excited[n₁₂] * coeff_gs[n₃₄] * 
                        calc_dipole_matrix_element(param, spstates, n₁, n₄, M) 
                    end

                    if n₁ === n₄ 
                        ME -= coeff_excited[n₁₂] * coeff_gs[n₃₄] *
                        calc_dipole_matrix_element(param, spstates, n₂, n₃, M)
                    end

                    # show progress
                    #next!(prog)

                end
            end

            @assert n₁₂ === dim_excited
        end
    end
    @assert n₃₄ === dim_gs

    BE1 = ME^2

    return BE1
end

function test_calc_BE1_strength(param; β=0.0, howmany=20, Γ=0.2)
    spbases = make_spbases(param)
    spstates = calc_single_particle_states(param, spbases, β)
    calc_occ!(spstates, param)

    show_spstates(spstates; Emax=1.0)

    # ground state
    Λ_gs = 0
    Π_gs = 1
    Hmat_3body = make_three_body_Hamiltonian(param, spstates, Λ_gs, Π_gs)
    @time Es_gs, coeffs_gs, info_gs = eigsolve(Hmat_3body, 1, :SR, eltype(Hmat_3body))
    E_gs = Es_gs[1]
    coeff_gs = coeffs_gs[1]
    @show E_gs

    Es = range(0, 10.0, step=0.01)
    fs = zeros(Float64, length(Es))
    p = plot(xlabel="E [MeV]", ylabel="B(E1)", title="β=$β  Γ=$Γ[MeV]")
    for M in -1:1
        Λ_excited = 2M
        Π_excited = -1 
        Hmat_3body = make_three_body_Hamiltonian(param, spstates, Λ_excited, Π_excited)
        @time Es_excited, coeffs_excited, info_excited = eigsolve(Hmat_3body, howmany, :SR, eltype(Hmat_3body))
        @show Es_excited[1:2]

        BE1s = zeros(Float64, length(Es_excited))
        @time for k in 1:length(Es_excited)
            BE1s[k] = calc_BE1_strength(param, spstates, coeff_gs, coeffs_excited[k], M)
        end
        @show BE1s[1:2]

        gs = zeros(Float64, length(Es))
        for k in 1:length(Es_excited)
            for iE in 1:length(Es)
                E = Es[iE]
                gs[iE] += (Γ/π) * 1/((E - Es_excited[k] + E_gs)^2 + Γ^2) * BE1s[k]
            end
        end
        plot!(p, Es, gs; label="M=$M")
        @. fs += gs
    end
    plot!(p, Es, fs; label="total")
    display(p)
end