
export test_calc_Vnn_matrix_element, test_make_three_body_Hamiltonian



function calc_Vnn_matrix_element(param, spstates, n₁, n₂, n₃, n₄)
    @unpack Nr, rs, Δr, lmax, v₀_nn, v_rho, R_rho, a_rho = param

    @unpack nstates, ψs, qnums = spstates 

    Vnn = zeros(Float64, Nr)
    @. Vnn = v₀_nn + v_rho/(1 + exp((rs - R_rho)/a_rho))

    i₁ = cld(n₁,2)
    i₂ = cld(n₂,2)
    i₃ = cld(n₃,2)
    i₄ = cld(n₄,2)

    Λ₁ = qnums[i₁].Λ
    Π₁ = qnums[i₁].Π
    if iseven(n₁)
        Λ₁ = -Λ₁
    end

    Λ₂ = qnums[i₂].Λ
    Π₂ = qnums[i₂].Π
    if iseven(n₂)
        Λ₂ = -Λ₂
    end

    Λ₃ = qnums[i₃].Λ
    Π₃ = qnums[i₃].Π
    if iseven(n₃)
        Λ₃ = -Λ₃
    end

    Λ₄ = qnums[i₄].Λ
    Π₄ = qnums[i₄].Π
    if iseven(n₄)
        Λ₄ = -Λ₄
    end

    if Λ₁+Λ₂ ≠ Λ₃+Λ₄ 
        return 0.0
    end
    M = div(Λ₁ - Λ₃, 2)

    ME_Vnn = 0.0
    for l₄ in 0:lmax, j₄ in 2l₄+1: -2: max(2l₄-1,0)
        if j₄ < abs(Λ₄) || (-1)^l₄ ≠ Π₄
            continue 
        end
        n₄_lj = calc_n_lj(l₄, j₄)

        for l₃ in 0:lmax, j₃ in 2l₃+1: -2: max(2l₃-1,0)
            if j₃ < abs(Λ₃) || (-1)^l₃ ≠ Π₃
                continue 
            end
            n₃_lj = calc_n_lj(l₃, j₃)

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
                        ME_rad += Vnn[ir] * 
                                    ψs[ir, n₁_lj, i₁] *
                                    ψs[ir, n₂_lj, i₂] *
                                    ψs[ir, n₃_lj, i₃] * 
                                    ψs[ir, n₄_lj, i₄] / rs[ir]^2
                    end
                    ME_rad *= Δr

                    # angular matrix element 
                    ME_ang = 0.0
                    Lmin = max(div(abs(j₁-j₃),2), div(abs(j₂-j₄),2))
                    Lmax = min(div(j₁+j₃,2), div(j₂+j₄,2))
                    for L in Lmin:Lmax
                        if L < abs(M) 
                            continue 
                        end
                        ME₁ = calc_angular_matrix_element(l₁,j₁,Λ₁,L,M,l₃,j₃,Λ₃)
                        ME₂ = calc_angular_matrix_element(l₂,j₂,Λ₂,L,M,l₄,j₄,Λ₄)
                        ME_ang += ME₁*ME₂
                    end

                    phase = 1.0
                    if isodd(div(j₂-Λ₂,2)) 
                        phase *= -1
                    end
                    if isodd(div(j₄-Λ₄,2))
                        phase *= -1
                    end

                    ME_Vnn += phase*ME_rad*ME_ang
                end
            end
        end
    end

    return ME_Vnn

end

function test_calc_Vnn_matrix_element(param; β=0.0, istate=1)
    @unpack Emax, lmax = param 

    spbases = make_spbases(param)
    spstates = calc_single_particle_states(param, spbases, β)
    @unpack nstates = spstates 

    @show Emax, lmax, nstates

    n₁ = 2istate-1
    n₂ = 2istate
    @time calc_Vnn_matrix_element(param, spstates, n₁, n₂, n₁, n₂)
end






function make_three_body_Hamiltonian(param, spstates, Λ, Π)
    @unpack Nr, rs, Δr, Emax, lmax = param
    num_lj = 2lmax+1

    @unpack nstates, ψs, spEs, qnums, occ = spstates 
    #show_spstates(spstates)

    # calculate the size of three-body Hamiltonian
    dim = 0
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

            if Λ ≠ Λ₁ + Λ₂ || Π ≠ Π₁*Π₂ || spE₁ + spE₂ > Emax
                continue 
            end
            dim += 1

        end
    end
    @show Emax lmax nstates dim
    
    Hmat_3body = zeros(Float64, dim, dim)

    prog = Progress(dim,1)

    n₃₄ = 0
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

            if Λ ≠ Λ₃ + Λ₄ || Π ≠ Π₃*Π₄ || spE₃ + spE₄ > Emax 
                continue 
            end
            n₃₄ += 1

            Hmat_3body[n₃₄, n₃₄] += spEs[i₃] + spEs[i₄]

            # show progress
            next!(prog)

            n₁₂ = 0
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

                    if Λ ≠ Λ₁ + Λ₂ || Π ≠ Π₁*Π₂ || spE₁ + spE₂ > Emax
                        continue 
                    end
                    n₁₂ += 1

                    Hmat_3body[n₁₂, n₃₄] += 
                        calc_Vnn_matrix_element(param, spstates, n₁, n₂, n₃, n₄)

                end
            end
        end
    end

    return Hmat_3body
end

function test_make_three_body_Hamiltonian(param; β=0.0, Λ=0, Π=1)
    spbases = make_spbases(param)
    spstates = calc_single_particle_states(param, spbases, β)
    calc_occ!(spstates, param)

    show_spstates(spstates)

    Hmat_3body = make_three_body_Hamiltonian(param, spstates, Λ, Π)

    @show eigmin(Hmat_3body)
    return 
end