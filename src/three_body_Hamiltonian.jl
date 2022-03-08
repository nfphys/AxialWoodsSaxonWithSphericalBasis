
export test_calc_Vnn_matrix_element, test_make_three_body_Hamiltonian, 
plot_ground_state_energy

#=
function calc_Vnn_matrix_element(param, spstates, n₁, n₂, n₃, n₄)
    @unpack Nr, rs, Δr, lmax, v₀_nn, v_rho, R_rho, a_rho = param

    @unpack nstates, ψs, qnums = spstates 

    Vnn = zeros(Float64, Nr)
    @. Vnn = v₀_nn + v_rho/(1 + exp((rs - R_rho)/a_rho))

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

    i₃ = cld(n₃,2)
    Λ₃ = qnums[i₃].Λ
    Π₃ = qnums[i₃].Π
    if iseven(n₃)
        Λ₃ = -Λ₃
    end

    i₄ = cld(n₄,2)
    Λ₄ = qnums[i₄].Λ
    Π₄ = qnums[i₄].Π
    if iseven(n₄)
        Λ₄ = -Λ₄
    end

    if Λ₁+Λ₂ ≠ Λ₃+Λ₄ 
        return 0.0
    end
    #M = div(Λ₃ - Λ₁, 2)
    M = div(Λ₁+Λ₂, 2)

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

                    phase = 1.0
                    if iseven(n₁) && isodd(div(j₁-Λ₁,2)) 
                        phase *= -1
                    end
                    if iseven(n₂) && isodd(div(j₂-Λ₂,2))
                        phase *= -1
                    end
                    if iseven(n₃) && isodd(div(j₃-Λ₃,2))
                        phase *= -1 
                    end
                    if iseven(n₄) && isodd(div(j₄-Λ₄,2))
                        phase *= -1 
                    end

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
                    Jmin = max(div(abs(j₁-j₂),2), div(abs(j₃-j₄),2), abs(M))
                    Jmax = min(div(j₁+j₂,2), div(j₃+j₄,2))
                    for J in Jmin:Jmax 
                        if isodd(l₁+l₂+J) || isodd(l₃+l₄+J)
                            continue 
                        end
                        ME_ang += (-1)^(l₁+l₃) * sqrt((j₁+1)*(j₃+1))/4π *
                        clebsch(j₁, 1, 2J, 0, j₂, 1) * clebsch(j₃, 1, 2J, 0, j₄, 1) * 
                        clebsch(j₁, Λ₁, j₂, Λ₂, 2J, 2M) * clebsch(j₃, Λ₃, j₄, Λ₄, 2J, 2M)
                    end

                    ME_Vnn += phase*ME_rad*ME_ang
                end
            end
        end
    end

    return ME_Vnn

end
=#


function calc_Vnn_matrix_element(param, spstates, β, n₁, n₂, n₃, n₄)
    @unpack Nr, rs, Δr, lmax, v₀_nn, v_rho, R_rho, a_rho = param

    @unpack nstates, ψs, qnums = spstates 

    Vnn = zeros(Float64, Nr)
    @. Vnn = v₀_nn + v_rho/(1 + exp((rs - R_rho)/a_rho))

    Vnn2 = zeros(Float64, Nr) 
    @. Vnn2 = β*(R_rho/a_rho)*v_rho*exp((rs - R_rho)/a_rho)/(1 + exp((rs - R_rho)/a_rho))^2

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

    i₃ = cld(n₃,2)
    Λ₃ = qnums[i₃].Λ
    Π₃ = qnums[i₃].Π
    if iseven(n₃)
        Λ₃ = -Λ₃
    end

    i₄ = cld(n₄,2)
    Λ₄ = qnums[i₄].Λ
    Π₄ = qnums[i₄].Π
    if iseven(n₄)
        Λ₄ = -Λ₄
    end

    if Λ₁+Λ₂ ≠ Λ₃+Λ₄ 
        return 0.0
    end
    #M = div(Λ₃ - Λ₁, 2)
    M = div(Λ₁+Λ₂, 2)

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

                    phase = 1.0
                    if iseven(n₁) && isodd(div(j₁-Λ₁,2)) 
                        phase *= -1
                    end
                    if iseven(n₂) && isodd(div(j₂-Λ₂,2))
                        phase *= -1
                    end
                    if iseven(n₃) && isodd(div(j₃-Λ₃,2))
                        phase *= -1 
                    end
                    if iseven(n₄) && isodd(div(j₄-Λ₄,2))
                        phase *= -1 
                    end

                    # radial matrix element (L=0)
                    ME_rad = 0.0
                    for ir in 1:Nr
                        ME_rad += Vnn[ir] * 
                                    ψs[ir, n₁_lj, i₁] *
                                    ψs[ir, n₂_lj, i₂] *
                                    ψs[ir, n₃_lj, i₃] * 
                                    ψs[ir, n₄_lj, i₄] / rs[ir]^2
                    end
                    ME_rad *= Δr

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
                        clebsch(j₁, Λ₁, j₂, Λ₂, 2J, 2M) * clebsch(j₃, Λ₃, j₄, Λ₄, 2J, 2M)
                    end

                    ME_Vnn += phase*ME_rad*ME_ang

                    # radial matrix element (L=2)
                    ME_rad = 0.0
                    for ir in 1:Nr
                        ME_rad += Vnn2[ir] * 
                                    ψs[ir, n₁_lj, i₁] *
                                    ψs[ir, n₂_lj, i₂] *
                                    ψs[ir, n₃_lj, i₃] * 
                                    ψs[ir, n₄_lj, i₄] / rs[ir]^2
                    end
                    ME_rad *= Δr

                    # angular matrix element (L=2)
                    ME_ang = 0.0 
                    for J₃₄ in max(div(abs(j₃-j₄),2), abs(M)) : div(j₃+j₄,2)
                        if isodd(l₃+l₄+J₃₄) 
                            continue 
                        end
                        for J₁₂ in max(div(abs(j₁-j₂),2), abs(M)) : div(j₁+j₂,2)
                            if isodd(l₁+l₂+J₁₂)
                                continue 
                            end
                            temp = sqrt((2J₁₂+1)*(2J₃₄+1)/20π) * 
                            clebsch(2J₁₂,   0, 2J₃₄,  0, 2*2, 0) * 
                            clebsch(2J₁₂, -2M, 2J₃₄, 2M, 2*2, 0)
                            if isapprox(temp, 0.0; rtol=1e-10)
                                continue 
                            end

                            ME_ang += (-1)^(l₁+l₃) * sqrt((j₁+1)*(j₃+1))/4π * temp * 
                            clebsch(j₁, 1, 2J₁₂, 0, j₂, 1) * clebsch(j₃, 1, 2J₃₄, 0, j₄, 1) * 
                            clebsch(j₁, Λ₁, j₂, Λ₂, 2J₁₂, 2M) * clebsch(j₃, Λ₃, j₄, Λ₄, 2J₃₄, 2M)
                        end
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





#=
function make_three_body_Hamiltonian(param, spstates, Λ, Π)
    @assert iseven(Λ)
    @assert Π === +1 || Π === -1 
    
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
    println("")
    @show Emax, lmax, nstates
    @show Λ, Π, dim
    
    Hmat_3body = zeros(Float64, dim, dim)

    prog = Progress(div(dim*(dim+1), 2), 1, "Making three-body Hamiltonian...")

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

                    if n₁₂ > n₃₄ # n₁₂ ≤ n₃₄
                        continue 
                    end

                    Hmat_3body[n₁₂, n₃₄] += 
                    calc_Vnn_matrix_element(param, spstates, n₁, n₂, n₃, n₄)

                    # show progress
                    next!(prog)

                end
            end
        end
    end

    return Symmetric(Hmat_3body)
end
=#

function make_three_body_Hamiltonian(param, spstates, β, Λ, Π)
    @assert iseven(Λ)
    @assert Π === +1 || Π === -1 
    
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
    @show Emax, lmax, nstates
    @show Λ, Π, dim
    
    Hmat_3body = zeros(Float64, dim, dim)

    prog = Progress(div(dim*(dim+1), 2), 1, "Making three-body Hamiltonian...")

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

                    if n₁₂ > n₃₄ # n₁₂ ≤ n₃₄
                        continue 
                    end

                    Hmat_3body[n₁₂, n₃₄] += 
                    calc_Vnn_matrix_element(param, spstates, β, n₁, n₂, n₃, n₄)

                    # show progress
                    next!(prog)

                end
            end
        end
    end
    println("")

    return Symmetric(Hmat_3body)
end

function test_make_three_body_Hamiltonian(param; β=0.0, Λ=0, Π=1, howmany=1)
    spbases = make_spbases(param)
    spstates = calc_single_particle_states(param, spbases, β)
    calc_occ!(spstates, param)

    #show_spstates(spstates)

    Hmat_3body = make_three_body_Hamiltonian(param, spstates, β, Λ, Π)

    @time Emin = eigmin(Hmat_3body)
    @time Es, coeffs, info = eigsolve(Hmat_3body, howmany, :SR, eltype(Hmat_3body))
    @show Emin Es

    return 
end

function plot_ground_state_energy(param; β_max=0.2, β_min=-0.2, Δβ=0.02)
    @unpack Z, N, Emax, lmax, Λmax = param

    spbases = make_spbases(param)
    @unpack nbases = spbases

    βs = range(β_min, β_max; step=Δβ)
    Nβ = length(βs)

    Λ = 0
    Π = 1

    Es_gs = zeros(Float64, Nβ)
    for iβ in 1:Nβ 
        β = βs[iβ]
        @show β

        spstates = calc_single_particle_states(param, spbases, β)
        calc_occ!(spstates, param) 

        Hmat_3body = make_three_body_Hamiltonian(param, spstates, β, Λ, Π)
        Es, coeffs, info = eigsolve(Hmat_3body, 1, :SR, eltype(Hmat_3body))
        Es_gs[iβ] = Es[1]
    end
    p = plot(ylim=(-5, 0), legend=false, 
    xlabel="β", ylabel="ground-state energy [MeV]", 
    title="Z=$Z  N=$N  Emax=$(Emax)  lmax=$(lmax)")
    plot!(p, βs, Es_gs)
    display(p)
end