
export test_calc_Vnn_matrix_element, test_make_three_body_Hamiltonian

#=
function Vnn_expectation_value(param; β=0.0, istate=1)
    @unpack Nr, rs, Δr = param

    # make single-particle bases
    spbases = make_spbases(param)
    @unpack nbases, ψs = spbases 

    # calculate single-particle states
    spstates = calc_single_particle_states(param, spbases, β)
    @unpack nstates, coeffs, qnums = spstates 
    @unpack Λ, Π = spstates.qnums[istate]
    
    @show nbases nstates
     
    # calculate radial matrix elements
    MEs_rad = zeros(Float64, nbases, nbases, nbases, nbases)
    @time for b₄ in 1:nbases, b₃ in 1:nbases, b₂ in 1:nbases, b₁ in 1:nbases
        ME_rad = 0.0
        for ir in 1:Nr
            ME_rad += ψs[ir,b₁]*ψs[ir,b₂]*ψs[ir,b₃]*ψs[ir,b₄]/rs[ir]^2
        end
        ME_rad *= Δr
        MEs_rad[b₁, b₂, b₃, b₄] = ME_rad
    end

    # calculate angular matrix elements 
    MEs_ang = zeros(Float64, nbases, nbases, nbases, nbases)
    @time for b₄ in 1:nbases
        l₄ = spbases.qnums[b₄].l
        j₄ = spbases.qnums[b₄].j
        if j₄ < abs(Λ) || (-1)^l₄ ≠ Π
            continue 
        end

        for b₃ in 1:nbases
            l₃ = spbases.qnums[b₃].l
            j₃ = spbases.qnums[b₃].j
            if j₃ < abs(Λ) || (-1)^l₃ ≠ Π
                continue 
            end

            for b₂ in 1:nbases
                l₂ = spbases.qnums[b₂].l
                j₂ = spbases.qnums[b₂].j
                if j₂ < abs(Λ) || (-1)^l₂ ≠ Π 
                    continue 
                end

                for b₁ in 1:nbases 
                    l₁ = spbases.qnums[b₁].l
                    j₁ = spbases.qnums[b₁].j
                    if j₁ < abs(Λ) || (-1)^l₁ ≠ Π 
                        continue 
                    end

                    ME_ang = 0.0
                    Lmin = max(div(abs(j₁-j₃),2), div(abs(j₂-j₄),2))
                    Lmax = min(div(j₁+j₃,2), div(j₂+j₄,2))
                    for L in Lmin:Lmax
                        ME₁ = calc_angular_matrix_element(l₁,j₁,+Λ,L,0,l₃,j₃,+Λ)
                        ME₂ = calc_angular_matrix_element(l₂,j₂,-Λ,L,0,l₄,j₄,-Λ)
                        ME_ang += ME₁*ME₂
                    end
                    MEs_ang[b₁, b₂, b₃, b₄] = ME_ang

                end
            end
        end
    end
    
    # calculate Vnn expectation value 
    V = 0.0
    @time for b₄ in 1:nbases
        l₄ = spbases.qnums[b₄].l
        j₄ = spbases.qnums[b₄].j
        if j₄ < abs(Λ) || (-1)^l₄ ≠ Π
            continue 
        end

        for b₃ in 1:nbases
            l₃ = spbases.qnums[b₃].l
            j₃ = spbases.qnums[b₃].j
            if j₃ < abs(Λ) || (-1)^l₃ ≠ Π
                continue 
            end

            for b₂ in 1:nbases
                l₂ = spbases.qnums[b₂].l
                j₂ = spbases.qnums[b₂].j
                if j₂ < abs(Λ) || (-1)^l₂ ≠ Π 
                    continue 
                end

                for b₁ in 1:nbases 
                    l₁ = spbases.qnums[b₁].l
                    j₁ = spbases.qnums[b₁].j
                    if j₁ < abs(Λ) || (-1)^l₁ ≠ Π 
                        continue 
                    end

                    temp = 1.0
        
                    temp *= coeffs[b₁,istate]*coeffs[b₂,istate]*
                    coeffs[b₃,istate]*coeffs[b₄,istate]

                    if isodd(div(j₂-Λ,2)) 
                        temp *= -1
                    end

                    if isodd(div(j₄-Λ,2))
                        temp *= -1
                    end
                    
                    if temp ≈ 0
                        continue
                    end
                    
                    V += temp*MEs_rad[b₁,b₂,b₃,b₄]*MEs_ang[b₁,b₂,b₃,b₄]
                end
            end
        end
    end
    
    return V
end
=#

#=
function Vnn_expectation_value(param; β=0.0, istate=1)
    @unpack Nr, rs, Δr, lmax = param
    num_lj = 2lmax+1

    # make single-particle bases
    spbases = make_spbases(param)

    # calculate single-particle states
    spstates = calc_single_particle_states(param, spbases, β)
    @unpack nstates, ψs, qnums = spstates 
    @unpack Λ, Π = spstates.qnums[istate]

    @show num_lj nstates

    n = 0
    for l in 0:lmax, j in 2l+1: -2: max(2l-1,0)
        n += j+1
    end
    @show n
    
    # calculate Vnn expectation value 
    V = 0.0
    @time for l₄ in 0:lmax, j₄ in 2l₄+1: -2: max(2l₄-1,0)
        if j₄ < abs(Λ) || (-1)^l₄ ≠ Π
            continue 
        end
        n₄_lj = calc_n_lj(l₄, j₄)

        for l₃ in 0:lmax, j₃ in 2l₃+1: -2: max(2l₃-1,0)
            if j₃ < abs(Λ) || (-1)^l₃ ≠ Π
                continue 
            end
            n₃_lj = calc_n_lj(l₃, j₃)

            for l₂ in 0:lmax, j₂ in 2l₂+1: -2: max(2l₂-1,0)
                if j₂ < abs(Λ) || (-1)^l₂ ≠ Π 
                    continue 
                end
                n₂_lj = calc_n_lj(l₂, j₂)

                for l₁ in 0:lmax, j₁ in 2l₁+1: -2: max(2l₁-1,0)
                    if j₁ < abs(Λ) || (-1)^l₁ ≠ Π 
                        continue 
                    end
                    n₁_lj = calc_n_lj(l₁, j₁)

                    # radial matrix element 
                    ME_rad = 0.0
                    for ir in 1:Nr
                        ME_rad += ψs[ir, n₁_lj, istate] *
                                    ψs[ir, n₂_lj, istate] *
                                    ψs[ir, n₃_lj, istate] * 
                                    ψs[ir, n₄_lj, istate] / rs[ir]^2
                    end
                    ME_rad *= Δr

                    # angular matrix element 
                    ME_ang = 0.0
                    Lmin = max(div(abs(j₁-j₃),2), div(abs(j₂-j₄),2))
                    Lmax = min(div(j₁+j₃,2), div(j₂+j₄,2))
                    for L in Lmin:Lmax
                        ME₁ = calc_angular_matrix_element(l₁,j₁,+Λ,L,0,l₃,j₃,+Λ)
                        ME₂ = calc_angular_matrix_element(l₂,j₂,-Λ,L,0,l₄,j₄,-Λ)
                        ME_ang += ME₁*ME₂
                    end

                    phase = 1.0
                    if isodd(div(j₂-Λ,2)) 
                        phase *= -1
                    end
                    if isodd(div(j₄-Λ,2))
                        phase *= -1
                    end

                    V += phase*ME_rad*ME_ang


                end
            end
        end
    end

    @time V2 = calc_Vnn_matrix_element(param, spstates, 1, 2, 1, 2)
    @show V V2

    return
end
=#


function calc_Vnn_matrix_element(param, spstates, n₁, n₂, n₃, n₄)
    @unpack Nr, rs, Δr, lmax = param

    @unpack nstates, ψs, qnums = spstates 

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

    V = 0.0
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
                        ME_rad += ψs[ir, n₁_lj, i₁] *
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

                    V += phase*ME_rad*ME_ang
                end
            end
        end
    end

    return V

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
    
    dim_max = div(nstates*(nstates+1), 2)
    @show Emax lmax dim_max

    Hmat_3body = zeros(Float64, dim_max, dim_max)

    n₃₄ = 0
    @time for n₄ in 1:2nstates
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

                    

                end
            end
        end
    end

    @show n₃₄

    return Hmat_3body[1:n₃₄, 1:n₃₄]
end

function test_make_three_body_Hamiltonian(param; β=0.0, Λ=0, Π=1)
    spbases = make_spbases(param)
    spstates = calc_single_particle_states(param, spbases, β)
    calc_occ!(spstates, param)

    make_three_body_Hamiltonian(param, spstates, Λ, Π)
end