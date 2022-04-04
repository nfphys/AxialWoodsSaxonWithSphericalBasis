
export test_calc_2pwfs_at_same_position, test_calc_greens_function

function calc_irJ(Nr, ir, J) 
    if iseven(J)
        return ir + Nr*div(J,2)
    else
        return ir + Nr*div(J-1,2)
    end
end

function calc_2pwf_at_same_position(param, spstates, n₁, n₂, ir, J)
    @unpack Nr, rs, Δr, Emax, lmax = param
    @unpack nstates, ψs, spEs, qnums, occ = spstates 

    @assert 1 ≤ ir ≤ Nr

    i₁ = cld(n₁, 2)
    Λ₁ = qnums[i₁].Λ
    Π₁ = qnums[i₁].Π
    if iseven(n₁)
        Λ₁ = -Λ₁
    end

    i₂ = cld(n₂, 2)
    Λ₂ = qnums[i₂].Λ
    Π₂ = qnums[i₂].Π
    if iseven(n₂)
        Λ₂ = -Λ₂
    end

    M = div(Λ₁+Λ₂, 2)

    ψ_2p = 0.0 

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

            ψ_2p += clebsch(j₁, Λ₁, j₂, Λ₂, 2J, 2M) *
            ψs[ir, n₁_lj, i₁]/rs[ir] * ψs[ir, n₂_lj, i₂]/rs[ir] * 
            (-1)^l₁ * sqrt((j₁+1)/4π) * clebsch(j₁, 1, 2J, 0, j₂, 1)
        end
    end

    return ψ_2p
end

function calc_2pwfs_at_same_position(param, spstates, Λ, Π)
    @unpack Nr, rs, Δr, Jmax, Emax, lmax = param
    @unpack nstates, ψs, spEs, qnums, occ = spstates 

    dim = 0
    ns_2p = zeros(Int64, 2, 2nstates*2nstates)
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
            ns_2p[1,dim] = n₁
            ns_2p[2,dim] = n₂
        end
    end
    @show dim

    prog = Progress(dim, 1, "Calculating two-particle wave functions...")

    dim_G = 0
    for J in 0:Jmax
        if (-1)^J !== Π 
            continue 
        end
        dim_G += Nr
    end

    ψs_2p = zeros(Float64, dim_G, dim)

    @threads for n₁₂ in 1:dim 
        n₁ = ns_2p[1,n₁₂]
        n₂ = ns_2p[2,n₁₂]

        for J in 0:Jmax, ir in 1:Nr 
            if (-1)^J !== Π 
                continue 
            end
            irJ = calc_irJ(Nr, ir, J)
            ψs_2p[irJ, n₁₂] = calc_2pwf_at_same_position(param, spstates, n₁, n₂, ir, J)
        end

        next!(prog)
    end

    return ψs_2p
end


function test_calc_2pwfs_at_same_position(param; β=0.0, Λ=0, Π=1)
    @unpack Nr, rs, Δr, Emax, lmax = param

    spbases = make_spbases(param)
    spstates = calc_single_particle_states(param, spbases, β)
    calc_occ!(spstates, param)
    @unpack nstates, ψs, spEs, qnums, occ = spstates 

    calc_2pwfs_at_same_position(param, spstates, Λ, Π)
    return 
end




function calc_greens_function(param, spstates, ψs_2p, Λ, Π, E, η)
    @unpack Nr, rs, Δr, Jmax, Emax, lmax = param
    @unpack nstates, ψs, spEs, qnums, occ = spstates 

    dim_G = 0
    for J in 0:Jmax
        if (-1)^J !== Π 
            continue 
        end
        dim_G += Nr
    end
    Gs = zeros(ComplexF64, dim_G, dim_G)

    prog = Progress(dim_G, 1, "Calculating Green's function...")

    @threads for irJ₂ in 1:dim_G 

        for irJ₁ in 1:dim_G 

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

                    Gs[irJ₁, irJ₂] += ψs_2p[irJ₁, n₁₂] * ψs_2p[irJ₂, n₁₂] / (spE₁ + spE₂ - E - im*η)
                end
            end

        end
    end
    
    return Gs
end



function test_calc_greens_function(param; β=0.0, Λ=0, Π=1, E=-3.0, η=0.2)
    @unpack Nr, rs, Δr, Jmax, Emax, lmax = param

    spbases = make_spbases(param)
    spstates = calc_single_particle_states(param, spbases, β)
    calc_occ!(spstates, param)
    @unpack nstates, ψs, spEs, qnums, occ = spstates 

    ψs_2p = calc_2pwfs_at_same_position(param, spstates, Λ, Π)

    Gs = calc_greens_function(param, spstates, ψs_2p, Λ, Π, E, η)
    @show det(Gs)
    return 
end