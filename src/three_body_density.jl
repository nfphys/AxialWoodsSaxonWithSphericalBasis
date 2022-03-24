
export test_calc_two_particle_density, test_calc_single_particle_density


function calc_uncorrelated_2pwf(param, spstates, n₁, n₂, ir, φ₁₂, σ₁, σ₂)
    @unpack Nr, rs, Δr, Emax, lmax = param
    @unpack nstates, ψs, spEs, qnums, occ = spstates 

    @assert 1 ≤ ir ≤ Nr
    r = rs[ir]

    P(l,m) = sqrt((2l+1)/4π * factorial(l-m)/factorial(l+m)) * legendre(l, m, 0)

    i₁ = cld(n₁, 2)
    spE₁ = spEs[i₁]
    Λ₁ = qnums[i₁].Λ
    Π₁ = qnums[i₁].Π
    if iseven(n₁)
        Λ₁ = -Λ₁
    end

    i₂ = cld(n₂, 2)
    spE₂ = spEs[i₂]
    Λ₂ = qnums[i₂].Λ
    Π₂ = qnums[i₂].Π
    if iseven(n₂)
        Λ₂ = -Λ₂
    end

    ψ₂ = 0.0 + 0.0im 

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

            if abs(Λ₁-σ₁) > 2l₁ || abs(Λ₂-σ₂) > 2l₂ 
                continue 
            end

            phase = 1.0
            if iseven(n₁) && isodd(div(j₁-Λ₁,2)-l₁) 
                phase *= -1
            end
            if iseven(n₂) && isodd(div(j₂-Λ₂,2)-l₂)
                phase *= -1
            end

            ψ₂ += phase * ψs[ir, n₁_lj, i₁]/rs[ir] * ψs[ir, n₂_lj, i₂]/rs[ir] * 
            clebsch_ls(l₁, j₁, Λ₁, σ₁) * clebsch_ls(l₂, j₂, Λ₂, σ₂) *
            P(l₁, div(Λ₁-σ₁,2)) * P(l₂, div(Λ₂-σ₂,2)) * exp(im*div(Λ₁-σ₁,2)*φ₁₂)
        end
    end

    return ψ₂
end


function calc_two_particle_density(param, spstates, Λ, Π, coeff, ir, φ₁₂, σ₁, σ₂)
    @unpack Nr, rs, Δr, Emax, lmax = param
    @unpack nstates, ψs, spEs, qnums, occ = spstates 

    @assert 1 ≤ ir ≤ Nr
    r = rs[ir]

    dim = length(coeff)
    #prog = Progress(dim, 1, "Calculating two-body density...")

    ψ₂ = 0.0 + 0.0im

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

            ψ₂ += coeff[n₁₂] * calc_uncorrelated_2pwf(param, spstates, n₁, n₂, ir,  φ₁₂, σ₁, σ₂)/√2
            ψ₂ -= coeff[n₁₂] * calc_uncorrelated_2pwf(param, spstates, n₂, n₁, ir,  φ₁₂, σ₁, σ₂)/√2

            # show progress
            #next!(prog)

        end
    end

    ρ₂ = abs2(ψ₂)
end


function test_calc_two_particle_density(param; Λ=0, Π=1, β=0.0, σ₁=1, σ₂=-1)
    @unpack Nr, rs, Δr, R₀, Emax, lmax = param

    spbases = make_spbases(param)
    spstates = calc_single_particle_states(param, spbases, β)
    calc_occ!(spstates, param)

    Hmat_3body = make_three_body_Hamiltonian(param, spstates, β, Λ, Π)
    @time Es, coeffs, info = eigsolve(Hmat_3body, 1, :SR, eltype(Hmat_3body))
    @show Es[1]
    coeff = coeffs[1]
    
    φs = range(0, π, length=100+1)
    Nφ = length(φs)

    ρ₁ = zeros(Float64, Nr)
    @time for ir in 1:Nr
        ρ₁[ir] = calc_single_particle_density(param, spstates, Λ, Π, coeff, ir, π/2)
    end

    ρ₂ = zeros(Float64, Nr, Nφ)
    prog = Progress(Nr*Nφ, 1, "Calculating two-body density...")
    @time for iφ in 1:Nφ, ir in 1:Nr 
        r = rs[ir]
        φ = φs[iφ]
        ρ₂[ir, iφ] = 2π*r^2 * 4π*r^2 * sin(φ) *
        calc_two_particle_density(param, spstates, Λ, Π, coeff, ir, φ, σ₁, σ₂)
        next!(prog)
    end

    ir = floor(Int, 5/Δr)
    p = plot(title="Emax=$(Emax)MeV  lmax=$(lmax)  β=$β", xlabel="φ/π", ylim=(0, 0.02))
    plot!(p, φs/π, ρ₂[ir,:]; label="ρ₂, r=$(rs[ir])fm")
    display(p)

    iφ = 1
    p = plot(title="Emax=$(Emax)MeV  lmax=$(lmax)  β=$β", xlabel="r [fm]")
    plot!(p, rs, ρ₂[:,iφ]; label="ρ₂, φ=$(φs[iφ])")
    #plot!(p, rs, ρ₁; label="ρ₁")
    display(p)

    return
    
    p = plot(xlabel="r [fm]", ylabel="φ/π", xlim=(0,20), 
    title="Emax=$Emax, lmax=$lmax, β=$β")
    heatmap!(p, rs, φs/π, ρ₂')
    display(p)
end









function calc_single_particle_density(param, spstates, Λ, Π, coeff, ir, θ)
    @unpack Nr, rs, Δr, Emax, lmax = param
    @unpack nstates, ψs, spEs, qnums, occ = spstates 

    @assert 1 ≤ ir ≤ Nr
    r = rs[ir]

    dim = length(coeff)

    P(l,m,θ) = sqrt((2l+1)/4π * factorial(l-m)/factorial(l+m)) * legendre(l, m, cos(θ))

    ρ₁ = 0.0

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

                    if n₂ === n₄ 
                        for l₃ in 0:lmax, j₃ in 2l₃+1: -2: max(2l₃-1,0)
                            if j₃ < abs(Λ₃) || (-1)^l₃ ≠ Π₃
                                continue 
                            end
                            n₃_lj = calc_n_lj(l₃, j₃)
                    
                            for l₁ in 0:lmax, j₁ in 2l₁+1: -2: max(2l₁-1,0)
                                if j₁ < abs(Λ₁) || (-1)^l₁ ≠ Π₁
                                    continue 
                                end
                                n₁_lj = calc_n_lj(l₁, j₁)

                                phase = 1.0
                                if iseven(n₁) && isodd(div(j₁-Λ₁,2)-l₁) 
                                    phase *= -1
                                end
                                if iseven(n₃) && isodd(div(j₃-Λ₃,2)-l₃)
                                    phase *= -1 
                                end

                                for σ in 1: -2: -1
                                    if abs(Λ₁-σ) > 2l₁ || abs(Λ₃-σ) > 2l₃
                                        continue 
                                    end
                                    ρ₁ += phase * coeff[n₁₂] * coeff[n₃₄] * 
                                    ψs[ir, n₁_lj, i₁]/rs[ir] * ψs[ir, n₃_lj, i₃]/rs[ir] * 
                                    clebsch_ls(l₁, j₁, Λ₁, σ) * 
                                    clebsch_ls(l₃, j₃, Λ₃, σ) *
                                    P(l₁, div(Λ₁-σ, 2), θ) * P(l₃, div(Λ₃-σ, 2), θ)
                                end
                                    
                            end
                        end
                    end

                    if n₂ === n₃ 
                        for l₄ in 0:lmax, j₄ in 2l₄+1: -2: max(2l₄-1,0)
                            if j₄ < abs(Λ₄) || (-1)^l₄ ≠ Π₄
                                continue 
                            end
                            n₄_lj = calc_n_lj(l₄, j₄)
                    
                            for l₁ in 0:lmax, j₁ in 2l₁+1: -2: max(2l₁-1,0)
                                if j₁ < abs(Λ₁) || (-1)^l₁ ≠ Π₁
                                    continue 
                                end
                                n₁_lj = calc_n_lj(l₁, j₁)

                                phase = 1.0
                                if iseven(n₁) && isodd(div(j₁-Λ₁,2)-l₁) 
                                    phase *= -1
                                end
                                if iseven(n₄) && isodd(div(j₄-Λ₄,2)-l₄)
                                    phase *= -1 
                                end

                                for σ in 1: -2: -1
                                    if abs(Λ₁-σ) > 2l₁ || abs(Λ₄-σ) > 2l₄
                                        continue 
                                    end
                                    ρ₁ -= phase * coeff[n₁₂] * coeff[n₃₄] * 
                                    ψs[ir, n₁_lj, i₁]/rs[ir] * ψs[ir, n₄_lj, i₄]/rs[ir] * 
                                    clebsch_ls(l₁, j₁, Λ₁, σ) * 
                                    clebsch_ls(l₄, j₄, Λ₄, σ) *
                                    P(l₁, div(Λ₁-σ, 2), θ) * P(l₄, div(Λ₄-σ, 2), θ)
                                end
                                    
                            end
                        end
                    end

                end
            end
        end
    end

    return ρ₁
end


function test_calc_single_particle_density(param; Λ=0, Π=1, β=0.0, θ=0.0)
    @unpack Nr, rs, Δr, R₀, Emax, lmax = param

    spbases = make_spbases(param)
    spstates = calc_single_particle_states(param, spbases, β)
    calc_occ!(spstates, param)
    show_spstates(spstates)

    Hmat_3body = make_three_body_Hamiltonian(param, spstates, β, Λ, Π)

    @time Es, coeffs, info = eigsolve(Hmat_3body, 1, :SR, eltype(Hmat_3body))
    @show Es[1]
    coeff = coeffs[1]

    #=
    ρ₁_rad = zeros(Float64, Nr)
    @time for ir in 1:Nr
        ρ₁_rad[ir] = calc_single_particle_density(param, spstates, Λ, Π, coeff, ir, θ)
    end

    sum_ρ₁_rad = sum(@. 4π*rs^2*ρ₁_rad)*Δr
    @show sum_ρ₁_rad

    p = plot()
    plot!(p, rs, ρ₁_rad)
    display(p)
    =#

    θs = range(0, π, length=100+1)
    Nθ = length(θs)
    Δθ = θs[2] - θs[1]
    ρ₁ = zeros(Float64, Nr, Nθ)
    prog = Progress(Nr*Nθ, 1, "Calculating single-particle density...")
    for iθ in 1:Nθ, ir in 1:Nr 
        ρ₁[ir, iθ] = calc_single_particle_density(param, spstates, Λ, Π, coeff, ir, θs[iθ])
        next!(prog)
    end

    sum_ρ₁ = 0.0
    for iθ in 1:Nθ, ir in 1:Nr 
        r = rs[ir]
        θ = θs[iθ]
        sum_ρ₁ += 2π*r^2*Δr * sin(θ)*Δθ  * ρ₁[ir, iθ]
    end
    @show sum_ρ₁

    p = plot(xlabel="r [fm]", ylabel="θ/π", xlim=(0,10))
    heatmap!(p, rs, θs/π, ρ₁')
    display(p)
end


