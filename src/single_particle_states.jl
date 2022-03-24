
export test_make_single_particle_Hamiltonian, test_calc_single_particle_states, 
plot_nilsson_diagram, plot_total_energy


@with_kw struct QuantumNumbers @deftype Int64
    Λ = 1; @assert isodd(Λ)
    Π = 1; @assert Π === 1 || Π === -1
end

@with_kw struct SingleParticleStates
    nstates::Int64
    coeffs::Matrix{Float64}; @assert size(coeffs, 2) === nstates
    ψs::Array{Float64, 3}; @assert size(ψs, 3) === nstates 
    spEs::Vector{Float64}; @assert length(spEs) === nstates 
    qnums::Vector{QuantumNumbers}; @assert length(qnums) === nstates 
    occ::Vector{Float64}; @assert length(occ) === nstates 
end



function make_single_particle_Hamiltonian(param, spbases, β, Λ, Π)
    @unpack V₀, R₀, a, Nr, Δr, rs = param 
    @unpack nbases, ψs, spEs, qnums = spbases 
    
    @assert isodd(Λ)
    @assert Π === 1 || Π === -1
    
    Hmat = zeros(Float64, nbases, nbases)
    
    Vs = zeros(Float64, Nr)
    @. Vs = (R₀/a)*V₀*exp((rs-R₀)/a)/(1+exp((rs-R₀)/a))^2
    
    #=
    Y(L,θ) = sqrt((2L+1)/4π) * legendre(L, 0, cos(θ)) # Y_LO(θ)

    Vs = zeros(Float64, Nr, Lmax_WS+1)
    θs = range(0, π, length=100+1); Nθ = length(θs); Δθ = θs[2]-θs[1]
    #p = plot()
    for L in 0:Lmax_WS
        for ir in 1:Nr 
            r = rs[ir]
            for iθ in 1:Nθ
                θ = θs[iθ]
                R = R₀*(1 + β*Y(2, θ))
                Vs[ir, 1+L] += sin(θ) * Y(L,θ) * V₀ / (1 + exp((r-R)/a))
            end
            Vs[ir, 1+L] *= 2π * Δθ 
        end
        #plot!(p, rs, Vs[:,1+L]; label="L=$L")
    end
    #display(p)
    =#
    
    n₂ = 0
    for i₂ in 1:nbases
        @views ψ₂ = ψs[:,i₂]
        l₂ = qnums[i₂].l
        j₂ = qnums[i₂].j
        if j₂ < abs(Λ) || (-1)^l₂ ≠ Π
            continue 
        end
        n₂ += 1
        
        Hmat[n₂, n₂] += spEs[i₂]
        
        n₁ = 0
        for i₁ in 1:nbases
            @views ψ₁ = ψs[:,i₁]
            l₁ = qnums[i₁].l
            j₁ = qnums[i₁].j
            if j₁ < abs(Λ) || (-1)^l₁ ≠ Π
                continue 
            end
            n₁ += 1
            
            # radial matrix element
            M_rad = 0.0
            for ir in 1:Nr
                M_rad += ψ₁[ir]*Vs[ir]*ψ₂[ir]
            end
            M_rad *= Δr
            
            # angular matrix element 
            M_ang = calc_angular_matrix_element(l₁,j₁,Λ,2,0,l₂,j₂,Λ)
            
            Hmat[n₁, n₂] += β*M_rad*M_ang

            #=
            for L in 1:Lmax_WS
                # radial matrix element
                M_rad = 0.0
                for ir in 1:Nr
                    M_rad += ψ₁[ir]*Vs[ir,1+L]*ψ₂[ir]
                end
                M_rad *= Δr
                
                # angular matrix element 
                M_ang = calc_angular_matrix_element(l₁,j₁,Λ,L,0,l₂,j₂,Λ)
                
                Hmat[n₁, n₂] += M_rad * M_ang
            end
            =#
        end
    end
    
    return Symmetric(Hmat[1:n₂,1:n₂])
end

function test_make_single_particle_Hamiltonian(param; β=0.0, Λ=1, Π=1)
    spbases = make_spbases(param)
    @time Hmat = make_single_particle_Hamiltonian(param, spbases, β, Λ, Π)
    @time eigvals(Hmat)
end



function calc_n_lj(l,j) 
    2l + (l===0) + (j===2l-1)
end

function test_calc_n_lj(lmax)
    n = 0
    for l in 0:lmax, j in 2l+1: -2: max(2l-1,0)
        n += 1
        n_lj = calc_n_lj(l,j)

        println("")
        @show n n_lj (n==n_lj)
    end
end

function calc_single_particle_states(param, spbases, β)
    @unpack Nr, Δr, lmax, Λmax = param
    @unpack nbases = spbases 
    
    nstates_max = nbases*2*cld(Λmax,2)
    
    coeffs = zeros(Float64, nbases, nstates_max)
    ψs = zeros(Float64, Nr, 2lmax+1, nstates_max)
    spEs = zeros(Float64, nstates_max)
    qnums = Vector{QuantumNumbers}(undef, nstates_max)
    occ = zeros(Float64, nstates_max)
    
    nstates = 0
    for Λ in 1:2:Λmax, Π in 1: -2: -1
        qnum = QuantumNumbers(Λ=Λ, Π=Π)
        
        Hmat = make_single_particle_Hamiltonian(param, spbases, β, Λ, Π)
        vals, vecs = eigen(Hmat)
        
        for ival in 1:length(vals)
            nstates += 1
            
            ibasis_active = 0
            for ibasis in 1:nbases
                @unpack l, j = spbases.qnums[ibasis]
                n_lj = calc_n_lj(l,j)

                if j < abs(Λ) || (-1)^l ≠ Π
                    continue
                end 
                ibasis_active += 1
                
                coeffs[ibasis, nstates] = vecs[ibasis_active, ival]
                @views @. ψs[:, n_lj, nstates] += spbases.ψs[:, ibasis]*coeffs[ibasis, nstates]
                spEs[nstates] = vals[ival]
                qnums[nstates] = qnum
            end
        end
    end
    p = sortperm(spEs[1:nstates])
    
    spstates = SingleParticleStates(
        nstates=nstates, 
        coeffs=coeffs[:,p], 
        ψs=ψs[:,:,p],
        spEs=spEs[p], 
        qnums=qnums[p], 
        occ=occ[p]
    )
end

function calc_occ!(spstates, param)
    @unpack N = param
    @unpack nstates, occ = spstates
    
    fill!(occ, 0)
    n = 0
    for i in 1:nstates
        if n + 2 ≤ N
            occ[i] = 1
            n += 2
        elseif n < N
            occ[i] = (N - n)/2
            n = N
        end
    end
    @assert n === N
    return
end

function show_spstates(spstates; Emax_show=0.0)
    @unpack nstates, spEs, qnums, occ = spstates 
    println("")
    for i in 1:nstates
        if spEs[i] > Emax_show 
            break 
        end
        println("i = ", i, ": ")
        @show spEs[i] occ[i] qnums[i]
    end
end

function plot_spstates(param, spstates, istate)
    @unpack rs, lmax = param
    @unpack ψs = spstates 
    p = plot(xlabel="r [fm]", ylabel="ψ_lj")
    for l in 0:lmax, j in 2l+1: -2: max(2l-1,0)
        n_lj = calc_n_lj(l,j)
        plot!(p, rs, ψs[:, n_lj, istate]; label="l=$l, j=$(j//2)")
    end
    display(p)
end


function test_calc_single_particle_states(param; β=0.0, Emax_show=1.0, istate=1)
    spbases = make_spbases(param)
    @time spstates = calc_single_particle_states(param, spbases, β)
    calc_occ!(spstates, param)
    show_spstates(spstates; Emax_show=Emax_show)
    #plot_spstates(param, spstates, istate)
end



function plot_nilsson_diagram(param; β_max=0.4, β_min=-0.4, Δβ=0.05)
    @unpack Z, N, Λmax = param

    spbases = make_spbases(param)
    @unpack nbases = spbases

    βs = range(β_min, β_max; step=Δβ)
    Nβ = length(βs)

    spEss = zeros(Float64, nbases, Nβ)

    p = plot(ylim=(-25,5), legend=false, 
    xlabel="β", ylabel="single-particle energy [MeV]", 
    title="nilsson diagram for Z=$Z, N=$N")

    for Λ in 1:2:Λmax, Π in 1: -2: -1
        spEss .= NaN
        for iβ in 1:Nβ
            β = βs[iβ]
            Hmat = make_single_particle_Hamiltonian(param, spbases, β, Λ, Π)
            vals = eigvals(Hmat)
            spEss[1:length(vals), iβ] = vals
        end
        plot!(βs, spEss')
    end
    display(p)
end



function plot_total_energy(param; β_max=0.4, β_min=-0.4, Δβ=0.05)
    @unpack Z, N, Λmax = param

    spbases = make_spbases(param)
    @unpack nbases = spbases

    βs = range(β_min, β_max; step=Δβ)
    Nβ = length(βs)

    Es_tot = zeros(Float64, Nβ)
    for iβ in 1:Nβ 
        β = βs[iβ]

        spstates = calc_single_particle_states(param, spbases, β)
        calc_occ!(spstates, param) 
        @unpack nstates, spEs, occ = spstates 

        for i in 1:nstates
            Es_tot[iβ] += 2occ[i]*spEs[i]
        end
        Es_tot[iβ] /= N 
    end
    p = plot(ylim=(-20, 0), legend=false, 
    xlabel="β", ylabel="total energy per nucleon [MeV]", title="Z=$Z, N=$N")
    plot!(p, βs, Es_tot)
    display(p)
end