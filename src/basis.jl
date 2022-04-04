
export test_make_spbases, test_integrate_SchEq!, test_calc_matching_condition

@with_kw struct BasisQuantumNumbers @deftype Int64
    l = 0
    j = 1
end

@with_kw struct SingleParticleBases
    nbases::Int64
    ψs::Matrix{Float64}; @assert size(ψs, 2) === nbases
    spEs::Vector{Float64}; @assert length(spEs) === nbases
    qnums::Vector{BasisQuantumNumbers}; @assert length(qnums) === nbases
end


"""
    make_basis_Hamiltonian(param, qnum)

Make Hamiltonian for single-particle bases.
"""
function make_basis_Hamiltonian(param, qnum)
    @unpack M, Nr, Δr, rs, V₀, V₁, r₀, R₀, a, V_gaus, μ_gaus = param
    @unpack l, j = qnum
    
    # potential
    Vs = zeros(Float64, Nr)
    
    # central part
    @. Vs += V₀/(1 + exp((rs-R₀)/a))
    
    # spin-orbit part
    ls = (j*(j+2) - 4l*(l+1) - 3)/8
    @. Vs += V₁*ls*(r₀*r₀/(a*rs))*exp((rs-R₀)/a)/(1 + exp((rs-R₀)/a))^2

    if l === 0
        @. Vs += V_gaus*exp(-μ_gaus*rs*rs)
    end
    
    # centrifugal part
    @. Vs += M*l*(l+1)/rs^2
    
    dv = zeros(Float64, Nr)
    ev = zeros(Float64, Nr-1)
    
    @. dv = 2M/Δr^2 + Vs
    @. ev = -M/Δr^2
    
    return SymTridiagonal(dv, ev)
end

function make_spbases2(param)
    @unpack Nr, Δr, Emax, lmax = param 
    
    ψs = zeros(Float64, Nr, Nr*(lmax+1))
    spEs = zeros(Float64, Nr*(lmax+1))
    qnums = Vector{BasisQuantumNumbers}(undef, Nr*(lmax+1))
    
    nbases = 0
    for l in 0:lmax, j in (2l+1): -2: max(2l-1, 0)
        qnum = BasisQuantumNumbers(l=l,j=j)
        
        Hmat = make_basis_Hamiltonian(param, qnum)
        vals, vecs = eigen(Hmat)
        #@show dot(vecs[:,1], vecs[:,1])
        
        @. vecs /= sqrt(Δr)
        #@show minimum(vals), maximum(vals)
        
        #@show dot(vecs[:,1], vecs[:,1])*Δr
        
        for i in 1:Nr
            if vals[i] > Emax
                break
            end
            nbases += 1
            ψs[:,nbases] = vecs[:,i]
            spEs[nbases] = vals[i]
            qnums[nbases] = qnum
        end
    end
    
    p = sortperm(spEs[1:nbases])
    
    SingleParticleBases(nbases, ψs[:,p], spEs[p], qnums[p])
end


function test_make_spbases(param)
    @time spbases = make_spbases(param)
    spbases.spEs
end



function calc_basis_potential!(Vs, param, qnum)
    @unpack M, Nr, Δr, rs, V₀, V₁, r₀, R₀, a, V_gaus, μ_gaus = param
    @unpack l, j = qnum

    fill!(Vs, 0.0)
    
    # central part
    @. Vs += V₀/(1 + exp((rs-R₀)/a))
    
    # spin-orbit part
    ls = (j*(j+2) - 4l*(l+1) - 3)/8
    @. Vs += V₁*ls*(r₀*r₀/(a*rs))*exp((rs-R₀)/a)/(1 + exp((rs-R₀)/a))^2

    # gaussian potential for s-wave
    if l === 0
        @. Vs += V_gaus*exp(-μ_gaus*rs*rs)
    end
    
    # centrifugal part
    @. Vs += M*l*(l+1)/rs^2

    return 
end

function normalize!(ys, rs)
    Δr = rs[2] - rs[1]
    norm = 0.0
    for ir in 1:length(rs)
        norm += ys[ir]^2
    end
    norm *= Δr 
    norm = sqrt(norm) 
    @. ys /= norm 
    return 
end

function integrate_SchEq!(ys, As, param, Vs, E, qnum)
    @unpack M, Nr, Δr, rs, R₀, ir_matching = param 
    @unpack l, j = qnum

    fill!(ys, 0.0)
    @. As = 1 - (Δr*Δr/12)*(Vs-E)/M

    # outward integration
    ys[1,1] = Δr^(l+1)
    G₀ = ifelse(l===1, -Δr*Δr/6, 0.0)
    G₁ = As[1]*ys[1,1]
    for ir in 1:ir_matching
        G₂ = 12ys[ir,1] - 10G₁ - G₀
        ys[ir+1,1] = As[ir+1]\G₂
        G₀ = G₁ 
        G₁ = G₂
    end

    @views normalize!(ys[:,1], rs)

    # inward integration 
    ys[Nr-1,2] = Δr 
    G₀ = 0.0
    G₁ = As[Nr-1]*ys[Nr-1,2]
    for ir in Nr-1: -1: ir_matching
        G₂ = 12ys[ir,2] - 10G₁ - G₀
        ys[ir-1,2] = As[ir-1]\G₂
        G₀ = G₁ 
        G₁ = G₂
    end

    @views normalize!(ys[:,2], rs)
    
    return 
end


function test_integrate_SchEq!(param; E=-20, l=0, j=1)
    @unpack Nr, rs = param

    qnum = QuantumNumbers(l=l, j=j)

    Vs = zeros(Float64, Nr)
    calc_basis_potential!(Vs, param, qnum)

    ys = zeros(Float64, Nr, 2)
    As = zeros(Float64, Nr)
    @time integrate_SchEq!(ys, As, param, Vs, E, qnum)

    p = plot()
    plot!(p, rs, ys[:,1]; label="outward")
    plot!(p, rs, ys[:,2]; label="inward")
    plot!(p, rs, Vs/10)
end


function calc_matching_condition(param, ys)
    @unpack Δr, ir_matching = param 

    ir = ir_matching 

    a =  ys[ir,1]
    b = -ys[ir,2]

    c =  (ys[ir+1,1] - ys[ir-1,1])/2Δr 
    d = -(ys[ir+1,2] - ys[ir-1,2])/2Δr 

    return a*d - b*c
end

function test_calc_matching_condition(param; lmax=3)
    @unpack Nr, rs, Δr = param 

    Es = range(-25, 5; step=0.1)
    fs = similar(Es)

    Vs = zeros(Float64, Nr)
    ys = zeros(Float64, Nr, 2)
    As = zeros(Float64, Nr)

    spEs = zeros(Float64, 100)
    
    p = plot(ylim=(-1, 1))
    nbases = 0
    @time for l in 0:lmax, j in 2l+1: -2: max(2l-1,0)
        qnum = QuantumNumbers(l=l,j=j)
        calc_basis_potential!(Vs, param, qnum)

        function f(E)
            integrate_SchEq!(ys, As, param, Vs, E, qnum)
            calc_matching_condition(param, ys)
        end

        @. fs = f(Es)
        for iE in 1:length(Es)-1 
            f₁ = fs[iE]
            f₂ = fs[iE+1]

            @views if f₁*f₂ < 0
                E = find_zero(f, (Es[iE], Es[iE+1]), Bisection())
                integrate_SchEq!(ys, As, param, Vs, E, qnum)
                nbases += 1
                spEs[nbases] = E
            end
        end

        plot!(p, Es, fs; label="l=$l, j=$j")
    end
    display(p)
    sort!(spEs[1:nbases])
end



function make_spbases(param)
    @unpack Nr, Δr, rs, ir_matching, Emax, lmax = param 

    Es = range(-25, Emax; step=0.1)
    
    ψs = zeros(Float64, Nr, Nr*(lmax+1))
    spEs = zeros(Float64, Nr*(lmax+1))
    qnums = Vector{BasisQuantumNumbers}(undef, Nr*(lmax+1))

    Vs = zeros(Float64, Nr)
    ys = zeros(Float64, Nr, 2)
    As = zeros(Float64, Nr)
    
    nbases = 0
    for l in 0:lmax, j in (2l+1): -2: max(2l-1, 0)
        qnum = BasisQuantumNumbers(l=l,j=j)
        calc_basis_potential!(Vs, param, qnum)

        function f(E)
            integrate_SchEq!(ys, As, param, Vs, E, qnum)
            calc_matching_condition(param, ys)
        end

        for iE in 1:length(Es)-1 
            f₁ = f(Es[iE])
            f₂ = f(Es[iE+1])

            @views if f₁*f₂ < 0
                E = find_zero(f, (Es[iE], Es[iE+1]), Bisection())
                integrate_SchEq!(ys, As, param, Vs, E, qnum)
                
                nbases += 1
                spEs[nbases] = E 
                qnums[nbases] = qnum 

                ratio = ys[ir_matching, 1]/ys[ir_matching, 2]
                @. ys[:,2] *= ratio 

                for ir in 1:ir_matching
                    ψs[ir,nbases] = ys[ir,1]
                end
                for ir in ir_matching+1:Nr 
                    ψs[ir,nbases] = ys[ir,2]
                end

                normalize!(ψs[:,nbases], rs)
            end
        end
    end
    
    p = sortperm(spEs[1:nbases])
    
    SingleParticleBases(nbases, ψs[:,p], spEs[p], qnums[p])
end