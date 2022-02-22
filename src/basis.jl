
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
    @unpack M, Nr, Δr, rs, V₀, r₀, R₀, a, κ = param
    @unpack l, j = qnum
    
    # potential
    Vs = zeros(Float64, Nr)
    
    # central part
    @. Vs += V₀/(1 + exp((rs-R₀)/a))
    
    # spin-orbit part
    ls = (j*(j+2) - 4l*(l+1) - 3)/8
    @. Vs += κ*V₀*ls*(r₀*r₀/(a*rs))*exp((rs-R₀)/a)/(1 + exp((rs-R₀)/a))^2
    
    # centrifugal part
    @. Vs += M*l*(l+1)/rs^2
    
    
    dv = zeros(Float64, Nr)
    ev = zeros(Float64, Nr-1)
    
    @. dv = 2M/Δr^2 + Vs
    @. ev = -M/Δr^2
    
    return SymTridiagonal(dv, ev)
end

function make_spbases(param)
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
    
    spbases = SingleParticleBases(nbases, ψs[:,p], spEs[p], qnums[p])
end


function test_make_spbases(param)
    @time spbases = make_spbases(param)
    spbases.spEs
end