
export Vnn_expectation_value

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