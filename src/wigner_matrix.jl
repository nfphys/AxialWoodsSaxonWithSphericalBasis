
export test_diagonalize_Jy, test_wigner_matrix

"""
    diagonalize_Jy(jmax)

Diagonalize Jy for j = 0, 1, ⋯, jmax
All spins are expressed as double their actual values.
"""
function diagonalize_Jy(jmax)
    vals_Jy = zeros(Float64, jmax+1, jmax+1)
    vecs_Jy = zeros(ComplexF64, jmax+1, jmax+1, jmax+1)
    
    Jy = zeros(ComplexF64, jmax+1, jmax+1)
    
    @views for j in 0:jmax
        n = 0
        for m in j-2: -2: -j
            n += 1
            Jy[n,n+1] = sqrt(j*(j+2) - m*(m+2))/4im
            Jy[n+1,n] = -Jy[n,n+1]
        end
        
        vals, vecs = eigen(Jy[1:j+1,1:j+1])
        
        vals_Jy[1:j+1,j+1] = vals
        vecs_Jy[1:j+1,1:j+1,j+1] = vecs
    end
    
    return vals_Jy, vecs_Jy
    
    #=
    i = 0
    for m in j-2: -2: -j
        i += 1
        Jy[i,i+1] = sqrt(j*(j+2) - m*(m+2))/4im
        Jy[i+1,i] = -Jy[i,i+1]
    end
    eigen(Jy)
    =#
end


function test_diagonalize_Jy(;jmax=20)
    @time vals_Jy, vecs_Jy = diagonalize_Jy(jmax)
    vals_Jy
end


function wigner_matrix(vals_Jy, vecs_Jy, j, m, mp, θ)
    n  = div(j-m ,2) + 1
    np = div(j-mp,2) + 1
    
    d = 0.0
    for k in 1:j+1
        d += vecs_Jy[n,k,j+1]*
            exp(-im*θ*vals_Jy[k,j+1])*
            conj(vecs_Jy[np,k,j+1])
    end
    return real(d)
end


function test_wigner_matrix(j, m, mp, d_exact;jmax=20)
    @time vals_Jy, vecs_Jy = diagonalize_Jy(jmax)
    
    θs = range(0, 2π, length=100+1)
    Nθ = length(θs)
    
    ds = zeros(Float64, Nθ)
    ds_exact = zeros(Float64, Nθ)
    
    @time for iθ in 1:Nθ
        θ = θs[iθ]
        ds[iθ] = wigner_matrix(vals_Jy, vecs_Jy, j, m, mp, θ)
        ds_exact[iθ] = d_exact(θ)
    end
    
    err = sum(@. (ds - ds_exact)^2)
    @show err
    
    p = plot(xlim=(0,1), xlabel="θ/2π")
    plot!(p, θs/2π, ds)
    plot!(p, θs/2π, ds_exact)
    display(p)
end