
using ForwardDiff
using Base.Threads

function example_equationHSE(sm, k, varm1::Vector{<:TT}, var00::Vector{<:TT}, varp1::Vector{<:TT},
                            eosm1::Vector{<:TT}, eos00::Vector{<:TT}, eosp1::Vector{<:TT},
                            κm1::TT, κ00::TT, κp1::TT)::TT where{TT<:Real}
    if k==sm.nz
        lnP₀ = var00[sm.vari[:lnP]]
        return lnP₀ + 10.0 #force very low surface density to resemble zero pressure condition
    end
    if k==1
        lnP₀ = var00[sm.vari[:lnP]]
        lnP₊ = varp1[sm.vari[:lnP]]
        return lnP₀ - lnP₊#make it simple, just force first to cells to have constant pressure
    end
    lnP₋ = varm1[sm.vari[:lnP]]
    lnP₊ = varp1[sm.vari[:lnP]]
    P₀ = exp(var00[sm.vari[:lnP]])
    r₀ = exp(var00[sm.vari[:lnr]])
    dm = (sm.m[k+1]-sm.m[k-1])
    
    return (P₀*(lnP₊ - lnP₋)/dm + CGRAV*sm.m[k]/(4π*r₀^4))/(CGRAV*sm.m[k]/(4π*r₀^4))
end

function example_equationT(sm, k, varm1::Vector{<:TT}, var00::Vector{<:TT}, varp1::Vector{<:TT},
                          eosm1::Vector{<:TT}, eos00::Vector{<:TT}, eosp1::Vector{<:TT},
                          κm1::TT, κ00::TT, κp1::TT)::TT where{TT<:Real}
    #if k==1
    #    lnT₀ = v00[sm.vari[:lnT]]
    #    lnT₊ = vp1[sm.vari[:lnT]]
    #    return lnT₀ - lnT₊#make it simple, just force first to cells to have constant temperature
    #end
    if k==sm.nz
        lnT₀ = var00[sm.vari[:lnT]]
        return lnT₀ + 10.0 #force very low surface temperature
    end
    ρ = 1 # in cm/g
    ρ00 = eos00[1]
    return log(ρ)-log(ρ00) # T should be such that the EOS returns the density we want
end

function example_equationContinuity(sm, k, varm1::Vector{<:TT}, var00::Vector{<:TT}, varp1::Vector{<:TT},
                                   eosm1::Vector{<:TT}, eos00::Vector{<:TT}, eosp1::Vector{<:TT},
                                   κm1::TT, κ00::TT, κp1::TT)::TT where{TT<:Real}
    ρ = 1 # in cm/g
    
    if k==1
        lnr₀ = var00[sm.vari[:lnr]]
        return lnr₀ + 10.0 #force very low central radiusw
    end
    r₋ = exp(varm1[sm.vari[:lnr]])
    r₀ = exp(var00[sm.vari[:lnr]])
    dm = (sm.m[k]-sm.m[k-1])

    #expected_r₀ = r₋ + dm/(4π*r₋^2*ρ)
    expected_dr³_dm = 3/(4π*ρ)
    actual_dr³_dm = (r₀^3-r₋^3)/dm
    
    #return log(r₀) - log(expected_r₀)
    return (expected_dr³_dm - actual_dr³_dm)*ρ
end

function example_equationComposition(sm, k, varm1::Vector{<:TT}, var00::Vector{<:TT}, varp1::Vector{<:TT},
                                    eosm1::Vector{<:TT}, eos00::Vector{<:TT}, eosp1::Vector{<:TT},
                                    κm1::TT, κ00::TT, κp1::TT)::TT where{TT<:Real}
    return var00[sm.vari[:H1]] - 1.0
end