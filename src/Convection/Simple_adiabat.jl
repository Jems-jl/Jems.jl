function simple_adiabatic(mlt::AdiabaticMLT, )
    κface = exp((sm.dm[k] * log(κ00) + sm.dm[k + 1] * log(κp1)) / (sm.dm[k] + sm.dm[k + 1]))
    L₀ = var00[sm.vari[:lum]] * LSUN
    r₀ = exp(var00[sm.vari[:lnr]])
    Pface = exp((sm.dm[k] * var00[sm.vari[:lnP]] + sm.dm[k + 1] * varp1[sm.vari[:lnP]]) /
                (sm.dm[k] + sm.dm[k + 1]))
    lnT₊ = varp1[sm.vari[:lnT]]
    lnT₀ = var00[sm.vari[:lnT]]
    Tface = exp((sm.dm[k] * lnT₀ + sm.dm[k + 1] * lnT₊) / (sm.dm[k] + sm.dm[k + 1]))
    ∇ᵣ = 3κface * L₀ * Pface / (16π * CRAD * CLIGHT * CGRAV * sm.m[k] * Tface^4)
    ∇ₐ = (sm.dm[k] * eos00[7] + sm.dm[k + 1] * eosp1[7]) / (sm.dm[k] + sm.dm[k + 1])

    if (∇ᵣ < ∇ₐ)
        return ∇ᵣ
    else
        return ∇ₐ
    end
end