function equationHSE!(results::AbstractVector{TT}, sm::StellarModel, k::Int,
                     varm1::AbstractVector{TT}, var00::AbstractVector{TT}, varp1::AbstractVector{TT},
                     eosm1::AbstractVector{TT}, eos00::AbstractVector{TT}, eosp1::AbstractVector{TT},
                     κm1::TT, κ00::TT, κp1::TT) where {TT<:Real}
    if k == sm.nz  # atmosphere boundary condition
        lnP₀ = var00[sm.vari[:lnP]]
        r₀ = exp(var00[sm.vari[:lnr]])
        g₀ = CGRAV * sm.mstar / r₀^2
        results[k] = lnP₀ - log(2g₀ / (3κ00))  # Eddington gray, ignoring radiation pressure term
        return
    end
    lnP₊ = varp1[sm.vari[:lnP]]
    lnP₀ = var00[sm.vari[:lnP]]
    lnPface = (sm.dm[k] * lnP₀ + sm.dm[k + 1] * lnP₊) / (sm.dm[k] + sm.dm[k + 1])
    r₀ = exp(var00[sm.vari[:lnr]])
    dm = (sm.m[k + 1] - sm.m[k])

    results[k] = (exp(lnPface) * (lnP₊ - lnP₀) / dm + CGRAV * sm.m[k] / (4π * r₀^4)) /
           (CGRAV * sm.m[k] / (4π * r₀^4))
end

function equationT!(results::AbstractVector{TT}, sm::StellarModel, k::Int,
                   varm1::AbstractVector{TT}, var00::AbstractVector{TT}, varp1::AbstractVector{TT},
                   eosm1::AbstractVector{TT}, eos00::AbstractVector{TT}, eosp1::AbstractVector{TT},
                   κm1::TT, κ00::TT, κp1::TT) where {TT<:Real}
    if k == sm.nz  # atmosphere boundary condition
        lnT₀ = var00[sm.vari[:lnT]]
        L₀ = var00[sm.vari[:lum]] * LSUN
        r₀ = exp(var00[sm.vari[:lnr]])
        results[k] = lnT₀ - log(L₀ / (BOLTZ_SIGMA * 4π * r₀^2)) / 4  # Eddington gray, ignoring radiation pressure term
        return
    end
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
        results[k] = (Tface * (lnT₊ - lnT₀) / sm.dm[k] +
                CGRAV * sm.m[k] * Tface / (4π * r₀^4 * Pface) * ∇ᵣ) / (CGRAV * sm.m[k] * Tface / (4π * r₀^4 * Pface))  # only radiative transport
    else  # should do convection here
        results[k] = (Tface * (lnT₊ - lnT₀) / sm.dm[k] +
                CGRAV * sm.m[k] * Tface / (4π * r₀^4 * Pface) * ∇ₐ) / (CGRAV * sm.m[k] * Tface / (4π * r₀^4 * Pface))  # only radiative transport
    end
end

function equationLuminosity!(results::AbstractVector{TT}, sm::StellarModel, k::Int,
                            varm1::AbstractVector{TT}, var00::AbstractVector{TT}, varp1::AbstractVector{TT},
                            eosm1::AbstractVector{TT}, eos00::AbstractVector{TT}, eosp1::AbstractVector{TT},
                            κm1::TT, κ00::TT, κp1::TT) where {TT<:Real}
    L₋::TT = 0  # central luminosity is zero at first cell
    if k > 1
        L₋ = varm1[sm.vari[:lum]] * LSUN  # change it if not at first cell
    end
    L₀ = var00[sm.vari[:lum]] * LSUN
    ρ₀ = eos00[1]
    cₚ = eos00[5]
    δ = eos00[6]
    dTdt = (exp(var00[sm.vari[:lnT]]) - exp(sm.ssi.lnT[k])) / sm.ssi.dt
    dPdt = (exp(var00[sm.vari[:lnP]]) - exp(sm.ssi.lnP[k])) / sm.ssi.dt

    ϵnuc = 0.1 * var00[sm.vari[:H1]]^2 * ρ₀ * (exp(var00[sm.vari[:lnT]]) / 1e6)^4 +
           0.1 * var00[sm.vari[:H1]] * ρ₀ * (exp(var00[sm.vari[:lnT]]) / 1e7)^18

    results[k] = ((L₀ - L₋) / sm.dm[k] - ϵnuc + cₚ * dTdt - (δ / ρ₀) * dPdt)  # no neutrinos
end

function equationContinuity!(results::AbstractVector{TT}, sm::StellarModel, k::Int,
                            varm1::AbstractVector{TT}, var00::AbstractVector{TT}, varp1::AbstractVector{TT},
                            eosm1::AbstractVector{TT}, eos00::AbstractVector{TT}, eosp1::AbstractVector{TT},
                            κm1::TT, κ00::TT, κp1::TT) where {TT<:Real}
    ρ₀ = eos00[1]
    r₀ = exp(var00[sm.vari[:lnr]])
    r₋::TT = 0  # central radius is zero at first cell
    if k > 1
        r₋ = exp(varm1[sm.vari[:lnr]])  # change it if not at first cell
    end

    dm = sm.m[k]  # this is only valid for k=1
    if k > 1
        dm = dm - sm.m[k - 1]
    end

    # expected_r₀ = r₋ + dm/(4π*r₋^2*ρ)
    expected_dr³_dm = 3 / (4π * ρ₀)
    actual_dr³_dm = (r₀^3 - r₋^3) / dm

    results[k] = (expected_dr³_dm - actual_dr³_dm) * ρ₀
end

# To test performance, include 8 isotopes similar to basic.net in MESA.
# of course we are keeping these fixed now, but it lets us test their impact on the
# computation of the jacobian
function equationH1!(results::AbstractVector{TT}, sm::StellarModel, k::Int,
                    varm1::AbstractVector{TT}, var00::AbstractVector{TT}, varp1::AbstractVector{TT},
                    eosm1::AbstractVector{TT}, eos00::AbstractVector{TT}, eosp1::AbstractVector{TT},
                    κm1::TT, κ00::TT, κp1::TT)::TT where {TT<:Real}
    ρ₀ = eos00[1]
    ϵnuc = 0.1 * var00[sm.vari[:H1]]^2 * ρ₀ * (exp(var00[sm.vari[:lnT]]) / 1e6)^4 +
           0.1 * var00[sm.vari[:H1]] * ρ₀ * (exp(var00[sm.vari[:lnT]]) / 1e7)^18
    rate_per_unit_mass = 4 * ϵnuc / ((4 * Chem.isotope_list[:H1].mass - Chem.isotope_list[:He4].mass) * AMU * CLIGHT^2)

    Xi = sm.ssi.ind_vars[(k - 1) * sm.nvars + sm.vari[:H1]]

    results[k] = (var00[sm.vari[:H1]] - Xi) / sm.ssi.dt +
           Chem.isotope_list[:H1].mass * AMU * rate_per_unit_mass
end

function equationHe4!(results::AbstractVector{TT}, sm::StellarModel, k::Int,
                     varm1::AbstractVector{TT}, var00::AbstractVector{TT}, varp1::AbstractVector{TT},
                     eosm1::AbstractVector{TT}, eos00::AbstractVector{TT}, eosp1::AbstractVector{TT},
                     κm1::TT, κ00::TT, κp1::TT)::TT where {TT<:Real}
    results[k] = var00[sm.vari[:He4]] + var00[sm.vari[:H1]] - 1.0
end
