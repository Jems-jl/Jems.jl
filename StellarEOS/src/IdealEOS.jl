export get_EOS_resultsTP, get_EOS_resultsTρ

struct IdealEOS <: AbstractEOS
    include_radiation::Bool
end

function get_μ_IdealEOS(isotope_data::Dict{Symbol, Isotope}, xa::Vector{<:TT},species::Vector{Symbol})::TT where{TT<:Real}
    # See section 4.2 of Kipp
    μ::TT = 0
    for i in eachindex(species)
        μ += xa[i]*(1+isotope_data[species[i]].Z)/isotope_data[species[i]].mass
    end
    return 1/μ
end

function get_EOS_resultsTP(eos::IdealEOS, isotope_data::Dict{Symbol, Isotope},
        lnT::TT, lnP::TT, xa::Vector{<:TT},species::Vector{Symbol})::Vector{<:TT} where{TT<:Real}
    # See section 13.2 of Kipp
    β::TT = 1
    T = exp(lnT)
    P = exp(lnP)
    μ = get_μ_IdealEOS(isotope_data, xa, species)
    Prad = CRAD*T^4/3
    ρ = μ/(CGAS*T)*(P-Prad)

    if eos.include_radiation
        β = 1-Prad/P
        if β < 0
            throw(DomainError(β, "Input values produce negative beta"))
        end
    end
    α = 1/β
    δ = (4-3*β)/β

    u = CGAS*T/μ*(3/2+3*(1-β)/β)
    cₚ = CGAS/μ*(3/2 + 3*(4+β)*(1-β)/β^2 + δ*α)
    ∇ₐ = CGAS*δ/(β*μ*cₚ)
    Γ₁::TT = 1
    if eos.include_radiation
        Γ₁ = 1/(α-δ*∇ₐ)
    else
        Γ₁ = 1/(1-∇ₐ)
    end

    return [ρ, μ, β, u, cₚ, δ, ∇ₐ, Γ₁]
end

