struct KippReactionRate{TT<:Real} <: ReactionRates.AbstractReactionRate
    name::Symbol
    # num_iso_in of isotopes iso_in are converted into num_iso_out of isotopes iso_out
    iso_in::Vector{Symbol}
    num_iso_in::Vector{Int}
    iso_out::Vector{Symbol}
    num_iso_out::Vector{Int}
    Qvalue::TT  # energy released per reaction of this type (i.e. the mass defect), in erg
end

reaction_list[:kipp_rates] = Dict(
    :kipp_pp => KippReactionRate(:kipp_pp, [:H1], [4], [:He4], [1],
                                 ((4 * Chem.isotope_list[:H1].mass - Chem.isotope_list[:He4].mass) * AMU * CLIGHT^2)),
    :kipp_cno => KippReactionRate(:kipp_cno, [:H1], [4], [:He4], [1],
                                  ((4 * Chem.isotope_list[:H1].mass - Chem.isotope_list[:He4].mass) * AMU * CLIGHT^2)),
    :kipp_3alphaCF88 => KippReactionRate(:kipp_3alphaCF88, [:He4], [3], [:C12], [1],
                                         ((3 * Chem.isotope_list[:He4].mass - Chem.isotope_list[:C12].mass) * AMU *
                                          CLIGHT^2)),
    :kipp_3alphaA99 => KippReactionRate(:kipp_3alphaA99, [:He4], [3], [:C12], [1],
                                        ((3 * Chem.isotope_list[:He4].mass - Chem.isotope_list[:C12].mass) * AMU *
                                         CLIGHT^2)),
    :kipp_C12alpha => KippReactionRate(:kipp_C12alpha, [:C12, :He4], [1, 1], [:O16], [1],
                                       ((1 * Chem.isotope_list[:He4].mass + 1 * Chem.isotope_list[:C12].mass -
                                         Chem.isotope_list[:O16].mass)
                                        * AMU * CLIGHT^2)),
    :kipp_O16alpha => KippReactionRate(:kipp_O16alpha, [:O16, :He4], [1, 1], [:Ne20], [1],
                                       ((Chem.isotope_list[:He4].mass + Chem.isotope_list[:O16].mass -
                                         Chem.isotope_list[:Ne20].mass)
                                        * AMU * CLIGHT^2)),
    :kipp_CC => KippReactionRate(:kipp_CC, [:C12], [2], [:O16, :He4], [1, 2],
                                 ((2 * Chem.isotope_list[:C12].mass - Chem.isotope_list[:O16].mass -
                                   2 * Chem.isotope_list[:He4].mass)
                                  * AMU * CLIGHT^2)),
    :kipp_OO => KippReactionRate(:kipp_OO, [:O16], [2], [:Mg24, :He4], [1, 2],
                                 ((2 * Chem.isotope_list[:O16].mass - Chem.isotope_list[:Mg24].mass -
                                   2 * Chem.isotope_list[:He4].mass)
                                  * AMU * CLIGHT^2)),
)

"""
    function get_reaction_rate(reaction::KippReactionRate, T::T1, ρ::T2, xa::AbstractVector{TT},
                                xa_index::Dict{Symbol,Int})::TT where {TT,T1,T2}

Input:
reaction: the reaction to evaluate for
T: the temperature
ρ: the density
xa: element mass fractions
xa_index: index of the elements

Output:
ϵ_nuc / Qvalue, has units s^-1 g^-1
"""
function get_reaction_rate(reaction::KippReactionRate, T::T1, ρ::T2, xa::AbstractVector{TT},
                           xa_index::Dict{Symbol,Int})::TT where {TT,T1,T2}
    ϵnuc::TT = 0
    if reaction.name == :kipp_pp
        phi = 1
        f_11 = 1
        T9 = T / 1e9
        X1 = xa[xa_index[:H1]]

        g_11 = (1 + 3.82 * T9 + 1.51 * T9^2 + 0.144 * T9^3 - 0.0114 * T9^4)
        ϵnuc = 2.57e4 * phi * f_11 * g_11 * ρ * X1^2 * cbrt(T9^(-2)) * exp(-3.381 * cbrt(T9^(-1)))

    elseif reaction.name == :kipp_cno
        T9 = T / 1e9
        X1 = xa[xa_index[:H1]]
        X_CNO = xa[xa_index[:C12]] + xa[xa_index[:N14]] + xa[xa_index[:O16]]

        g_14 = (1 - 2.00 * T9 + 3.41 * T9^2 - 2.43 * T9^3)
        ϵnuc = 8.24e25 * g_14 * X_CNO * X1 * ρ * cbrt(T9^(-2)) * exp(-15.231 * cbrt(T9^(-1)) - (T9 / 0.8)^2)

    elseif reaction.name == :kipp_3alphaCF88
        f_3alpha = 1
        X4 = xa[xa_index[:He4]]
        T8 = T / 1e8

        ϵnuc = 5.09e11 * f_3alpha * ρ^2 * X4^3 * T8^-3 * exp(-44.027 / T8)

    elseif reaction.name == :kipp_3alphaA99
        f_3alpha = 1
        X4 = xa[xa_index[:He4]]
        T9 = T / 1e9

        ϵnuc = 6.272 * ρ^2 * X4^3 * (1 + 0.0158 * T9^(-0.65)) *
               (2.43e9 * cbrt(T9^(-2)) * exp(-13.490 * cbrt(T9^(-1)) - (T9 / 0.15)^2) * (1 + 74.5 * T9) +
                6.09e5 * sqrt(T9^(-3)) * exp(-1.054 / T9)) *
               (2.76e7 * cbrt(T9^(-2)) * exp(-23.570 * cbrt(T9^(-1)) - (T9 / 0.4)^2)
                * (1 + 5.47 * T9 + 326 * T9^2) + 130.7 * sqrt(T9^(-3)) * exp(-3.338 / T9)
                + 2.51e4 * sqrt(T9^(-3)) * exp(-20.307 / T9))

    elseif reaction.name == :kipp_C12alpha
        f_12alpha = 1
        X4 = xa[xa_index[:He4]]
        X12 = xa[xa_index[:C12]]
        T8 = T / 1e8

        ϵnuc = 1.3e27 * f_12alpha * ρ * X4 * X12 * T8^(-2) *
               ((1 + 0.134 * cbrt(T8^(2))) / (1 + 0.017 * cbrt(T8^(2))))^2 * exp(-69.20 / cbrt(T8))

    elseif reaction.name == :kipp_O16alpha
        f_16alpha = 1
        X4 = xa[xa_index[:He4]]
        X16 = xa[xa_index[:O16]]
        T9 = T / 1e9

        # ϵnuc = 1.91e27 * cbrt(T9^(-2)) * X16 * X4 * ρ * f_16alpha *
        #        exp(-39.76 * T9^(-1/3) - (T9/1.6)^2) 
        #        + 3.64 * 10^18 * sqrt(T9^(-3)) * exp(-10.32 / T9)
        #        + 4.39 * 10^19 * sqrt(T9^(-3)) * exp(-12.20 / T9)
        #        + 2.92 * 10^16 * T9^(2.966) * exp(-11.90 / T9)
        # Pablo: I think the above one has a typo in the Kippenhahn book itself,
        # otherwise it does not make sense that part of the rate is independent
        # of density and mass fractions. I think the one below would make sense.
        ϵnuc = X16 * X4 * ρ * f_16alpha *
               (1.91e27 * cbrt(T9^(-2)) * exp(-39.76 * T9^(-1 / 3) - (T9 / 1.6)^2)
                + 3.64e18 * sqrt(T9^(-3)) * exp(-10.32 / T9)
                + 4.39e19 * sqrt(T9^(-3)) * exp(-12.20 / T9)
                + 2.92e16 * T9^(2.966) * exp(-11.90 / T9))

    elseif reaction.name == :kipp_CC
        f_CC = 1
        T9 = T / 1e9
        T_9a = T9 / (1 + 0.0396 * T9)
        X12 = xa[xa_index[:C12]]

        ϵnuc = 1.86e43 * f_CC * ρ * X12^2 * sqrt(T9^(-3)) * T_9a^(5 / 6) * exp(-84.165 / cbrt(T_9a) - 2.12e-3 * T9^3)

    elseif reaction.name == :kipp_OO
        f_OO = 1
        T9 = T / 1e9
        X16 = xa[xa_index[:O16]]

        exp_func = -135.93 / cbrt(T9) - 0.629 * cbrt(T9^(2)) - 0.445 * cbrt(T9^(4)) + 0.0103 * T9^2
        ϵnuc = 2.14e53 * f_OO * ρ * X16^2 * cbrt(T9^(-2)) * exp(exp_func)

    else
        throw(ArgumentError("No method to compute rate for $(reaction.name)"))
    end
    return ϵnuc / reaction.Qvalue
end