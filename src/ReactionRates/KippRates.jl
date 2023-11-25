struct KippReactionRate{TT<:Real}<:ReactionRates.AbstractReactionRate
    name::Symbol
    iso_in::Vector{Symbol}
    num_iso_in::Vector{Int64}
    iso_out::Vector{Symbol}
    num_iso_out::Vector{Int64}
    Qvalue::TT
end

reaction_list[:kipp_rates] = Dict(
    :kipp_pp => KippReactionRate(:kipp_pp, [:H1], [4], [:He4], [1],
        ((4 * Chem.isotope_list[:H1].mass - Chem.isotope_list[:He4].mass) * AMU * CLIGHT^2)),
        
    :kipp_cno => KippReactionRate(:kipp_cno, [:H1], [4], [:He4], [1],
        ((4 * Chem.isotope_list[:H1].mass - Chem.isotope_list[:He4].mass) * AMU * CLIGHT^2)),

    :kipp_3alpha => KippReactionRate(:kipp_3alpha, [:He4], [3], [:C12], [1],
        ((3 * Chem.isotope_list[:He4].mass - Chem.isotope_list[:C12].mass) * AMU * CLIGHT^2)),

    :kipp_12alpha => KippReactionRate(:kipp_12alpha, [:C12, :He4], [1,1], [:O16], [1],
        ((1 * Chem.isotope_list[:He4].mass + 1 * Chem.isotope_list[:C12].mass - Chem.isotope_list[:O16].mass) * AMU * CLIGHT^2)),

    :kipp_16alpha => KippReactionRate(:kipp_16alpha, [:O16, :He4], [1,1], [:Ne20], [1],
        ((Chem.isotope_list[:He4].mass + Chem.isotope_list[:O16].mass - Chem.isotope_list[:Ne20].mass) * AMU * CLIGHT^2)),

    :kipp_CC => KippReactionRate(:kipp_CC, [:C12], [2], [:O16, :He4], [1,2],
        ((2 * Chem.isotope_list[:C12].mass - Chem.isotope_list[:O16].mass - 2 * Chem.isotope_list[:He4].mass) * AMU * CLIGHT^2)),

    :kipp_OO => KippReactionRate(:kipp_OO, [:O16], [2], [:Mg24, :He4], [1,2],
        ((2 * Chem.isotope_list[:O16].mass - Chem.isotope_list[:Mg24].mass - 2 * Chem.isotope_list[:He4].mass) * AMU * CLIGHT^2))   
)


function get_reaction_rate(reaction::KippReactionRate, eos00::EOSResults{TT}, xa::AbstractVector{TT}, xa_index::Dict{Symbol,Int})::TT where{TT}
    
    """
    Input:
    reaction: the reactions dictionary that is being used
    eos00: results equation of state
    xa: element fractions in the star
    xa_index: index of the elements

    Output:
    for each reaction, the ϵnuc value is calculated in the function.
        
    """
    
    if reaction.name == :kipp_pp 

        chi  = 1
        f_11 = 1
        T9   = (eos00.T / 1e9)
        X1   = xa[xa_index[:H1]]

        g_11 = (1 + 3.82 * T9 + 1.51 * T9^2 + 0.144 * T9^3 - 0.0114 * T9^4)
        ϵnuc = 2.57 * 10^4 * chi * f_11 * g_11 * eos00.ρ * X1^2 * T9^(-2/3) * exp(-3.381 * T9^(-1/3))

        return ϵnuc / reaction.Qvalue

    elseif reaction.name == :kipp_cno 

        T9    = (eos00.T / 1e9)
        X1    = xa[xa_index[:H1]]
        X_CNO = xa[xa_index[:C12]] + xa[xa_index[:N14]] + xa[xa_index[:O16]]

        g_14  = (1 - 2 * T9 + 3.41 * T9^2 - 2.43 * T9^3 )
        ϵnuc  = 8.24 * 10^(25) * g_14 * X_CNO * X1 * eos00.ρ * 
                T9^(-2/3) * exp(-15.231 * T9^(-1/3) - (T9/0.8)^2)

        return ϵnuc / reaction.Qvalue

    elseif reaction.name == :kipp_3alpha

        f_3alpha = 1
        X4   = xa[xa_index[:He4]]
        T8   = eos00.T / 1e8

        ϵnuc = 5.09 * 10^(11) * f_3alpha * (eos00.ρ)^2 * X4^3 * 
               T8^-3 * exp(-44.027 / T8)

    elseif reaction.name == :kipp_12alpha

        f_12alpha = 1
        X4   = xa[xa_index[:He4]]
        X12  = xa[xa_index[:C12]]
        T8   = eos00.T / 1e8
        
        ϵnuc = 1.3 * 10^(27) * f_12alpha * eos00.ρ * X4 * X12 * T8^(-2) *
                ((1 + 0.134 * T8^(2/3))/(1 + 0.017 * T8^(2/3)))^2 * exp(-69.20 / T8^(1/3))

    elseif reaction.name == :kipp_12alpha

        f_16alpha = 1
        X4   = xa[xa_index[:He4]]
        X16  = xa[xa_index[:016]]
        T9   = eos00.T / 1e9
                
        ϵnuc = 1.91 * 10^(27) * T9^(-2/3) * X16 * X4 * eos00.ρ * f_16alpha *
                exp(-39.76 * T9^(-1/3) - (T9/1.6)^2) 
                + 3.64 * 10^18 * T9^(-2/3) * exp(-10.32 / T9)
                + 4.39 * 10^19 * T9^(-2/3) * exp(-12.20 / T9)
                + 2.92 * 10^16 * T9^(2.966) * exp(-11.90 / T9)


    elseif reaction.name == :kipp_CC

        f_CC = 1
        T9   = (eos00.T / 1e9)
        T_9a = T9 / (1 + 0.0396 * T9) 
        X12  = xa[xa_index[:C12]]

        ϵnuc = 1.86 * 10^43 * f_CC * eos00.ρ * X12^2 * T9^(-3/2) * T_9a^(5/6) *
               exp(-84.165 / T_9a^(1/3) - 2.12 * 10^(-3) * T9^3)

    elseif reaction.name == :kipp_OO

        f_OO = 1
        T9   = (eos00.T / 1e9)
        X16  = xa[xa_index[:O16]]

        exp_func = -135.93 / (T9^(1/3)) - 0.629 * T9^(2/3) - 0.445 * T9^(4/3) + 0.0103 * T9^2
        ϵnuc = 2.14 * 10^53 * f_OO * eos00.ρ * X16^2 * T9^(-2/3) * exp(exp_func)

    else
        throw(ArgumentError("No method to compute rate for $(reaction.name)"))
    end
end