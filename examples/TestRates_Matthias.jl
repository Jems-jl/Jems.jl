using Jems.ReactionRates
using Jems.Chem
using Jems.EOS
using Jems.Constants
using CairoMakie

##
include("../src/ReactionRates/JinaRates.jl")

jinarates       = ReactionRates.jina_rates
jinareferences  = ReactionRates.jina_references
# ReactionRates.read_set(file_contents, jinarates, ref_dict)
kipprates = ReactionRates.reaction_list[:kipp_rates]

##

logts = LinRange(6.0, 9.0, 30)
eosses = [EOS.EOSResults{Float64}() for i in 1:length(logts)]
for i in eachindex(logts)
    eosses[i].ρ = 1
    eosses[i].T = 10^logts[i]
end
xa = [0.0, 1.0, 0.0, 0.0, 0.0]
xa_index = Dict(:H1 => 1, :He4 => 2, :C12 => 3, :O16 => 4, :N14 => 5)

##  

# choose a reaction

### O16 α ###

# j_reactions = [jinarates[:He4_O16_to_Ne20_co10_r_x_0],
#                jinarates[:He4_O16_to_Ne20_co10_r_x_1],
#                jinarates[:He4_O16_to_Ne20_co10_n_x_0]
#              ]
# k_reaction = kipprates[:kipp_O16alpha]

### C12 α ###

# j_reactions = [jinarates[:He4_C12_to_O16_nac2_x_x_0],
#                jinarates[:He4_C12_to_O16_nac2_x_x_1]
#                ]
# k_reaction = kipprates[:kipp_C12alpha]

### triple α ###

j_reactions = [# jinarates[:He4_He4_He4_to_C12_fy05_r_x_0],
               jinarates[:He4_He4_He4_to_C12_fy05_r_x_1],
               jinarates[:He4_He4_He4_to_C12_fy05_n_x_0]]
              
               
k_reaction = kipprates[:kipp_3alphaA99]

### PP ###

# j_reactions = [jinarates[:H1_H1_to_D2_betplus_w_x_0],
#                jinarates[:H1_H1_to_D2_xxec_w_x_0]
#                ]
# k_reaction = kipprates[:kipp_pp]


### CNO ###

# j_reactions = [jinarates[:H1_N14_to_O15_im05_r_x_0],
#                jinarates[:H1_N14_to_O15_im05_n_x_0],
#                jinarates[:H1_N14_to_O15_im05_n_x_1],
#                jinarates[:H1_N14_to_O15_im05_r_x_1]
#                ]

# k_reaction = kipprates[:kipp_cno]

##

rates = zeros(length(logts))
# println(rates)
for i in eachindex(logts)
    for reaction in j_reactions
        rates[i] += ReactionRates.get_reaction_rate(reaction, eosses[i], xa, xa_index) # * Constants.AVO^(-1)
    end
end
# println(rates)
j_rates_old = rates.*(Constants.AVO)
j_rates     = log10.(j_rates_old)


##

function angulo_pp(eosr, xa, xa_index)
    phi = 1
    f_11 = 1
    T9 = (eosr.T / 1e9)
    X1 = xa[xa_index[:H1]]

    g_11 = (1 + 3.82 * T9 + 1.51 * T9^2 + 0.144 * T9^3 - 0.0114 * T9^4)
    nasigmav = 4.08e-15 * f_11 * phi * g_11 * cbrt(T9^(-2)) * exp(-3.381 * cbrt(T9^(-1)))
    return nasigmav * eosr.ρ * 
            (X1 / (Chem.isotope_list[:H1].mass * Constants.AMU))^2 /
            Constants.AVO / 2
end

function angulo_cno(eosr, xa, xa_index)

    T9     = (eosr.T / 1e9)
    X14    = xa[xa_index[:N14]]
    X1     = xa[xa_index[:H1]]
    # X_CNO = xa[xa_index[:C12]] + xa[xa_index[:N14]] + xa[xa_index[:O16]]

    g_14  = (1 - 2.00 * T9 + 3.41 * T9^2 - 2.43 * T9^3 )
    nasigmav  = 4.83e7 * g_14 * cbrt(T9^(-2)) * exp(-15.231 * cbrt(T9^(-1)) - (T9/0.8)^2)
                + 2.36e3 * T9^(-3/2) * exp(-3.010 / T9)
                + 6.72e3 * T9^(0.380) * exp(-9530 / T9)

    return nasigmav * eosr.ρ * 
            (X14 / (Chem.isotope_list[:N14].mass * Constants.AMU)) *
            (X1  / (Chem.isotope_list[:H1].mass  * Constants.AMU)) /
            (Constants.AVO)

end


function angulo_3α(eosr, xa, xa_index)

    # f_3alpha = 1
    X4   = xa[xa_index[:He4]]
    T9   = eosr.T / 1e9

    fact1 = 2.76e7 * cbrt(T9^(-2)) * exp(-23.570*cbrt(T9^(-1)) - (T9/0.4)^2) *
            (1 + 5.47 * T9 + 326 * T9^2) + 130.7 * sqrt(T9^(-3)) * exp(-3.338/T9)
            + 2.51e4 * sqrt(T9^(-3)) * exp(-20.307/T9)

    
    fact2 = 2.43e9 * cbrt(T9^(-2)) * exp(-13.490 * cbrt(T9^(-1)) - (T9/0.15)^2) * (1 + 74.5 * T9)
            + 6.09e5 * sqrt(T9^(-3)) * exp(-1.054/T9)

    # nasigmav = 3.44e-16 * (1 + 0.0158 * T9^(-0.65)) * fact1 * fact2

    if T9 <= 0.03

        na2sigmav = 3.07e-16 * ( 1- 29.1 * T9 + 1308 * T9^2) * fact1 * fact2

    else

        na2sigmav = 3.44e-16 * (1 + 0.0158 * T9^(-0.65)) * fact1 * fact2

    end


    return na2sigmav * (eosr.ρ)^2 *
            (X4  / (Chem.isotope_list[:He4].mass  * Constants.AMU))^3 / 
            (Constants.AVO)^2 / 2  
    
end


for i in eachindex(logts)
    # rates[i] = ReactionRates.get_reaction_rate(k_reaction, eosses[i], xa, xa_index)
    # rates[i] = angulo_pp(eosses[i], xa, xa_index)  # matches exactly with Jina
    # rates[i] = angulo_cno(eosses[i], xa, xa_index)
    rates[i] = angulo_3α(eosses[i], xa, xa_index)

end
k_rates = log10.(rates)


##


f = Figure();
ax = Axis(f[1, 1], xlabel=L"\log T/\textrm{K}", ylabel=L"\log R / {\textrm{g}}^{-1} {\textrm{s}}^{-1}");
lines!(ax, logts, j_rates, label="jina")
lines!(ax, logts, k_rates, label="Kipp")
axislegend(position=:lt)
f



