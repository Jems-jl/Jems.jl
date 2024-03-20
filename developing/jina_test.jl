using Jems.ReactionRates
using Jems.Chem
using Jems.EOS
using Jems.Constants
using CairoMakie

##
include("../src/ReactionRates/JinaRates.jl")

jinarates = Dict()
ref_dict = Dict()
read_set(file_contents, jinarates, ref_dict)
kipprates = ReactionRates.reaction_list[:kipp_rates]

##
logts = LinRange(7.0, 9.0, 30)
eosses = [EOS.EOSResults{Float64}() for i in 1:length(logts)]
for i in eachindex(logts)
    eosses[i].ρ = 1
    eosses[i].T = 10^logts[i]
end
xa = [0.2, 0.2, 0.2, 0.2, 0.2]
xa_index = Dict(:H1 => 1, :He4 => 2, :C12 => 3, :O16 => 4, :N14 => 5)

##  choose a reaction
# j_reactions = [jinarates[:He4_O16_to_Ne20],
#              jinarates[:He4_O16_to_Ne20_co10_n_x],
#              jinarates[:He4_O16_to_Ne20_co10_r_x]
#              ]
# k_reaction = kipprates[:kipp_O16alpha]
# j_reactions = [jinarates[:He4_C12_to_O16],
#                jinarates[:He4_C12_to_O16_nac2_x_x]
#                ]
# k_reaction = kipprates[:kipp_C12alpha]
# j_reactions = [jinarates[:He4_He4_He4_to_C12],
#                jinarates[:He4_He4_He4_to_C12_fy05_n_x],
#                jinarates[:He4_He4_He4_to_C12_fy05_r_x]]
# k_reaction = kipprates[:kipp_3alphaA99]

rates = zeros(length(logts))
for i in eachindex(logts)
    for reaction in j_reactions
        rates[i] += get_reaction_rate(reaction, eosses[i], xa, xa_index)
    end
end
j_rates = log10.(rates)
for i in eachindex(logts)
    rates[i] = ReactionRates.get_reaction_rate(k_reaction, eosses[i], xa, xa_index)
end
k_rates = log10.(rates)
## 
f = Figure();
ax = Axis(f[1, 1], xlabel=L"\log T/\textrm{K}", ylabel=L"\log R / {\textrm{g}}^{-1} {\textrm{s}}^{-1}");
lines!(ax, logts, j_rates, label="jina")
lines!(ax, logts, k_rates, label="Kipp")
axislegend(position=:lt)
f