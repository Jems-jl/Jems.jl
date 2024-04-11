# Script to plot reaction rates easily
# Nog onder constructie

using Jems.ReactionRates
using Jems.Chem
using Jems.EOS
using Jems.Constants
using CairoMakie

##

toyrates  = ReactionRates.reaction_list[:toy_rates]
jinarates = ReactionRates.reaction_list[:jina_rates]
kipprates = ReactionRates.reaction_list[:kipp_rates]

##

logts = LinRange(3.0, 9.0, 30)
eosses = [EOS.EOSResults{Float64}() for i in 1:length(logts)]
for i in eachindex(logts)
    eosses[i].ρ = 0.0003053927091989808
    eosses[i].T = 10^logts[i]
end


xa = [0.8, 0.2, 0.0, 0.0, 0.0, 0.0]
xa_index = Dict(:H1 => 1, :D2 => 2, :He4 => 3, :C12 => 4, :O16 => 5, :N14 => 6)

##

# choose sets of reactions that you want to plot

reactions =       [[jinarates[:H1_H1_to_D2_betplus_w_x_0],
                    jinarates[:H1_H1_to_D2_xxec_w_x_0]],

                    [kipprates[:kipp_pp]],
              
                   [toyrates[:toy_pp]]
                   
                   ]


##


rates_list = []
for set in reactions

    rates = zeros(length(logts));
    for i in eachindex(logts)
        for reaction in set
            rates[i] += ReactionRates.get_reaction_rate(reaction, eosses[i], xa, xa_index)
        end

    end

    log_rates = log10.(rates);
    push!(rates_list, log_rates)

end

##


##
f = Figure();
ax = Axis(f[1, 1], xlabel=L"\log T/\textrm{K}", ylabel=L"\log R / {\textrm{g}}^{-1} {\textrm{s}}^{-1}");
lines!(ax, logts, rates_list[1], label="Jina PP")
lines!(ax, logts, rates_list[2], label="Kipp PP")
lines!(ax, logts, rates_list[3], label="Toy PP")
# lines!(ax, logts, k_rates, label="Kipp")
axislegend(position=:lt)
f

