using BenchmarkTools
using Jems.Chem
using Jems.Constants
using Jems.EOS
using Jems.Opacity
using Jems.NuclearNetworks
using Jems.Turbulence
using Jems.StellarModels
using Jems.Evolution
using Jems.ReactionRates


##

References = ReactionRates.jina_references
Jina_Rates = ReactionRates.jina_rates
# read_set(file_contents, Jina_Rates, References)

##

r = EOSResults{Float64}();
eos = EOS.IdealEOS(false);     
logT = LinRange(8,9,100);      
rates_Kipp_O16alpha = zeros(100);            
rates_Jina_O16alpha = zeros(100); 
rates_Kipp_3alphaCF88 = zeros(100)
rates_Kipp_3alphaA99  = zeros(100)
rates_Jina_3alpha = zeros(100)
rates_Kipp_C12alpha = zeros(100)
rates_Jina_C12alpha = zeros(100)


##

for i in eachindex(logT)                    

    # kipp_C12alpha

    # set_EOS_resultsTρ!(eos, r, logT[i]*log(10), log(1), [0.2, 0.8, 0.0], [:C12, :He4, :O16])              
    # rates_Kipp_C12alpha[i] = ReactionRates.get_reaction_rate(ReactionRates.reaction_list[:kipp_rates][:kipp_C12alpha], r, [0.2, 0.8, 0.0], Dict{Symbol, Int64}(:He4 => 1, :C12 => 2, :O16 => 3))
    # rates_Jina_C12alpha[i] = ReactionRateS.get_reaction_rate(Jina_Rates[:He4_C12_to_O16], r, [0.2, 0.8, 0.0], Dict{Symbol, Int64}(:He4 => 1, :C12 => 2, :O16 => 3))                            

    # 016 alpha

    # set_EOS_resultsTρ!(eos, r, logT[i]*log(10), log(1), [0.2, 0.8, 0.0], [:O16, :He4, :Ne20])              
    # rates_Kipp_O16alpha[i] = ReactionRates.get_reaction_rate(ReactionRates.reaction_list[:kipp_rates][:kipp_O16alpha], r, [0.2, 0.8, 0.0], Dict{Symbol, Int64}(:O16 => 1, :He4 => 2, :Ne20 => 3))
    # rates_Jina_O16alpha[i] = ReactionRates.get_reaction_rate(Jina_Rates[:He4_O16_to_Ne20_co10_n_x], r, [0.2, 0.8, 0.0], Dict{Symbol, Int64}(:O16 => 1, :He4 => 2, :Ne20 => 3))                            
    # rates_Jina_O16alpha[i] += ReactionRates.get_reaction_rate(Jina_Rates[:He4_O16_to_Ne20_co10_r_x], r, [0.2, 0.8, 0.0], Dict{Symbol, Int64}(:O16 => 1, :He4 => 2, :Ne20 => 3))                            
    # rates_Jina_O16alpha[i] += ReactionRates.get_reaction_rate(Jina_Rates[:He4_O16_to_Ne20], r, [0.2, 0.8, 0.0], Dict{Symbol, Int64}(:O16 => 1, :He4 => 2, :Ne20 => 3))                            

    # H + H

    # rates_Jina_HH[i] = ReactionRates.get_reaction_rate(Jina_Rates[:H2_H2_to_], r, [0.2, 0.8, 0.0], Dict{Symbol, Int64}(:O16 => 1, :He4 => 2, :Ne20 => 3))                            
    

    # Triple alpha

    set_EOS_resultsTρ!(eos, r, logT[i]*log(10), log(1), [0.2, 0.8, 0.0], [:C12, :He4])              
    rates_Kipp_3alphaCF88[i] = ReactionRates.get_reaction_rate(ReactionRates.reaction_list[:kipp_rates][:kipp_3alphaCF88], r, [0.2, 0.8], Dict{Symbol, Int64}(:He4 => 1, :C12 => 2))
    rates_Kipp_3alphaA99[i]  = ReactionRates.get_reaction_rate(ReactionRates.reaction_list[:kipp_rates][:kipp_3alphaA99],  r, [0.2, 0.8], Dict{Symbol, Int64}(:He4 => 1, :C12 => 2))
    
    rates_Jina_3alpha[i] =  ReactionRates.get_reaction_rate(Jina_Rates[:He4_He4_He4_to_C12_fy05_r_x_0], r, [0.2, 0.8], Dict{Symbol, Int64}(:He4 => 1, :C12 => 2))                            
    rates_Jina_3alpha[i] += ReactionRates.get_reaction_rate(Jina_Rates[:He4_He4_He4_to_C12_fy05_r_x_1], r, [0.2, 0.8], Dict{Symbol, Int64}(:He4 => 1, :C12 => 2))                            
    rates_Jina_3alpha[i] += ReactionRates.get_reaction_rate(Jina_Rates[:He4_He4_He4_to_C12_fy05_n_x_0], r, [0.2, 0.8], Dict{Symbol, Int64}(:He4 => 1, :C12 => 2))                            


end


##

# rates_Kipp_C12alpha .= rates_Kipp_C12alpha .* (Constants.MEV_TO_ERGS)^(-1)
# rates_Jina_C12alpha .= rates_Jina_C12alpha .* (Constants.AVO)

# rates_Kipp_O16alpha .= rates_Kipp_O16alpha .* (Constants.MEV_TO_ERGS)^(-1)
# rates_Jina_O16alpha .= rates_Jina_O16alpha .* (Constants.AVO)

rates_Kipp_3alphaCF88 .= rates_Kipp_3alphaCF88 .* (Constants.MEV_TO_ERGS)^(-1)
rates_Kipp_3alphaA99 .= rates_Kipp_3alphaA99 .* (Constants.MEV_TO_ERGS)^(-1)
rates_Jina_3alpha .= rates_Jina_3alpha .* (Constants.AVO)



##



using CairoMakie

f = Figure();
ax = Axis(f[1, 1]; xlabel=L"\log_{10}(T_\mathrm{eff}/[K])", ylabel="Reaction rate Kipp [mol / g * s]", xreversed=false)
lines!(ax, logT, log10.(rates_Kipp_3alphaCF88), label = "rates_Kipp_3alphaCF88")
lines!(ax, logT, log10.(rates_Kipp_3alphaA99), label = "rates_Kipp_3alphaA99")
lines!(ax, logT, log10.(rates_Jina_3alpha), label = "rates_Jina_3alpha")

# lines!(ax, logT, log10.(rates_Kipp_O16alpha), label = "rates_Kipp_O16alpha")
# lines!(ax, logT, log10.(rates_Jina_O16alpha), label = "rates_Jina_O16alpha")

# lines!(ax, logT, log10.(rates_Kipp_C12alpha), label = "rates_Kipp_C12log10")
# lines!(ax, logT, log10.(rates_Jina_C12alpha) .+ 0.9, label = "rates_Jina_C12alpha")

axislegend(ax; position=:rb)
f

##




