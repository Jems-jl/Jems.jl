export AbundanceList, get_mass_fractions


##
using Jems.Chem #doe terug weg
Chem.isotope_list[:Au197]

##
@kwdef struct AbundanceList
    massfractions::Dict{Symbol, Float64}
    abundance_sources::Dict{Symbol, String}
    species_names::Vector{Symbol}
end

function AbundanceList(path) #returns object of type AbundanceList
    #count the number of elements
    Nb_elements = 0
    X = 0.0
    Y = 0.0
    Z = 0.0
    open(path) do io
        readline(io); readline(io)
        X = parse(Float64,split(readline(io), "X =")[2]);
        Y = parse(Float64,split(readline(io), "Y =")[2]);
        Z = parse(Float64,split(readline(io), "Z =")[2]);
        readline(io); readline(io)
        lines = readlines(io)
        Nb_elements = length(lines)
    end
    massfractions = Dict{Symbol, Float64}()
    species_names = Vector{Symbol}()
    abundance_sources = Dict{Symbol, String}()
    open(path) do io
        readline(io); readline(io); readline(io); 
        readline(io); readline(io); readline(io); readline(io)
        for i = 1:Nb_elements
            information = split(readline(io), r"\s+")
            Z = information[1]
            chemsymbol = information[2]
            abundance = parse(Float64, information[3])
            source = information[4]
            isotope_symbol = Symbol(chemsymbol)
            isotope = isotope_list[isotope_symbol]
            X_i = (isotope.mass / isotope_list[:H1].mass) * X * 10^(abundance - 12)
            massfractions[isotope_symbol] = X_i
            abundance_sources[isotope_symbol] = source
            push!(species_names, isotope_symbol)
            @show isotope_symbol, isotope, X_i, abundance
        end
    end
    return AbundanceList(massfractions, abundance_sources, species_names)
end
#prepare a dictionary that contains mixtures
abundance_lists = Dict{Symbol, AbundanceList}()
#add the Asplund 2009 mixture
abundance_lists[:ASG_09] = AbundanceList(pkgdir(Chem, "data/ChemData", "asplund2009.data"))
#check the species included
abundance_lists[:ASG_09].species_names
#get the mass fraction of a specific isotope of the :ASG_09 mixture
abundance_lists[:ASG_09].massfractions[:Li7]
#get the source of the abundance of a specific isotope of the :ASG_09 mixture
abundance_lists[:ASG_09].abundance_sources[:Li7]


##
function get_mass_fractions(abundance_list::AbundanceList, network, X, Z, Dfraction)
    #count the sum of all metals
    sum_of_metals = 0.0
    for i in eachindex(network.species_names)
        species = network.species_names[i]
        if species ≠ :H1 && species ≠ :He4
            if species in keys(abundance_list.massfractions)
                sum_of_metals += abundance_list.massfractions[species]
            end
        end
    end
    #calculate the mass fractions
    massfractions = Dict()
    for i in eachindex(network.species_names)
        species = network.species_names[i]
        if species ≠ :H1 && species ≠ :He4
            if species in abundance_list.species_names
                fraction = Z / sum_of_metals
                massfractions[species] = abundance_list.massfractions[species] * fraction
            else
                massfractions[species] = 0.0 #put to zero if species not in abundance list
            end
        end
    end
    #now put hydrogen and helium mass fractions
    massfractions[:H1] = X * (1-Dfraction)
    massfractions[:D2] = X * Dfraction #choose some fraction of hydrogen to be deuterium
    massfractions[:He4] = 1-X-Z
    return massfractions
end


##

###
using Jems
using Jems.NuclearNetworks
#creating a NuclearNetwork, also including some exotic species
net = NuclearNetwork([:H1,:He4,:C12, :N14, :O16, :Al27, :Ba138, :U238], [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
massfractions_for_network = get_mass_fractions(abundance_lists[:ASG_09], net, 0.7, 0.05,0.03)
@show massfractions_for_network