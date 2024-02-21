export AbundanceList, get_mass_fractions

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
        end
    end
    return AbundanceList(massfractions, abundance_sources, species_names)
end

#################################################################################################
#prepare a dictionary that contains mixtures
abundance_lists = Dict{Symbol, AbundanceList}()
#add the Asplund 2009 mixture
abundance_lists[:ASG_09] = AbundanceList(pkgdir(Chem, "data/ChemData", "asplund2009.data"))

#Some functionalities:
#Check the species included
#    abundance_lists[:ASG_09].species_names
#get the mass fraction of a specific isotope of the :ASG_09 mixture
#    abundance_lists[:ASG_09].massfractions[:Li7]
#get the source of the abundance of a specific isotope of the :ASG_09 mixture
#    abundance_lists[:ASG_09].abundance_sources[:Li7]
#
##################################################################################################

"""
    get_mass_fractions(abundance_list::AbundanceList, network, X, Z, Dfraction)

Given an AbundanceList (= a certain mixture), a nuclear network, fractions X and Z and the deuterium fraction Dfraction, this function
rescales the mass fraction of the mixture to the disired X and Z values, keeping the relative abundances of the metals fixed. All species
not in the network are given a mass fraction of zero. The function returns a dictionary with the mass fractions of the species in the network.
"""

function get_mass_fractions(abundance_list::AbundanceList, species_names, X, Z, Dfraction)
    #count the sum of all metals
    sum_of_metals = 0.0
    for i in eachindex(species_names)
        species = species_names[i]
        if species ≠ :H1 && species ≠ :He4
            if species in abundance_list.species_names
                sum_of_metals += abundance_list.massfractions[species]
            end
        end
    end
    #calculate the mass fractions
    massfractions = Dict()
    for i in eachindex(species_names)
        species = species_names[i]
        if species ≠ :H1 && species ≠ :He4
            if species in abundance_list.species_names
                fraction = Z / sum_of_metals
                massfractions[species] = abundance_list.massfractions[species] * fraction
            elseif species ≠ :D2
                println("Species $(species) is not in the abundance list, setting its mass fraction to 0.0")
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