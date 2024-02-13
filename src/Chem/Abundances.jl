using Jems.Chem
using Jems.NuclearNetworks
println("RESTART ##################################")
Chem.isotope_list[:Zn64]

##
@kwdef struct AbundanceList
    massfractions::Dict{Symbol, Float64}
    species_names::Vector{Symbol}
end

function get_abundance_list(path)
    #count the number of elements
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
    open(path) do io
        readline(io); readline(io); readline(io); 
        readline(io); readline(io); readline(io); readline(io)
        for i = 1:30
            information = split(readline(io), r"\s+")
            Z = information[1]
            chemsymbol = information[2]
            abun_phot = parse(Float64, information[3])
            abun_met = information[4]
            isotope_symbol = Symbol(chemsymbol)
            isotope = Chem.isotope_list[isotope_symbol]
            X_i = (isotope.mass / Chem.isotope_list[:H1].mass) * X * 10^(abun_phot - 12)
            massfractions[isotope_symbol] = X_i
            push!(species_names, isotope_symbol)
        end
    end
    return AbundanceList(massfractions, species_names)
end
#prepare dictionary
abundance_lists = Dict{Symbol, AbundanceList}()
#get the abundances from Asplund2009
abundance_lists[:ASG_09] = get_abundance_list("src/Chem/asplund2009.txt")
abundance_lists[:ASG_09]
#check the species included
abundance_lists[:ASG_09].species_names
#check the mass fraction of a specific isotope
abundance_lists[:ASG_09].massfractions[:Li7]


##
function get_mass_fractions(abundance_list::AbundanceList, network::NuclearNetwork, X, Z,Dfraction)
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
    @show sum_of_metals
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
    massfractions[:H2] = X * Dfraction #choose some fraction of hydrogen to be deuterium
    massfractions[:He4] = 1-X-Z
    return massfractions
end

##
println("Restart ######################################################################################################")
net = NuclearNetwork([:H1,:He4,:C12, :N14, :O16], [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
a = get_mass_fractions(abundance_lists[:ASG_09], net, 0.7, 0.05,0.03)
@show a
#println(a[:O16] / a[:C12])
#abundance_lists[:ASG_09].massfractions[:O16] / abundance_lists[:ASG_09].massfractions[:C12]

##








##
################################################################"
println("RESTART ##################################")
open(path) do io
    readline(io); readline(io)
    X = parse(Float64,split(readline(io), "X =")[2]);
    Y = parse(Float64,split(readline(io), "Y =")[2]);
    Z = parse(Float64,split(readline(io), "Z =")[2]);
    println(X,Y,Z)
    readline(io); readline(io)
    lines = readlines(io)
    Nb_elements = length(lines)
    println("Number of elements = $Nb_elements")
end
#make an empty array
massfractions = zeros(Nb_elements)
open(path) do io
    println("begonnen aan loop")
    readline(io); readline(io); readline(io); 
    readline(io); readline(io); readline(io); readline(io)
    for i = 1:30
        information = split(readline(io), r"\s+")
        Z = information[1]; chemsymbol = information[2]; abun_phot = information[3]; abun_met = information[4]
        isotope_symbol = Symbol(chemsymbol)
        isotope = Chem.isotope_list[isotope_symbol]
        X_i = (isotope.mass / Chem.isotope_list[:H1].mass) * X * 10^(abun_phot - 12)
        massfractions[i] = X_i
        println(information)
        println("Z = $Z")
        println("symbol = $symbol")
        println("abun_phot = $abun_phot")
        println("abun_met = $abun_met")
        println("Xi = $Xi")
        break
    end
end

##