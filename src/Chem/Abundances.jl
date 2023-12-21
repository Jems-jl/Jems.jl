@kwdef struct AbundanceList
    # elements::Vector{}  [:H1, :He4]
    #also store the numbers
    #indicate which of the two
    #also store the mass fractions
    massfractions::Dict{Symbol, Float64}
end


function get_abundance_list(path)
    #count the number of elements
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
    massfractions = Dict{Symbol, Float64}()
    open(path) do io
        readline(io); readline(io); readline(io); 
        readline(io); readline(io); readline(io); readline(io)
        for i = 1:30
            information = split(readline(io), r"\s+")
            Z = information[1]; chemsymbol = information[2]; abun_phot = information[3]; abun_met = information[4]
            isotope_symbol = Symbol(chemsymbol)
            isotope = Chem.isotope_list[isotope_symbol]
            X_i = (isotope.mass / Chem.isotope_list[:H1].mass) * X * 10^(abun_phot - 12)
            massfractions[isotope_symbol] = X_i
            break
        end
    end
    return AbundanceList(massfractions)
end
#make dictionary
abundance_lists = Dict{Symbol, AbundanceList}()
abundance_list[:ASG_09] = get_abundance_list("src/Chem/asplund2009.txt")

function get_mass_fractions(abundance_list::AbundanceList, network::NuclearNetwork, X, Z)
    #count the sum of all metals
    sum_of_metals = 0.0
    for i in eachindex(network.species_names)
        species = network.species_names[i]
        if species ≠ :H1 && species ≠ :He4
            if species in abundance_list.massfractions
                sum_of_metals = sum_of_metals + abundance_list.massfractions[species]
            end
        end
    end
    #calculate the mass fractions
    for i in eachindex(network.species_names)
        species = network.species_names[i]
        if species == :H1
            massfractions[species] = X
        elseif species == :He4
            massfractions[species] = 1-X-Z
        else
            if species in abundance_list.massfractions
                fraction = Z / sum_of_metals
                massfractions[species] = abundance_list.massfractions[species] * fraction
            else
                massfractions[species] = 0.0 #put to zero if species not in abundance list
            end
        end
    end
    return massfractions
end













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