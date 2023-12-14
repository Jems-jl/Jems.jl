struct Abundancelist
    path::String
end

function get_abundance_list(abundancelist::Abundancelist,X,Y,Z)
    if X+Y+Z ≠ 1
        throw(ErrorException("X+Y+Z must be 1"))
    end
    Nb_elements = 0
    open(abundancelist.path) do io
        readline(io); readline(io)
        lines = readlines(io)
        Nb_elements = length(lines)
        println("Number of elements = $Nb_elements")
    end
    abundance_list = Dict{}()

    open(abundancelist.path) do io
        readline(io); readline(io)
        for i = 1:Nb_elements
            information = split(readline(io), "\t")
            #println(information)
            Zelement = information[1]
            symbol = information[2]
            abun_phot = parse(Float16,information[3])
            abun_met = parse(Float16,information[4])
            println("abun = $abun_phot")
            if symbol ≠ "H" && symbol ≠ "He"
                abundance_list[symbol] = abun_phot*Z;
            end
            if symbol == "H"
                abundance_list[symbol] = abun_phot*X
            end
            if symbol == "He"
                abundance_list[symbol] = abun_phot*Y
            end
        end
    end
    return abundance_list
end

##
path = "src/Chem/asplund2009.txt"
asplund2009 = Abundancelist(path)
test_dict = get_abundance_list(asplund2009,1/3,1/3,1/3)
test_dict["Ar"]

##

testDict = Dict{}()
testDict["test"]=5
testDict
