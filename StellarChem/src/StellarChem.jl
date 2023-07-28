module StellarChem

export Isotope

struct Isotope
    Z::Int64
    A::Int64
    name::String
    mass::Float64 #in atomic mass units
end

function get_isotope_list()
    Niso = 0
    # First pass to count isotopes
    open(pkgdir(StellarChem,"data","isotope.data")) do io
        readline(io) # throw out the first line
        lines = readlines(io)
        Niso = (length(lines)+1)/8 #last isotope is missing an empty line
    end
    isotope_list = Dict{Symbol,Isotope}()
    # second pass to read isotopes
    open(pkgdir(StellarChem,"data","isotope.data")) do io
        readline(io) # throw out the first line
        for i in 1:Niso
            Z = parse(Int64, split(readline(io)," = ")[2]) # Read atomic number
            name = split(readline(io)," = ")[2] # Read name
            A = parse(Int64, split(readline(io)," = ")[2]) # Read mass number
            name_symbol = Symbol(name,A) # create symbol for dictionary
            mass = parse(Float64, split(split(readline(io)," = ")[2],"(")[1]) # Read atomic number
            isotope_list[name_symbol] = Isotope(Z,A,name,mass)
            #skip remaining lines
            readline(io)
            readline(io)
            readline(io)
            readline(io)
        end
    end
    return isotope_list
end

end # module StellarChem
