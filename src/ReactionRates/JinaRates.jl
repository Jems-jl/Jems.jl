file_contents = open(pkgdir(Chem, "data/ReactionRatesData", "Jina_reactionrates.data")) do io
    read(io, String)
end

"""
    JinaReactionRate{TT<:Real}<:ReactionRates.AbstractReactionRate

Struct that holds the following information for a given reaction rate:
    name: name of the reaction as a symbol
    iso_in: vector that contains all elements on the LHS of the reaction
    iso_out: vector that contains all elements on the RHS of the reaction
    Qvalue: Q-value of the reaction (MeV)
    coeff: different a_i values of the reaction. Contains a vector of 7 values
    set_label: Symbol containing set label of the reaction
    res_rate: A 1 character flag symbol:
        when blank or n it is a non-resonant rate
        when r it is a resonant rate
        when w it is a weak rate.
    rev_rate: a 1 character flag symbol which is set to 'v' when it is a reverse rate.
    chapter: chapter this reaction is in

"""

struct JinaReactionRate{TT<:Real}<:ReactionRates.AbstractReactionRate
    name::Symbol
    iso_in::Vector{Symbol}
    iso_out::Vector{Symbol}
    Qvalue::TT
    coeff::Vector{TT}
    set_label::Symbol
    res_rate::Symbol
    rev_rate::Symbol
    chapter::Int64
end


"""

    add_to_references(main_dict, ref_dict, reaction, new_info::JinaReactionRate)

Function to identify rates with the same reaction equation
Evaluates if a reaction rate is already in the reference dictionary ref_dict

If the reaction rate does not exist allready in the reference dictionary:
    added as a new key to the reference dictionary
    the value of the key is a list containing all variations of the specific reaction
    the reaction will be added to the main dictionary 

If the reaction rate allready exists in the reference dictionary:
    keys in the main dictionary update so they have unique keys
    value of the key of the reaction in ref_dict is updated so all the unique versions of the rate are in

"""
function add_to_references(main_dict, ref_dict, reaction, new_info::JinaReactionRate)

    # main_dict = general dictionary containing all JINA Reaction rates
    # ref_dict  = dictionary containing all unique versions of each reaction rates
    # reaction  = Symbol of the reaction that has to be added to the main dictionary
    # new_info  = JinaReactionRate of the new rate

    if haskey(ref_dict, reaction)

        if ref_dict[reaction] == []

            cur_info = main_dict[reaction]

            cur_set_label = cur_info.set_label
            cur_res_rate  = cur_info.res_rate
            cur_rev_rate = cur_info.rev_rate

            new_set_label = new_info.set_label
            new_res_rate  = new_info.res_rate
            new_rev_rate = new_info.rev_rate

            reaction_string_cur = "$(reaction)_$(cur_set_label)_$(cur_res_rate)_$(cur_rev_rate)"
            reaction_string_new = "$(reaction)_$(new_set_label)_$(new_res_rate)_$(new_rev_rate)"

            reaction_symbol_cur = Symbol(replace(reaction_string_cur, ' ' => 'x'))
            reaction_symbol_new = Symbol(replace(reaction_string_new, ' ' => 'x'))

            list = []
            push!(list, reaction_symbol_cur)
            push!(list, reaction_symbol_new)
            ref_dict[reaction] = list

            main_dict[reaction_symbol_cur] = cur_info
            main_dict[reaction_symbol_new] = new_info

        else

            new_set_label = new_info.set_label
            new_res_rate  = new_info.res_rate
            new_rev_rate = new_info.rev_rate

            reaction_string_new = "$(reaction)_$(new_set_label)_$(new_res_rate)_$(new_rev_rate)"
            reaction_symbol_new = Symbol(replace(reaction_string_new, ' ' => 'x'))

            list = ref_dict[reaction]
            push!(list, reaction_symbol_new)
            ref_dict[reaction] = list


            main_dict[reaction_symbol_new] = new_info

        end

    else 
        
        ref_dict[reaction]  = []        
        main_dict[reaction] = new_info 

    end
    
end

"""
    correct_names(JINA_name)

This function will return the name that corresponds with the JEMS isotope database

JINA_name is the name of the element as it is given in the JINA library (without the extra spaces) as a string
RETURN_name is the corrected name given as a string

"""
function correct_names(JINA_name)
    change_name = Dict("p" => "H1", "d" => "D2", "t" => "T3", "n" => "n")

    if haskey(change_name, JINA_name)
        RETURN_name = change_name[JINA_name]
    else
        RETURN_name = uppercase(JINA_name[1]) * lowercase(JINA_name[2:end])
    end

    return RETURN_name
end

"""

    read_set(dataset, dictionary, reference_dictionary)

    * explanation *

"""
function read_set(dataset, dictionary, reference_dictionary)

    chap = 0
    n = 0

    while n <= lastindex(dataset) - 225       
    
        if dataset[(n + 1)] == ' '

            set_label = Symbol(dataset[(n + 44):(n + 47)])
            res_rate  = Symbol(dataset[(n + 48)])
            rev_rate = Symbol(dataset[(n + 49)])

            a0 = parse(Float64, dataset[(n + 76) : (n + 88) ])
            a1 = parse(Float64, dataset[(n + 89) : (n + 101)])
            a2 = parse(Float64, dataset[(n + 102): (n + 114)])
            a3 = parse(Float64, dataset[(n + 115): (n + 127)])
    
            a4 = parse(Float64, dataset[(n + 151): (n + 163)])
            a5 = parse(Float64, dataset[(n + 164): (n + 176)])
            a6 = parse(Float64, dataset[(n + 177): (n + 189)])

            a  = [a0, a1, a2, a3, a4, a5, a6]

            if chap == 1  # A -> B

                char_1_JINA = strip(dataset[(n + 6) : (n + 10)]); char_1 = correct_names(char_1_JINA)
                char_2_JINA = strip(dataset[(n + 11): (n + 15)]); char_2 = correct_names(char_2_JINA)

                reaction_symbol = Symbol(char_1 * "_to_" * char_2)

                elem_1 = [Symbol(char_1)];
                elem_2 = [Symbol(char_2)];

            elseif chap == 2  # A -> B + C

                char_1_JINA = strip(dataset[(n + 6) : (n + 10)]); char_1 = correct_names(char_1_JINA)
                char_2_JINA = strip(dataset[(n + 11): (n + 15)]); char_2 = correct_names(char_2_JINA)
                char_3_JINA = strip(dataset[(n + 16): (n + 20)]); char_3 = correct_names(char_3_JINA)
            
                reaction_symbol = Symbol(char_1 * "_to_" * char_2 * "_" * char_3)
                
                elem_1 = [Symbol(char_1)];
                elem_2 = [Symbol(char_2), Symbol(char_3)];

            elseif chap == 3  # A -> B + C + D

                char_1_JINA = strip(dataset[(n + 6) : (n + 10)]); char_1 = correct_names(char_1_JINA)
                char_2_JINA = strip(dataset[(n + 11): (n + 15)]); char_2 = correct_names(char_2_JINA)
                char_3_JINA = strip(dataset[(n + 16): (n + 20)]); char_3 = correct_names(char_3_JINA)
                char_4_JINA = strip(dataset[(n + 21): (n + 25)]); char_4 = correct_names(char_4_JINA)
            
                reaction_symbol = Symbol(char_1 * "_to_" * char_2 * "_" * char_3 * "_" *char_4)
                
                elem_1 = [Symbol(char_1)];
                elem_2 = [Symbol(char_2), Symbol(char_3), Symbol(char_4)];

            elseif chap == 4  # A + B -> C

                char_1_JINA = strip(dataset[(n + 6) : (n + 10)]); char_1 = correct_names(char_1_JINA)
                char_2_JINA = strip(dataset[(n + 11): (n + 15)]); char_2 = correct_names(char_2_JINA)
                char_3_JINA = strip(dataset[(n + 16): (n + 20)]); char_3 = correct_names(char_3_JINA)
            
                reaction_symbol = Symbol(char_1 * "_" * char_2 * "_to_" * char_3)
                
                elem_1 = [Symbol(char_1), Symbol(char_2)];
                elem_2 = [Symbol(char_3)];

            elseif chap == 5  # A + B -> C + D

                char_1_JINA = strip(dataset[(n + 6) : (n + 10)]); char_1 = correct_names(char_1_JINA)
                char_2_JINA = strip(dataset[(n + 11): (n + 15)]); char_2 = correct_names(char_2_JINA)
                char_3_JINA = strip(dataset[(n + 16): (n + 20)]); char_3 = correct_names(char_3_JINA)
                char_4_JINA = strip(dataset[(n + 21): (n + 25)]); char_4 = correct_names(char_4_JINA)
            
                reaction_symbol = Symbol(char_1 * "_" * char_2 * "_to_" * char_3 * "_" *char_4)
                
                elem_1 = [Symbol(char_1), Symbol(char_2)];
                elem_2 = [Symbol(char_3), Symbol(char_4)];

            elseif chap == 6  # A + B -> C + D + E

                char_1_JINA = strip(dataset[(n + 6) : (n + 10)]); char_1 = correct_names(char_1_JINA)
                char_2_JINA = strip(dataset[(n + 11): (n + 15)]); char_2 = correct_names(char_2_JINA)
                char_3_JINA = strip(dataset[(n + 16): (n + 20)]); char_3 = correct_names(char_3_JINA)
                char_4_JINA = strip(dataset[(n + 21): (n + 25)]); char_4 = correct_names(char_4_JINA)
                char_5_JINA = strip(dataset[(n + 26): (n + 30)]); char_5 = correct_names(char_5_JINA)

                reaction_symbol = Symbol(char_1 * "_" * char_2 * "_to_" * char_3 * "_" * char_4 * "_" * char_5)
                
                elem_1 = [Symbol(char_1), Symbol(char_2)];
                elem_2 = [Symbol(char_3), Symbol(char_4), Symbol(char_5)];

            elseif chap == 7  # A + B -> C + D + E + F

                char_1_JINA = strip(dataset[(n + 6) : (n + 10)]); char_1 = correct_names(char_1_JINA)
                char_2_JINA = strip(dataset[(n + 11): (n + 15)]); char_2 = correct_names(char_2_JINA)
                char_3_JINA = strip(dataset[(n + 16): (n + 20)]); char_3 = correct_names(char_3_JINA)
                char_4_JINA = strip(dataset[(n + 21): (n + 25)]); char_4 = correct_names(char_4_JINA)
                char_5_JINA = strip(dataset[(n + 26): (n + 30)]); char_5 = correct_names(char_5_JINA)
                char_6_JINA = strip(dataset[(n + 31): (n + 35)]); char_6 = correct_names(char_6_JINA)

                reaction_symbol = Symbol(char_1 * "_" * char_2 * "_to_" *
                                         char_3 * "_" * char_4 * "_" * char_5 * "_" * char_6)
                
                elem_1 = [Symbol(char_1), Symbol(char_2)];
                elem_2 = [Symbol(char_3), Symbol(char_4), Symbol(char_5), Symbol(char_6)];

            elseif chap == 8  # A + B + C -> D

                char_1_JINA = strip(dataset[(n + 6) : (n + 10)]); char_1 = correct_names(char_1_JINA)
                char_2_JINA = strip(dataset[(n + 11): (n + 15)]); char_2 = correct_names(char_2_JINA)
                char_3_JINA = strip(dataset[(n + 16): (n + 20)]); char_3 = correct_names(char_3_JINA)
                char_4_JINA = strip(dataset[(n + 21): (n + 25)]); char_4 = correct_names(char_4_JINA)

                reaction_symbol = Symbol(char_1 * "_" * char_2 * "_" * char_3 * "_to_" * char_4)
                
                elem_1 = [Symbol(char_1), Symbol(char_2), Symbol(char_3)];
                elem_2 = [Symbol(char_4)];

             elseif chap == 9  # A + B + C -> D + E

                char_1_JINA = strip(dataset[(n + 6) : (n + 10)]); char_1 = correct_names(char_1_JINA)
                char_2_JINA = strip(dataset[(n + 11): (n + 15)]); char_2 = correct_names(char_2_JINA)
                char_3_JINA = strip(dataset[(n + 16): (n + 20)]); char_3 = correct_names(char_3_JINA)
                char_4_JINA = strip(dataset[(n + 21): (n + 25)]); char_4 = correct_names(char_4_JINA)
                char_5_JINA = strip(dataset[(n + 26): (n + 30)]); char_5 = correct_names(char_5_JINA)

                reaction_symbol = Symbol(char_1 * "_" * char_2 * "_" * char_3 * "_to_" * char_4 * "_" * char_5)
                
                elem_1 = [Symbol(char_1), Symbol(char_2), Symbol(char_3)];
                elem_2 = [Symbol(char_4), Symbol(char_5)];

            end

            Q_value = parse(Float64, dataset[(n + 53):(n + 64)])

            reaction_info = JinaReactionRate(reaction_symbol, elem_1, elem_2, Q_value, a, set_label, res_rate, rev_rate,
                                             chap)
            add_to_references(dictionary, reference_dictionary, reaction_symbol, reaction_info)
        else

            chap += 1

        end

        n += 225

    end
        
end


"""

    get_reaction_rate(reaction::JinaReactionRate, eos00::EOSResults{TT}, xa::AbstractVector{TT}, xa_index::Dict{Symbol,Int})

    * explanation *

"""
function get_reaction_rate(reaction::JinaReactionRate, eos00::EOSResults{TT}, xa::AbstractVector{TT}, xa_index::Dict{Symbol,Int})::TT where{TT}

    # determine λ

    T_9 = (eos00.T / 1e9)
    a = reaction.coeff

    x = a[1] +  a[7] * log(T_9)

    for i in 2:6
        x += a[i] * T_9^((2(i-1) - 5)/3)
    end

    λ = exp(x)

    println(λ)

    # determine elements and how many times they occur
    # code gives a dictionary with the elements and how many times they occur

    LHS_elements = reaction.iso_in

    sorted_elements = Dict{Symbol, Int}()
    
    for element in LHS_elements
        if haskey(sorted_elements, element)
            sorted_elements[element] += 1
        else
            sorted_elements[element] = 1
        end
    end

    elements   = collect(keys(sorted_elements))      # elements that occur on the LHS
    N_elements = collect(values(sorted_elements))    # how many times they occur on the LHS

    println(elements, N_elements)
    # determine all needed parameters for every element

    ν = -1
    factors = []

    for index in eachindex(elements)

        elem = elements[index]                                       # A
        N_elem  = N_elements[index]                                  # N_A
        X_elem  = xa[xa_index[elem]]                                 # X_A
        m_elem  = Chem.isotope_list[elem].mass * Constants.AMU       # m_A

        Y_elem = X_elem / (m_elem * Constants.AVO)                   # Y_A
        ν += N_elem
        factor_elem = Y_elem^(N_elem) / factorial(N_elem)

        push!(factors, factor_elem)

    end

    # Calculate the reaction rate
    println(ν)
    ρ = eos00.ρ
    RR = ρ^ν * λ

    for factor in factors
        RR *= factor
    end

    return RR * Constants.AVO

end

















