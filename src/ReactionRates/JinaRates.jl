"""
    JinaReactionRate{TT<:Real}<:ReactionRates.AbstractReactionRate

Struct that holds the following information for a given reaction rate:
    name: name of the reaction as a symbol
    iso_in: vector that contains all elements on the LHS of the reaction
    iso_out: vector that contains all elements on the RHS of the reaction
    Qvalue: Q-value of the reaction
    coeff: different a_i values of the reaction. Contains a vector of 7 values
    set_label: Symbol containing set label of the reaction
    res_rate: A 1 character flag symbol:
        when blank or n it is a non-resonant rate
        when r it is a resonant rate
        when w it is a weak rate.
    rev_rate: a 1 character flag symbol which is set to 'v' when it is a reverse rate.
    chapter: chapter this reaction is in

"""

struct JinaReactionRate{TT<:Real} <: ReactionRates.AbstractReactionRate
    name::Symbol
    iso_in::Vector{Symbol}
    num_iso_in::Vector{Int}
    iso_out::Vector{Symbol}
    num_iso_out::Vector{Int}
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
Evaluates if a reaction rate is already in the reference dictionary `ref_dict`

If the reaction rate does not exist already in the reference dictionary:
added as a new key to the reference dictionary
the value of the key is a list containing all variations of the specific reaction
the reaction will be added to the main dictionary

If the reaction rate already exists in the reference dictionary:
keys in the main dictionary update so they have unique keys
value of the key of the reaction in ref_dict is updated so all the unique versions of the rate are in
"""
function add_to_references(main_dict, ref_dict, reaction, new_info::JinaReactionRate)

    # main_dict = general dictionary containing all JINA Reaction rates
    # ref_dict  = dictionary containing all unique versions of each reaction rates
    # reaction  = Symbol of the reaction that has to be added to the main dictionary
    # new_info  = JinaReactionRate of the new rate

    new_set_label = new_info.set_label
    new_res_rate = new_info.res_rate
    new_rev_rate = new_info.rev_rate

    reaction_string_new = "$(reaction)_$(new_set_label)_$(new_res_rate)_$(new_rev_rate)_0"
    reaction_string_new = replace(reaction_string_new, ' ' => 'x')
    reaction_symbol_new = Symbol(replace(reaction_string_new, '+' => "plus"))

    if haskey(ref_dict, reaction)

        # build new reaction name
        reaction_string_short = "$(reaction)_$(new_set_label)_$(new_res_rate)_$(new_rev_rate)"
        reaction_string_short = replace(reaction_string_short, ' ' => 'x')
        reaction_symbol_short = Symbol(replace(reaction_string_short, '+' => "plus"))

        # check if this name is in the sublist of ref_dict
        list = ref_dict[reaction]

        if reaction_symbol_new in list
            new_list = [Symbol(string(symbol)[1:(end - 2)]) for symbol in list]
            curr_amount = count(x -> x == reaction_symbol_short, new_list)  # count how many match the short name
            reaction_symbol_new = Symbol("$(reaction_symbol_short)_$(curr_amount)")  # reaction_symbol_new is appended with this number
        end

        push!(list, reaction_symbol_new)   # add to the ref_dist sublist
        main_dict[reaction_symbol_new] = new_info  # add to main dict with the long key name

    else  # reference dictionary doesn't have reaction yet --> first of its kind
        ref_dict[reaction] = [reaction_symbol_new]
        main_dict[reaction_symbol_new] = new_info
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
    function sort_reaction(elements)

Extracts the unique elements and the amount of times they appear in the list `elements`.
"""
function sort_reaction(elements)
    sorted_elements = Dict{Symbol,Int}()  # dict mapping how may times an element occurs on the LHS of reaction

    for element in elements  # do the counting
        if haskey(sorted_elements, element)
            sorted_elements[element] += 1
        else
            sorted_elements[element] = 1
        end
    end

    elements_return = collect(keys(sorted_elements))  # unique elements that occur on the LHS
    N_elements_return = collect(values(sorted_elements))  # will equal reaction.num_iso_in, how many times they occur
    return [elements_return, N_elements_return]
end

"""
    function read_dataset(dataset, dictionary, reference_dictionary)

Parses a large string object `dataset`, containing JINA reaction rate data, into a main dictionary and a reference
dictionary. The main dictionary has entries with unique keys for every reaction rate in the database (double entries
are distinguised with "_0", "_1" at the end of their key Symbol), while the reference dictionary has one entry per
reaction equation, but has as value a list of all JINA rates for that reaction equation.

_Example:_ if H1 + H1 -> D2 has two rates in the JINA database, main dictionary will contain:
```
    dictionary[:H1_H1_to_D2_0] = _rate_info_(...)
    dictionary[:H1_H1_to_D2_1] = _rate_info_(...)
```
while the reference dict will contain:
```
    reference_dictionary[:H1_H1_to_D2] = [:H1_H1_to_D2_0, :H1_H1_to_D2_1]
```
For rates with only one entry in the JINA database, only an "_0" entry is made in the main dictionary, and the reference
dictionary contains the list with only that one key.
"""
function read_dataset(dataset, dictionary, reference_dictionary)
    chap = 0
    n = 0

    while n <= lastindex(dataset) - 225
        if dataset[(n + 1)] == ' '
            reaction = true

            set_label = Symbol(dataset[(n + 44):(n + 47)])
            res_rate = Symbol(dataset[(n + 48)])
            rev_rate = Symbol(dataset[(n + 49)])

            a0 = parse(Float64, dataset[(n + 76):(n + 88)])
            a1 = parse(Float64, dataset[(n + 89):(n + 101)])
            a2 = parse(Float64, dataset[(n + 102):(n + 114)])
            a3 = parse(Float64, dataset[(n + 115):(n + 127)])

            a4 = parse(Float64, dataset[(n + 151):(n + 163)])
            a5 = parse(Float64, dataset[(n + 164):(n + 176)])
            a6 = parse(Float64, dataset[(n + 177):(n + 189)])

            a = [a0, a1, a2, a3, a4, a5, a6]

            if chap == 1
                char_1_JINA = strip(dataset[(n + 6):(n + 10)])
                char_1 = correct_names(char_1_JINA)
                char_2_JINA = strip(dataset[(n + 11):(n + 15)])
                char_2 = correct_names(char_2_JINA)

                reaction_symbol = Symbol(char_1 * "_to_" * char_2)

                elem_1_u = [Symbol(char_1)]
                elem_2_u = [Symbol(char_2)]

                num_elem_1 = sort_reaction(elem_1_u)[2]
                num_elem_2 = sort_reaction(elem_2_u)[2]

                elem_1 = sort_reaction(elem_1_u)[1]
                elem_2 = sort_reaction(elem_2_u)[1]

                # num_elem_1 = zeros(Int, length(elem_1))
                # num_elem_2 = zeros(Int, length(elem_2))

                Q_value = parse(Float64, dataset[(n + 53):(n + 64)]) * Constants.MEV_TO_ERGS

                reaction_info = JinaReactionRate(reaction_symbol, elem_1, num_elem_1, elem_2, num_elem_2, Q_value, a,
                                                 set_label, res_rate, rev_rate, chap)
                add_to_references(dictionary, reference_dictionary, reaction_symbol, reaction_info)

            elseif chap == 2
                char_1_JINA = strip(dataset[(n + 6):(n + 10)])
                char_1 = correct_names(char_1_JINA)
                char_2_JINA = strip(dataset[(n + 11):(n + 15)])
                char_2 = correct_names(char_2_JINA)
                char_3_JINA = strip(dataset[(n + 16):(n + 20)])
                char_3 = correct_names(char_3_JINA)

                reaction_symbol = Symbol(char_1 * "_to_" * char_2 * "_" * char_3)

                elem_1_u = [Symbol(char_1)]
                elem_2_u = [Symbol(char_2), Symbol(char_3)]

                num_elem_1 = sort_reaction(elem_1_u)[2]
                num_elem_2 = sort_reaction(elem_2_u)[2]

                elem_1 = sort_reaction(elem_1_u)[1]
                elem_2 = sort_reaction(elem_2_u)[1]

                Q_value = parse(Float64, dataset[(n + 53):(n + 64)]) * Constants.MEV_TO_ERGS

                reaction_info = JinaReactionRate(reaction_symbol, elem_1, num_elem_1, elem_2, num_elem_2, Q_value, a,
                                                 set_label, res_rate, rev_rate, chap)
                add_to_references(dictionary, reference_dictionary, reaction_symbol, reaction_info)

            elseif chap == 3
                char_1_JINA = strip(dataset[(n + 6):(n + 10)])
                char_1 = correct_names(char_1_JINA)
                char_2_JINA = strip(dataset[(n + 11):(n + 15)])
                char_2 = correct_names(char_2_JINA)
                char_3_JINA = strip(dataset[(n + 16):(n + 20)])
                char_3 = correct_names(char_3_JINA)
                char_4_JINA = strip(dataset[(n + 21):(n + 25)])
                char_4 = correct_names(char_4_JINA)

                reaction_symbol = Symbol(char_1 * "_to_" * char_2 * "_" * char_3 * "_" * char_4)

                elem_1_u = [Symbol(char_1)]
                elem_2_u = [Symbol(char_2), Symbol(char_3), Symbol(char_4)]

                num_elem_1 = sort_reaction(elem_1_u)[2]
                num_elem_2 = sort_reaction(elem_2_u)[2]

                elem_1 = sort_reaction(elem_1_u)[1]
                elem_2 = sort_reaction(elem_2_u)[1]

                Q_value = parse(Float64, dataset[(n + 53):(n + 64)]) * Constants.MEV_TO_ERGS

                reaction_info = JinaReactionRate(reaction_symbol, elem_1, num_elem_1, elem_2, num_elem_2, Q_value, a,
                                                 set_label, res_rate, rev_rate, chap)
                add_to_references(dictionary, reference_dictionary, reaction_symbol, reaction_info)

            elseif chap == 4
                char_1_JINA = strip(dataset[(n + 6):(n + 10)])
                char_1 = correct_names(char_1_JINA)
                char_2_JINA = strip(dataset[(n + 11):(n + 15)])
                char_2 = correct_names(char_2_JINA)
                char_3_JINA = strip(dataset[(n + 16):(n + 20)])
                char_3 = correct_names(char_3_JINA)

                reaction_symbol = Symbol(char_1 * "_" * char_2 * "_to_" * char_3)

                elem_1_u = [Symbol(char_1), Symbol(char_2)]
                elem_2_u = [Symbol(char_3)]

                num_elem_1 = sort_reaction(elem_1_u)[2]
                num_elem_2 = sort_reaction(elem_2_u)[2]

                elem_1 = sort_reaction(elem_1_u)[1]
                elem_2 = sort_reaction(elem_2_u)[1]

                Q_value = parse(Float64, dataset[(n + 53):(n + 64)]) * Constants.MEV_TO_ERGS

                reaction_info = JinaReactionRate(reaction_symbol, elem_1, num_elem_1, elem_2, num_elem_2, Q_value, a,
                                                 set_label, res_rate, rev_rate, chap)
                add_to_references(dictionary, reference_dictionary, reaction_symbol, reaction_info)

            elseif chap == 5
                char_1_JINA = strip(dataset[(n + 6):(n + 10)])
                char_1 = correct_names(char_1_JINA)
                char_2_JINA = strip(dataset[(n + 11):(n + 15)])
                char_2 = correct_names(char_2_JINA)
                char_3_JINA = strip(dataset[(n + 16):(n + 20)])
                char_3 = correct_names(char_3_JINA)
                char_4_JINA = strip(dataset[(n + 21):(n + 25)])
                char_4 = correct_names(char_4_JINA)

                reaction_symbol = Symbol(char_1 * "_" * char_2 * "_to_" * char_3 * "_" * char_4)

                elem_1_u = [Symbol(char_1), Symbol(char_2)]
                elem_2_u = [Symbol(char_3), Symbol(char_4)]

                num_elem_1 = sort_reaction(elem_1_u)[2]
                num_elem_2 = sort_reaction(elem_2_u)[2]

                elem_1 = sort_reaction(elem_1_u)[1]
                elem_2 = sort_reaction(elem_2_u)[1]

                Q_value = parse(Float64, dataset[(n + 53):(n + 64)]) * Constants.MEV_TO_ERGS

                reaction_info = JinaReactionRate(reaction_symbol, elem_1, num_elem_1, elem_2, num_elem_2, Q_value, a,
                                                 set_label, res_rate, rev_rate, chap)
                add_to_references(dictionary, reference_dictionary, reaction_symbol, reaction_info)

            elseif chap == 6
                char_1_JINA = strip(dataset[(n + 6):(n + 10)])
                char_1 = correct_names(char_1_JINA)
                char_2_JINA = strip(dataset[(n + 11):(n + 15)])
                char_2 = correct_names(char_2_JINA)
                char_3_JINA = strip(dataset[(n + 16):(n + 20)])
                char_3 = correct_names(char_3_JINA)
                char_4_JINA = strip(dataset[(n + 21):(n + 25)])
                char_4 = correct_names(char_4_JINA)
                char_5_JINA = strip(dataset[(n + 26):(n + 30)])
                char_5 = correct_names(char_5_JINA)

                reaction_symbol = Symbol(char_1 * "_" * char_2 * "_to_" * char_3 * "_" * char_4 * "_" * char_5)

                elem_1_u = [Symbol(char_1), Symbol(char_2)]
                elem_2_u = [Symbol(char_3), Symbol(char_4), Symbol(char_5)]

                num_elem_1 = sort_reaction(elem_1_u)[2]
                num_elem_2 = sort_reaction(elem_2_u)[2]

                elem_1 = sort_reaction(elem_1_u)[1]
                elem_2 = sort_reaction(elem_2_u)[1]

                Q_value = parse(Float64, dataset[(n + 53):(n + 64)]) * Constants.MEV_TO_ERGS

                reaction_info = JinaReactionRate(reaction_symbol, elem_1, num_elem_1, elem_2, num_elem_2, Q_value, a,
                                                 set_label, res_rate, rev_rate, chap)
                add_to_references(dictionary, reference_dictionary, reaction_symbol, reaction_info)

            elseif chap == 7
                char_1_JINA = strip(dataset[(n + 6):(n + 10)])
                char_1 = correct_names(char_1_JINA)
                char_2_JINA = strip(dataset[(n + 11):(n + 15)])
                char_2 = correct_names(char_2_JINA)
                char_3_JINA = strip(dataset[(n + 16):(n + 20)])
                char_3 = correct_names(char_3_JINA)
                char_4_JINA = strip(dataset[(n + 21):(n + 25)])
                char_4 = correct_names(char_4_JINA)
                char_5_JINA = strip(dataset[(n + 26):(n + 30)])
                char_5 = correct_names(char_5_JINA)
                char_6_JINA = strip(dataset[(n + 31):(n + 35)])
                char_6 = correct_names(char_6_JINA)

                reaction_symbol = Symbol(char_1 * "_" * char_2 * "_to_" * char_3 * "_" * char_4 * "_" * char_5 * "_" *
                                         char_6)

                elem_1_u = [Symbol(char_1), Symbol(char_2)]
                elem_2_u = [Symbol(char_3), Symbol(char_4), Symbol(char_5), Symbol(char_6)]

                num_elem_1 = sort_reaction(elem_1_u)[2]
                num_elem_2 = sort_reaction(elem_2_u)[2]

                elem_1 = sort_reaction(elem_1_u)[1]
                elem_2 = sort_reaction(elem_2_u)[1]

                Q_value = parse(Float64, dataset[(n + 53):(n + 64)]) * Constants.MEV_TO_ERGS

                reaction_info = JinaReactionRate(reaction_symbol, elem_1, num_elem_1, elem_2, num_elem_2, Q_value, a,
                                                 set_label, res_rate, rev_rate, chap)
                add_to_references(dictionary, reference_dictionary, reaction_symbol, reaction_info)

            elseif chap == 8
                char_1_JINA = strip(dataset[(n + 6):(n + 10)])
                char_1 = correct_names(char_1_JINA)
                char_2_JINA = strip(dataset[(n + 11):(n + 15)])
                char_2 = correct_names(char_2_JINA)
                char_3_JINA = strip(dataset[(n + 16):(n + 20)])
                char_3 = correct_names(char_3_JINA)
                char_4_JINA = strip(dataset[(n + 21):(n + 25)])
                char_4 = correct_names(char_4_JINA)

                reaction_symbol = Symbol(char_1 * "_" * char_2 * "_" * char_3 * "_to_" * char_4)

                elem_1_u = [Symbol(char_1), Symbol(char_2), Symbol(char_3)]
                elem_2_u = [Symbol(char_4)]

                num_elem_1 = sort_reaction(elem_1_u)[2]
                num_elem_2 = sort_reaction(elem_2_u)[2]

                elem_1 = sort_reaction(elem_1_u)[1]
                elem_2 = sort_reaction(elem_2_u)[1]

                Q_value = parse(Float64, dataset[(n + 53):(n + 64)]) * Constants.MEV_TO_ERGS

                reaction_info = JinaReactionRate(reaction_symbol, elem_1, num_elem_1, elem_2, num_elem_2, Q_value, a,
                                                 set_label, res_rate, rev_rate, chap)
                add_to_references(dictionary, reference_dictionary, reaction_symbol, reaction_info)
            end

        else
            reaction = false
            chap += 1
        end

        n += 225
    end
end

"""
    get_reaction_rate(reaction::JinaReactionRate, eos00::EOSResults{TT}, xa::AbstractVector{TT}, xa_index::Dict{Symbol,Int})

Evaluates the reaction rate, in s^{-1}g^{-1}, given an equation of state result for the relevant cell, is abundances and index array,
by computing Eqs. 1 and 2 from Cyburt+2010.
"""
function get_reaction_rate(reaction::JinaReactionRate, eos00::EOSResults{TT}, xa::AbstractVector{TT},
                           xa_index::Dict{Symbol,Int})::TT where {TT}

    # determine λ

    T_9 = (eos00.T / 1e9)
    a = reaction.coeff

    x = a[1] + a[7] * log(T_9)

    for i = 2:6
        x += a[i] * cbrt(T_9^((2(i - 1) - 5)))
    end

    λ = exp(x)

    # println("Lambda is equal to", λ)
    # determine elements and how many times they occur
    # code gives a dictionary with the elements and how many times they occur
    elements = reaction.iso_in
    N_elements = reaction.num_iso_in
    # determine all needed parameters for every element
    ν = -1
    factors = 1

    for index in eachindex(elements)
        elem = elements[index]                                       # A
        N_elem = N_elements[index]                                  # N_A
        # println("N_elem = ", N_elem)
        X_elem = xa[xa_index[elem]]                                 # X_A
        # println("X_elem = ", X_elem)
        m_elem = Chem.isotope_list[elem].mass * Constants.AMU       # m_A
        # println("m_elem = ", m_elem)
        Y_elem = X_elem / (m_elem * Constants.AVO)                   # Y_A
        # println("Y_elem = ", Y_elem)
        ν += N_elem
        # println("ν = ", ν)
        factor_elem = Y_elem^(N_elem) / factorial(N_elem)
        # println("factor_elem = ", factor_elem)
        factors *= factor_elem

    end

    # Calculate the reaction rate
    ρ = eos00.ρ
    RR = ρ^ν * λ * factors

    return RR * Constants.AVO  # times avo since Eq. 2 is molar rate, we want in g^-1 s^-1
end

# executes when laoding in ReactionRates
file_contents = open(pkgdir(Chem, "data", "Jina_reactionrates.data")) do io
    read(io, String)
end
reaction_list = Dict()
jina_references = Dict()
jina_rates = Dict()
read_dataset(file_contents, jina_rates, jina_references)
reaction_list[:jina_rates] = jina_rates
