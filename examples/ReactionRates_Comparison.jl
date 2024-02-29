##

using BenchmarkTools
using Jems.Chem
using Jems.Constants
using Jems.EOS
using Jems.Opacity
using Jems.NuclearNetworks
using Jems.StellarModels
using Jems.Evolution
using Jems.ReactionRates

### Kippenhahn Rates ###

struct KippReactionRate{TT<:Real}<:ReactionRates.AbstractReactionRate
    name::Symbol
    iso_in::Vector{Symbol}
    num_iso_in::Vector{Int64}
    iso_out::Vector{Symbol}
    num_iso_out::Vector{Int64}
    Qvalue::TT
end
reaction_list = Dict()
reaction_list[:kipp_rates] = Dict(
    :kipp_pp => KippReactionRate(:kipp_pp, [:H1], [4], [:He4], [1],
        ((4 * Chem.isotope_list[:H1].mass - Chem.isotope_list[:He4].mass) * AMU * CLIGHT^2)),
        
    :kipp_cno => KippReactionRate(:kipp_cno, [:H1], [4], [:He4], [1],
        ((4 * Chem.isotope_list[:H1].mass - Chem.isotope_list[:He4].mass) * AMU * CLIGHT^2)),

    :kipp_3alphaCF88 => KippReactionRate(:kipp_3alphaCF88, [:He4], [3], [:C12], [1],
        ((3 * Chem.isotope_list[:He4].mass - Chem.isotope_list[:C12].mass) * AMU * CLIGHT^2)),

    :kipp_3alphaA99 => KippReactionRate(:kipp_3alphaA99, [:He4], [3], [:C12], [1],
        ((3 * Chem.isotope_list[:He4].mass - Chem.isotope_list[:C12].mass) * AMU * CLIGHT^2)),

    :kipp_C12alpha => KippReactionRate(:kipp_C12alpha, [:C12, :He4], [1,1], [:O16], [1],
        ((1 * Chem.isotope_list[:He4].mass + 1 * Chem.isotope_list[:C12].mass - Chem.isotope_list[:O16].mass) * AMU * CLIGHT^2)),

    :kipp_O16alpha => KippReactionRate(:kipp_O16alpha, [:O16, :He4], [1,1], [:Ne20], [1],
        ((Chem.isotope_list[:He4].mass + Chem.isotope_list[:O16].mass - Chem.isotope_list[:Ne20].mass) * AMU * CLIGHT^2)),

    :kipp_CC => KippReactionRate(:kipp_CC, [:C12], [2], [:O16, :He4], [1,2],
        ((2 * Chem.isotope_list[:C12].mass - Chem.isotope_list[:O16].mass - 2 * Chem.isotope_list[:He4].mass) * AMU * CLIGHT^2)),

    :kipp_OO => KippReactionRate(:kipp_OO, [:O16], [2], [:Mg24, :He4], [1,2],
        ((2 * Chem.isotope_list[:O16].mass - Chem.isotope_list[:Mg24].mass - 2 * Chem.isotope_list[:He4].mass) * AMU * CLIGHT^2))   
)
function get_reaction_rate_Kipp(reaction::KippReactionRate, eos00::EOSResults{TT}, xa::AbstractVector{TT}, xa_index::Dict{Symbol,Int})::TT where{TT}
    
    """
    Input:
    reaction: the reactions dictionary that is being used
    eos00: results equation of state
    xa: element fractions in the star
    xa_index: index of the elements

    Output:
    for each reaction, the ϵnuc value is calculated in the function.
        
    """
    
    if reaction.name == :kipp_pp 

        phi  = 1
        f_11 = 1
        T9   = (eos00.T / 1e9)
        X1   = xa[xa_index[:H1]]

        g_11 = (1 + 3.82 * T9 + 1.51 * T9^2 + 0.144 * T9^3 - 0.0114 * T9^4)
        ϵnuc = 2.57e4 * phi * f_11 * g_11 * eos00.ρ * X1^2 * cbrt(T9^(-2)) * exp(-3.381 * cbrt(T9^(-1)))

        return ϵnuc / reaction.Qvalue

    elseif reaction.name == :kipp_cno 

        T9    = (eos00.T / 1e9)
        X1    = xa[xa_index[:H1]]
        X_CNO = xa[xa_index[:C12]] + xa[xa_index[:N14]] + xa[xa_index[:O16]]

        g_14  = (1 - 2.00 * T9 + 3.41 * T9^2 - 2.43 * T9^3 )
        ϵnuc  = 8.24e25 * g_14 * X_CNO * X1 * eos00.ρ * 
                cbrt(T9^(-2)) * exp(-15.231 * cbrt(T9^(-1)) - (T9/0.8)^2)

        return ϵnuc / reaction.Qvalue

    elseif reaction.name == :kipp_3alphaCF88

        f_3alpha = 1
        X4   = xa[xa_index[:He4]]
        T8   = eos00.T / 1e8

        ϵnuc = 5.09e11 * f_3alpha * (eos00.ρ)^2 * X4^3 * 
               T8^-3 * exp(-44.027 / T8)

    elseif reaction.name == :kipp_3alphaA99

        f_3alpha = 1
        X4   = xa[xa_index[:He4]]
        T9   = eos00.T / 1e9

        ϵnuc = 6.272*(eos00.ρ)^2*X4^3*(1+0.0158*T9^(-0.65))*
                (2.43e9*cbrt(T9^(-2))*exp(-13.490*cbrt(T9^(-1))-(T9/0.15)^2)*(1+74.5*T9)
                    + 6.09e5*sqrt(T9^(-3))*exp(-1.054/T9))*
                (2.76e7*cbrt(T9^(-2))*exp(-23.570*cbrt(T9^(-1))-(T9/0.4)^2)
                    * (1+5.47*T9+326*T9^2) + 130.7*sqrt(T9^(-3))*exp(-3.338/T9)
                    + 2.51e4*sqrt(T9^(-3))*exp(-20.307/T9))

    elseif reaction.name == :kipp_C12alpha

        f_12alpha = 1
        X4   = xa[xa_index[:He4]]
        X12  = xa[xa_index[:C12]]
        T8   = eos00.T / 1e8
        
        ϵnuc = 1.3e27 * f_12alpha * eos00.ρ * X4 * X12 * T8^(-2) *
                ((1 + 0.134 * cbrt(T8^(2)))/(1 + 0.017 * cbrt(T8^(2))))^2 * exp(-69.20 / cbrt(T8))

    elseif reaction.name == :kipp_O16alpha

        f_16alpha = 1
        X4   = xa[xa_index[:He4]]
        X16  = xa[xa_index[:O16]]
        T9   = eos00.T / 1e9
                
        #ϵnuc = 1.91e27 * cbrt(T9^(-2)) * X16 * X4 * eos00.ρ * f_16alpha *
        #        exp(-39.76 * T9^(-1/3) - (T9/1.6)^2) 
        #        + 3.64 * 10^18 * sqrt(T9^(-3)) * exp(-10.32 / T9)
        #        + 4.39 * 10^19 * sqrt(T9^(-3)) * exp(-12.20 / T9)
        #        + 2.92 * 10^16 * T9^(2.966) * exp(-11.90 / T9)
        # Pablo: I think the above one has a typo in the Kippenhahn book itself,
        # otherwise it does not make sense that part of the rate is independent
        # of density and mass fractions. I think the one below would make sense.
        ϵnuc = X16 * X4 * eos00.ρ * f_16alpha * (
                1.91e27 * cbrt(T9^(-2)) * exp(-39.76 * T9^(-1/3) - (T9/1.6)^2) 
                + 3.64e18 * sqrt(T9^(-3)) * exp(-10.32 / T9)
                + 4.39e19 * sqrt(T9^(-3)) * exp(-12.20 / T9)
                + 2.92e16 * T9^(2.966) * exp(-11.90 / T9)
            )


    elseif reaction.name == :kipp_CC

        f_CC = 1
        T9   = (eos00.T / 1e9)
        T_9a = T9 / (1 + 0.0396 * T9) 
        X12  = xa[xa_index[:C12]]

        ϵnuc = 1.86e43 * f_CC * eos00.ρ * X12^2 * sqrt(T9^(-3)) * T_9a^(5/6) *
               exp(-84.165 / cbrt(T_9a) - 2.12e-3 * T9^3)

    elseif reaction.name == :kipp_OO

        f_OO = 1
        T9   = (eos00.T / 1e9)
        X16  = xa[xa_index[:O16]]

        exp_func = -135.93 / cbrt(T9) - 0.629 * cbrt(T9^(2)) - 0.445 * cbrt(T9^(4)) + 0.0103 * T9^2
        ϵnuc = 2.14e53 * f_OO * eos00.ρ * X16^2 * cbrt(T9^(-2)) * exp(exp_func)

    else
        throw(ArgumentError("No method to compute rate for $(reaction.name)"))
    end
end


##

### Jina Rates ###

# open file --> Question: when this is on github, to what change the path name?
file_path = "/Users/evakuipers/Documents/Master sterrenkunde/Thesis/results11271425"
file = open(file_path, "r")
file_contents = read(file, String)
close(file)
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
function correct_names(JINA_name)
    change_name = Dict("p" => "H1", "d" => "D2", "t" => "T3", "n" => "n")

    if haskey(change_name, JINA_name)
        RETURN_name = change_name[JINA_name]
    else
        RETURN_name = uppercase(JINA_name[1]) * lowercase(JINA_name[2:end])
    end

    return RETURN_name
end
function read_set(dataset, dictionary, reference_dictionary)

    chap = 0 
    n = 0

    while n <= lastindex(dataset) - 225       
    
        if dataset[(n + 1)] == ' '

            reaction = true

            set_label = Symbol(dataset[(n + 44): (n + 47)])
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

            if chap == 1

                char_1_JINA = strip(dataset[(n + 6) : (n + 10)]); char_1 = correct_names(char_1_JINA) 
                char_2_JINA = strip(dataset[(n + 11): (n + 15)]); char_2 = correct_names(char_2_JINA) 

                reaction_symbol = Symbol(char_1 * "_to_" * char_2)

                elem_1 = [Symbol(char_1)];
                elem_2 = [Symbol(char_2)];

                Q_value = parse(Float64, dataset[(n + 53): (n+ 64)])
                
                reaction_info = JinaReactionRate(reaction_symbol, elem_1, elem_2, Q_value, a, set_label, res_rate, rev_rate, chap)
                add_to_references(dictionary, reference_dictionary, reaction_symbol, reaction_info)
                

            elseif chap == 2

                char_1_JINA = strip(dataset[(n + 6) : (n + 10)]); char_1 = correct_names(char_1_JINA) 
                char_2_JINA = strip(dataset[(n + 11): (n + 15)]); char_2 = correct_names(char_2_JINA) 
                char_3_JINA = strip(dataset[(n + 16): (n + 20)]); char_3 = correct_names(char_3_JINA) 
            
                reaction_symbol = Symbol(char_1 * "_to_" * char_2 * "_" * char_3)
                
                elem_1 = [Symbol(char_1)];
                elem_2 = [Symbol(char_2), Symbol(char_3)];

                Q_value = parse(Float64, dataset[(n + 53): (n+ 64)])                         
                
                reaction_info = JinaReactionRate(reaction_symbol, elem_1, elem_2, Q_value, a, set_label, res_rate, rev_rate, chap)
                add_to_references(dictionary, reference_dictionary, reaction_symbol, reaction_info)
                
            elseif chap == 3

                char_1_JINA = strip(dataset[(n + 6) : (n + 10)]); char_1 = correct_names(char_1_JINA)
                char_2_JINA = strip(dataset[(n + 11): (n + 15)]); char_2 = correct_names(char_2_JINA)
                char_3_JINA = strip(dataset[(n + 16): (n + 20)]); char_3 = correct_names(char_3_JINA)
                char_4_JINA = strip(dataset[(n + 21): (n + 25)]); char_4 = correct_names(char_4_JINA)
            
                reaction_symbol = Symbol(char_1 * "_to_" * char_2 * "_" * char_3 * "_" *char_4)
                
                elem_1 = [Symbol(char_1)];
                elem_2 = [Symbol(char_2), Symbol(char_3), Symbol(char_4)];

                Q_value = parse(Float64, dataset[(n + 53): (n+ 64)])                          
                
                reaction_info = JinaReactionRate(reaction_symbol, elem_1, elem_2, Q_value, a, set_label, res_rate, rev_rate, chap)
                add_to_references(dictionary, reference_dictionary, reaction_symbol, reaction_info)
                
            elseif chap == 4

                char_1_JINA = strip(dataset[(n + 6) : (n + 10)]); char_1 = correct_names(char_1_JINA)
                char_2_JINA = strip(dataset[(n + 11): (n + 15)]); char_2 = correct_names(char_2_JINA)
                char_3_JINA = strip(dataset[(n + 16): (n + 20)]); char_3 = correct_names(char_3_JINA)
            
                reaction_symbol = Symbol(char_1 * "_" * char_2 * "_to_" * char_3)
                
                elem_1 = [Symbol(char_1), Symbol(char_2)];
                elem_2 = [Symbol(char_3)];

                Q_value = parse(Float64, dataset[(n + 53): (n+ 64)])                        
                
                reaction_info = JinaReactionRate(reaction_symbol, elem_1, elem_2, Q_value, a, set_label, res_rate, rev_rate, chap)
                add_to_references(dictionary, reference_dictionary, reaction_symbol, reaction_info)
                
            elseif chap == 5

                char_1_JINA = strip(dataset[(n + 6) : (n + 10)]); char_1 = correct_names(char_1_JINA)
                char_2_JINA = strip(dataset[(n + 11): (n + 15)]); char_2 = correct_names(char_2_JINA)
                char_3_JINA = strip(dataset[(n + 16): (n + 20)]); char_3 = correct_names(char_3_JINA)
                char_4_JINA = strip(dataset[(n + 21): (n + 25)]); char_4 = correct_names(char_4_JINA)
            
                reaction_symbol = Symbol(char_1 * "_" * char_2 * "_to_" * char_3 * "_" *char_4)
                
                elem_1 = [Symbol(char_1), Symbol(char_2)];
                elem_2 = [Symbol(char_3), Symbol(char_4)];

                Q_value = parse(Float64, dataset[(n + 53): (n+ 64)])                          
                
                reaction_info = JinaReactionRate(reaction_symbol, elem_1, elem_2, Q_value, a, set_label, res_rate, rev_rate, chap)
                add_to_references(dictionary, reference_dictionary, reaction_symbol, reaction_info)
                
            elseif chap == 6

                char_1_JINA = strip(dataset[(n + 6) : (n + 10)]); char_1 = correct_names(char_1_JINA)
                char_2_JINA = strip(dataset[(n + 11): (n + 15)]); char_2 = correct_names(char_2_JINA)
                char_3_JINA = strip(dataset[(n + 16): (n + 20)]); char_3 = correct_names(char_3_JINA)
                char_4_JINA = strip(dataset[(n + 21): (n + 25)]); char_4 = correct_names(char_4_JINA)
                char_5_JINA = strip(dataset[(n + 26): (n + 30)]); char_5 = correct_names(char_5_JINA)

                reaction_symbol = Symbol(char_1 * "_" * char_2 * "_to_" * char_3 * "_" * char_4 * "_" * char_5)
                
                elem_1 = [Symbol(char_1), Symbol(char_2)];
                elem_2 = [Symbol(char_3), Symbol(char_4), Symbol(char_5)];

                Q_value = parse(Float64, dataset[(n + 53): (n+ 64)])                        
                
                reaction_info = JinaReactionRate(reaction_symbol, elem_1, elem_2, Q_value, a, set_label, res_rate, rev_rate, chap)
                add_to_references(dictionary, reference_dictionary, reaction_symbol, reaction_info)
                
            elseif chap == 7

                char_1_JINA = strip(dataset[(n + 6) : (n + 10)]); char_1 = correct_names(char_1_JINA)
                char_2_JINA = strip(dataset[(n + 11): (n + 15)]); char_2 = correct_names(char_2_JINA)
                char_3_JINA = strip(dataset[(n + 16): (n + 20)]); char_3 = correct_names(char_3_JINA)
                char_4_JINA = strip(dataset[(n + 21): (n + 25)]); char_4 = correct_names(char_4_JINA)
                char_5_JINA = strip(dataset[(n + 26): (n + 30)]); char_5 = correct_names(char_5_JINA)
                char_6_JINA = strip(dataset[(n + 31): (n + 35)]); char_6 = correct_names(char_6_JINA)

                reaction_symbol = Symbol(char_1 * "_" * char_2 * "_to_" * char_3 * "_" * char_4 * "_" * char_5 * "_" * char_6)
                
                elem_1 = [Symbol(char_1), Symbol(char_2)];
                elem_2 = [Symbol(char_3), Symbol(char_4), Symbol(char_5), Symbol(char_6)];

                Q_value = parse(Float64, dataset[(n + 53): (n+ 64)])                       
                
                reaction_info = JinaReactionRate(reaction_symbol, elem_1, elem_2, Q_value, a, set_label, res_rate, rev_rate, chap)
                add_to_references(dictionary, reference_dictionary, reaction_symbol, reaction_info)
                
            elseif chap == 8

                char_1_JINA = strip(dataset[(n + 6) : (n + 10)]); char_1 = correct_names(char_1_JINA)
                char_2_JINA = strip(dataset[(n + 11): (n + 15)]); char_2 = correct_names(char_2_JINA)
                char_3_JINA = strip(dataset[(n + 16): (n + 20)]); char_3 = correct_names(char_3_JINA)
                char_4_JINA = strip(dataset[(n + 21): (n + 25)]); char_4 = correct_names(char_4_JINA)

                reaction_symbol = Symbol(char_1 * "_" * char_2 * "_" * char_3 * "_to_" * char_4)
                
                elem_1 = [Symbol(char_1), Symbol(char_2), Symbol(char_3)];
                elem_2 = [Symbol(char_4)];

                Q_value = parse(Float64, dataset[(n + 53): (n+ 64)])                       
                
                reaction_info = JinaReactionRate(reaction_symbol, elem_1, elem_2, Q_value, a, set_label, res_rate, rev_rate, chap)
                add_to_references(dictionary, reference_dictionary, reaction_symbol, reaction_info)
                
            end
            
        else

            reaction = false
            chap += 1

        end

        n += 225

    end
        
end
function get_reaction_rate_Jina(reaction::JinaReactionRate, eos00::EOSResults{TT}, xa::AbstractVector{TT}, xa_index::Dict{Symbol,Int})::TT where{TT}

    # determine λ

    T_9 = (eos00.T / 1e9)
    a = reaction.coeff

    x = a[1] +  a[7] * log(T_9)

    for i in 2:6
        x += a[i] * T_9^((2(i-1) - 5)/3)
        
    end

    λ = exp(x)


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

    ρ = eos00.ρ
    RR = ρ^ν * λ

    for factor in factors
        RR = RR * factor

    end

    return RR

end



##

# TEST #

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
    # rates_Kipp_C12log10.(rates_Kipp_O16alpha)alpha[i] = get_reaction_rate_Kipp(reaction_list[:kipp_rates][:kipp_C12alpha], r, [0.2, 0.8, 0.0], Dict{Symbol, Int64}(:He4 => 1, :C12 => 2, :O16 => 3))
    # rates_Jina_C12alpha[i] = get_reaction_rate_Jina(Jina_Rates[:He4_C12_to_O16], r, [0.2, 0.8, 0.0], Dict{Symbol, Int64}(:He4 => 1, :C12 => 2, :O16 => 3))                            

    # 016 alpha

    set_EOS_resultsTρ!(eos, r, logT[i]*log(10), log(1), [0.2, 0.8, 0.0], [:O16, :He4, :Ne20])              
    rates_Kipp_O16alpha[i] = get_reaction_rate_Kipp(reaction_list[:kipp_rates][:kipp_O16alpha], r, [0.2, 0.8, 0.0], Dict{Symbol, Int64}(:O16 => 1, :He4 => 2, :Ne20 => 3))
    rates_Jina_O16alpha[i] = get_reaction_rate_Jina(Jina_Rates[:He4_O16_to_Ne20], r, [0.2, 0.8, 0.0], Dict{Symbol, Int64}(:O16 => 1, :He4 => 2, :Ne20 => 3))                            

    # Triple alpha

    # set_EOS_resultsTρ!(eos, r, logT[i]*log(10), log(1), [0.2, 0.8, 0.0], [:C12, :He4])              
    # rates_Kipp_3alphaCF88[i] = get_reaction_rate_Kipp(reaction_list[:kipp_rates][:kipp_3alphaCF88], r, [0.2, 0.8], Dict{Symbol, Int64}(:He4 => 1, :C12 => 2))
    # rates_Kipp_3alphaA99[i]  = get_reaction_rate_Kipp(reaction_list[:kipp_rates][:kipp_3alphaA99],  r, [0.2, 0.8], Dict{Symbol, Int64}(:He4 => 1, :C12 => 2))
    # rates_Jina_3alpha[i] = get_reaction_rate_Jina(Jina_Rates[:He4_He4_He4_to_C12], r, [0.2, 0.8], Dict{Symbol, Int64}(:He4 => 1, :C12 => 2))                            


end


##

rates_Kipp_C12alpha .= rates_Kipp_C12alpha .* (Constants.MEV_TO_ERGS)^(-1)
rates_Jina_C12alpha .= rates_Jina_C12alpha .* (Constants.AVO)

rates_Kipp_O16alpha .= rates_Kipp_O16alpha .* (Constants.MEV_TO_ERGS)^(-1)
rates_Jina_O16alpha .= rates_Jina_O16alpha .* (Constants.AVO)

rates_Kipp_3alphaCF88 .= rates_Kipp_3alphaCF88 .* (Constants.MEV_TO_ERGS)^(-1)
rates_Kipp_3alphaA99 .= rates_Kipp_3alphaA99 .* (Constants.MEV_TO_ERGS)^(-1)
rates_Jina_3alpha .= rates_Jina_3alpha .* (Constants.AVO)



##



using CairoMakie

f = Figure();
ax = Axis(f[1, 1]; xlabel=L"\log_{10}(T_\mathrm{eff}/[K])", ylabel="Reaction rate Kipp [mol / g * s]", xreversed=false)
# lines!(ax, logT, log10.(rates_Kipp_3alphaCF88), label = "rates_Kipp_3alphaCF88")
# lines!(ax, logT, log10.(rates_Kipp_3alphaA99), label = "rates_Kipp_3alphaA99")
# lines!(ax, logT, log10.(rates_Jina_3alpha) .+ 0.9, label = "rates_Jina_3alpha")

# lines!(ax, logT, log10.(rates_Kipp_O16alpha), label = "rates_Kipp_O16alpha")
# lines!(ax, logT, log10.(rates_Jina_O16alpha) .+ 0.9, label = "rates_Jina_O16alpha")

lines!(ax, logT, log10.(rates_Kipp_C12alpha), label = "rates_Kipp_C12log10")
lines!(ax, logT, log10.(rates_Jina_C12alpha) .+ 0.9, label = "rates_Jina_C12alpha")

axislegend(ax; position=:rb)
f

##





