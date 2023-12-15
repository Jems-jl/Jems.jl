# JINA database - test file

########################################
############## TEXT ECXPL ##############
########################################

# test single reaction rate:
# p + d -> he3
# https://reaclib.jinaweb.org/d(p,g)he3/de04/

# Format file:
# per chapter:
#          n    p                            wc12w     7.82300e-01          
# -6.781610e+00 0.000000e+00 0.000000e+00 0.000000e+00                      
# 0.000000e+00 0.000000e+00 0.000000e+00                                   

# The REACLIB 1 (R1) file format is the default output of the REACLIB Database.
# Every three lines in the R1 file comprises of an entry can be of two different types.
# It is either a header entry or a set entry. Every line in a REACLIB file is
# 74 characters long and is padding with spaces. Because of this a blank line is
# actually 74 space characters.


# MEANING DATA

# The first line of the set entry is as follows: --> OK

# Five(5) space characters. ok
# Six(6) fields of five(5) characters that contain the nuclides for the rate padded with spaces. ok
# Eight(8) spaces. ok 
# Four(4) characters that contains the set label. ok
# A one(1) character flag denoting:
# blank or n: A non-resonant rate.
# r: A resonant rate.
# w: A weak rate. ok
# A one(1) character flag that is set to 'v' when this is a reverse rate. ok
# The "v" flag shows reverse rates in which detailed balance (without partition functions) was used to derive the reaction rate. Therefore, rates with this flag must be corrected to include partition function modifications.
# Three(3) space characters. ok
# A twelve(12) character field that contains the Q value for the reaction. 0.00000e+00. ok
# Ten(10) space characters. ok

# The second line of the set entry contains the first four ai (a0, a1, a2, a3)
# values in that set. Each of these fields are 13 characters (0.000000e+00).
# Followed by 22 space characters.  --> OK

# The third line of the set entry contains the last three ai (a4, a5, a6) values 
# in that set. Each of these fields are 13 characters (0.000000e+00).
# Followed by 35 space characters. --> OK

# Load in file

# load in a file from my downloads


########################################
############## START CODE ##############
########################################

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

##

# open file --> Question: when this is on github, to what change the path name?
file_path = "/Users/evakuipers/Documents/Master sterrenkunde/Thesis/results11271425"
file = open(file_path, "r")
file_contents = read(file, String)
close(file)

##

# code for a single element

n = 3 * 75       # start of the selected reaction
                 # each line is 74 tokens long            

char_1 = file_contents[(n + 6) : (n + 10)]; println("character 1: $char_1")
char_2 = file_contents[(n + 11): (n + 15)]; println("character 2: $char_2")
char_3 = file_contents[(n + 16): (n + 20)]; println("character 3: $char_3")
char_4 = file_contents[(n + 21): (n + 25)]; println("character 4: $char_4")
char_5 = file_contents[(n + 26): (n + 30)]; println("character 5: $char_5")
char_6 = file_contents[(n + 31): (n + 35)]; println("character 6: $char_6")

set_label = file_contents[(n + 44): (n + 47)]; println("set label: $set_label")

res_rate  = file_contents[(n + 48)]; println("resonance rate: $res_rate")
char_flag = file_contents[(n + 49)]; println("character flag: $char_flag")

Q_value   = file_contents[(n + 53): (n + 64)]; println("Q value: $Q_value")

# second line

a0 = file_contents[(n + 76) : (n + 88) ]; println("a0: $a0")
a1 = file_contents[(n + 89) : (n + 101)]; println("a1: $a1")
a2 = file_contents[(n + 102): (n + 114)]; println("a2: $a2")
a3 = file_contents[(n + 115): (n + 127)]; println("a3: $a3")

# third line

a4 = file_contents[(n + 151): (n + 163)]; println("a4: $a4")
a5 = file_contents[(n + 164): (n + 176)]; println("a5: $a5")
a6 = file_contents[(n + 177): (n + 189)]; println("a6: $a6")                            

##

# define the start of a chapter

start_char = file_contents[(n + 1)]; print("start character: $start_char")

# when a chapter is started, the start_char is a model_number
# what we want = to have a separate condition when this token is a number than the above
# for a reaction, this is a space

if start_char == ' '
    reaction = true
else
    reaction = false
end

# if reaction = true --> doe het bovenstaande
# if reaction = false --> markeer de chapter met + 1 & ga 3 lijnen verder


###########################
###### FULL FUNCTION ######
###########################

##

# define the Jina Reaction rate structure 

struct JinaReactionRate{TT<:Real}<:ReactionRates.AbstractReactionRate
    name::Symbol
    iso_in::Vector{Symbol}
    iso_out::Vector{Symbol}
    Qvalue::TT
    coeff::Vector{TT}
    set_label::Symbol
    res_rate::Symbol
    char_flag::Symbol
end


## 



# function om de dubbele rates te identificeren

function add_to_references(main_dict, ref_dict, reaction, new_info::JinaReactionRate)

    # main_dict = de algemene dictionary met alle reaction rates in
    # ref_dict  = de dictionary waarin alle references naar de reacties instaan
    # reaction  = Symbol van de reactie die toegevoegd moet worden aan de algemene dictionary
    # new_info  = JinaReactionRate van de nieuwe reactie

    if haskey(ref_dict, reaction) # de reaction rate bestaat al

        # current_reactions = ref_dict[reaction]

        if ref_dict[reaction] == []

            cur_info = main_dict[reaction]

            cur_set_label = cur_info.set_label
            cur_res_rate  = cur_info.res_rate
            cur_char_flag = cur_info.char_flag

            new_set_label = new_info.set_label
            new_res_rate  = new_info.res_rate
            new_char_flag = new_info.char_flag

            reaction_string_cur = "$(reaction)_$(cur_set_label)_$(cur_res_rate)_$(cur_char_flag)"
            reaction_string_new = "$(reaction)_$(new_set_label)_$(new_res_rate)_$(new_char_flag)"

            reaction_symbol_cur = Symbol(replace(reaction_string_cur, ' ' => 'x'))
            reaction_symbol_new = Symbol(replace(reaction_string_new, ' ' => 'x'))

            list = []
            push!(list, reaction_symbol_cur)
            push!(list, reaction_symbol_new)
            ref_dict[reaction] = list

            # ref_dict[reaction].

            main_dict[reaction_symbol_cur] = cur_info
            main_dict[reaction_symbol_new] = new_info

        else

            new_set_label = new_info.set_label
            new_res_rate  = new_info.res_rate
            new_char_flag = new_info.char_flag

            reaction_string_new = "$(reaction)_$(new_set_label)_$(new_res_rate)_$(new_char_flag)"
            reaction_symbol_new = Symbol(replace(reaction_string_new, ' ' => '/'))

            # ref_dict[reaction].append(reaction_symbol_new)

            list = ref_dict[reaction]
            push!(list, reaction_symbol_new)
            # list.append(reaction_symbol_new)
            ref_dict[reaction] = list


            main_dict[reaction_symbol_new] = new_info

        end

    else 
        
        ref_dict[reaction]  = []        # reaction rate bestaat nog niet
        main_dict[reaction] = new_info  # JinaReactionRate(reaction_symbol, elem_1, elem_2, Q_value, a, set_label, res_rate, char_flag)

    end
    
end


##

# functie om namen aan te passen

##

function read_set(dataset, dictionary, reference_dictionary)

    chap = 0 

    for x in 0:80030    # let it go through file to count / while loop
        n = 225 * x      
    
        if dataset[(n + 1)] == ' '

            reaction = true

            set_label = Symbol(dataset[(n + 44): (n + 47)])
            res_rate  = Symbol(dataset[(n + 48)])
            char_flag = Symbol(dataset[(n + 49)])

            a0 = parse(Float64, dataset[(n + 76) : (n + 88) ])
            a1 = parse(Float64, dataset[(n + 89) : (n + 101)])
            a2 = parse(Float64, dataset[(n + 102): (n + 114)])
            a3 = parse(Float64, dataset[(n + 115): (n + 127)])
    
            a4 = parse(Float64, dataset[(n + 151): (n + 163)])
            a5 = parse(Float64, dataset[(n + 164): (n + 176)])
            a6 = parse(Float64, dataset[(n + 177): (n + 189)])

            a  = [a0, a1, a2, a3, a4, a5, a6]

            if chap == 1

                char_1_l = strip(dataset[(n + 6) : (n + 10)]); char_1 = uppercase(first(char_1_l)) * lowercase(char_1_l[2:end]) # functie voor schrijven
                char_2_l = strip(dataset[(n + 11): (n + 15)]); char_2 = uppercase(first(char_2_l)) * lowercase(char_2_l[2:end])

                reaction_symbol = Symbol(char_1 * "_to_" * char_2)

                elem_1 = [Symbol(char_1)];
                elem_2 = [Symbol(char_2)];

                Q_value = parse(Float64, dataset[(n + 53): (n+ 64)])
                
                reaction_info = JinaReactionRate(reaction_symbol, elem_1, elem_2, Q_value, a, set_label, res_rate, char_flag)
                add_to_references(dictionary, reference_dictionary, reaction_symbol, reaction_info)
                

            elseif chap == 2

                char_1_l = strip(dataset[(n + 6) : (n + 10)]); char_1 = uppercase(first(char_1_l)) * lowercase(char_1_l[2:end])
                char_2_l = strip(dataset[(n + 11): (n + 15)]); char_2 = uppercase(first(char_2_l)) * lowercase(char_2_l[2:end])
                char_3_l = strip(dataset[(n + 16): (n + 20)]); char_3 = uppercase(first(char_3_l)) * lowercase(char_3_l[2:end])
            
                reaction_symbol = Symbol(char_1 * "_to_" * char_2 * "_" * char_3)
                
                elem_1 = [Symbol(char_1)];
                elem_2 = [Symbol(char_2), Symbol(char_3)];

                Q_value = parse(Float64, dataset[(n + 53): (n+ 64)])                         
                
                reaction_info = JinaReactionRate(reaction_symbol, elem_1, elem_2, Q_value, a, set_label, res_rate, char_flag)
                add_to_references(dictionary, reference_dictionary, reaction_symbol, reaction_info)
                
            elseif chap == 3

                char_1_l = strip(dataset[(n + 6) : (n + 10)]); char_1 = uppercase(first(char_1_l)) * lowercase(char_1_l[2:end])
                char_2_l = strip(dataset[(n + 11): (n + 15)]); char_2 = uppercase(first(char_2_l)) * lowercase(char_2_l[2:end])
                char_3_l = strip(dataset[(n + 16): (n + 20)]); char_3 = uppercase(first(char_3_l)) * lowercase(char_3_l[2:end])
                char_4_l = strip(dataset[(n + 21): (n + 25)]); char_4 = uppercase(first(char_4_l)) * lowercase(char_4_l[2:end])
            
                reaction_symbol = Symbol(char_1 * "_to_" * char_2 * "_" * char_3 * "_" *char_4)
                
                elem_1 = [Symbol(char_1)];
                elem_2 = [Symbol(char_2), Symbol(char_3), Symbol(char_4)];

                Q_value = parse(Float64, dataset[(n + 53): (n+ 64)])                          
                
                reaction_info = JinaReactionRate(reaction_symbol, elem_1, elem_2, Q_value, a, set_label, res_rate, char_flag)
                add_to_references(dictionary, reference_dictionary, reaction_symbol, reaction_info)
                
            elseif chap == 4

                char_1_l = strip(dataset[(n + 6) : (n + 10)]); char_1 = uppercase(first(char_1_l)) * lowercase(char_1_l[2:end])
                char_2_l = strip(dataset[(n + 11): (n + 15)]); char_2 = uppercase(first(char_2_l)) * lowercase(char_2_l[2:end])
                char_3_l = strip(dataset[(n + 16): (n + 20)]); char_3 = uppercase(first(char_3_l)) * lowercase(char_3_l[2:end])
            
                reaction_symbol = Symbol(char_1 * "_" * char_2 * "_to_" * char_3)
                
                elem_1 = [Symbol(char_1), Symbol(char_2)];
                elem_2 = [Symbol(char_3)];

                Q_value = parse(Float64, dataset[(n + 53): (n+ 64)])                        
                
                reaction_info = JinaReactionRate(reaction_symbol, elem_1, elem_2, Q_value, a, set_label, res_rate, char_flag)
                add_to_references(dictionary, reference_dictionary, reaction_symbol, reaction_info)
                
            elseif chap == 5

                char_1_l = strip(dataset[(n + 6) : (n + 10)]); char_1 = uppercase(first(char_1_l)) * lowercase(char_1_l[2:end])
                char_2_l = strip(dataset[(n + 11): (n + 15)]); char_2 = uppercase(first(char_2_l)) * lowercase(char_2_l[2:end])
                char_3_l = strip(dataset[(n + 16): (n + 20)]); char_3 = uppercase(first(char_3_l)) * lowercase(char_3_l[2:end])
                char_4_l = strip(dataset[(n + 21): (n + 25)]); char_4 = uppercase(first(char_4_l)) * lowercase(char_4_l[2:end])
            
                reaction_symbol = Symbol(char_1 * "_" * char_2 * "_to_" * char_3 * "_" *char_4)
                
                elem_1 = [Symbol(char_1), Symbol(char_2)];
                elem_2 = [Symbol(char_3), Symbol(char_4)];

                Q_value = parse(Float64, dataset[(n + 53): (n+ 64)])                          
                
                reaction_info = JinaReactionRate(reaction_symbol, elem_1, elem_2, Q_value, a, set_label, res_rate, char_flag)
                add_to_references(dictionary, reference_dictionary, reaction_symbol, reaction_info)
                
            elseif chap == 6

                char_1_l = strip(dataset[(n + 6) : (n + 10)]); char_1 = uppercase(first(char_1_l)) * lowercase(char_1_l[2:end])
                char_2_l = strip(dataset[(n + 11): (n + 15)]); char_2 = uppercase(first(char_2_l)) * lowercase(char_2_l[2:end])
                char_3_l = strip(dataset[(n + 16): (n + 20)]); char_3 = uppercase(first(char_3_l)) * lowercase(char_3_l[2:end])
                char_4_l = strip(dataset[(n + 21): (n + 25)]); char_4 = uppercase(first(char_4_l)) * lowercase(char_4_l[2:end])
                char_5_l = strip(dataset[(n + 26): (n + 30)]); char_5 = uppercase(first(char_5_l)) * lowercase(char_5_l[2:end])

                reaction_symbol = Symbol(char_1 * "_" * char_2 * "_to_" * char_3 * "_" * char_4 * "_" * char_5)
                
                elem_1 = [Symbol(char_1), Symbol(char_2)];
                elem_2 = [Symbol(char_3), Symbol(char_4), Symbol(char_5)];

                Q_value = parse(Float64, dataset[(n + 53): (n+ 64)])                        
                
                reaction_info = JinaReactionRate(reaction_symbol, elem_1, elem_2, Q_value, a, set_label, res_rate, char_flag)
                add_to_references(dictionary, reference_dictionary, reaction_symbol, reaction_info)
                
            elseif chap == 7

                char_1_l = strip(dataset[(n + 6) : (n + 10)]); char_1 = uppercase(first(char_1_l)) * lowercase(char_1_l[2:end])
                char_2_l = strip(dataset[(n + 11): (n + 15)]); char_2 = uppercase(first(char_2_l)) * lowercase(char_2_l[2:end])
                char_3_l = strip(dataset[(n + 16): (n + 20)]); char_3 = uppercase(first(char_3_l)) * lowercase(char_3_l[2:end])
                char_4_l = strip(dataset[(n + 21): (n + 25)]); char_4 = uppercase(first(char_4_l)) * lowercase(char_4_l[2:end])
                char_5_l = strip(dataset[(n + 26): (n + 30)]); char_5 = uppercase(first(char_5_l)) * lowercase(char_5_l[2:end])
                char_6_l = strip(dataset[(n + 31): (n + 35)]); char_6 = uppercase(first(char_6_l)) * lowercase(char_6_l[2:end])

                reaction_symbol = Symbol(char_1 * "_" * char_2 * "_to_" * char_3 * "_" * char_4 * "_" * char_5 * "_" * char_6)
                
                elem_1 = [Symbol(char_1), Symbol(char_2)];
                elem_2 = [Symbol(char_3), Symbol(char_4), Symbol(char_5), Symbol(char_6)];

                Q_value = parse(Float64, dataset[(n + 53): (n+ 64)])                       
                
                reaction_info = JinaReactionRate(reaction_symbol, elem_1, elem_2, Q_value, a, set_label, res_rate, char_flag)
                add_to_references(dictionary, reference_dictionary, reaction_symbol, reaction_info)
                
            elseif chap == 8

                char_1_l = strip(dataset[(n + 6) : (n + 10)]); char_1 = uppercase(first(char_1_l)) * lowercase(char_1_l[2:end])
                char_2_l = strip(dataset[(n + 11): (n + 15)]); char_2 = uppercase(first(char_2_l)) * lowercase(char_2_l[2:end])
                char_3_l = strip(dataset[(n + 16): (n + 20)]); char_3 = uppercase(first(char_3_l)) * lowercase(char_3_l[2:end])
                char_4_l = strip(dataset[(n + 21): (n + 25)]); char_4 = uppercase(first(char_4_l)) * lowercase(char_4_l[2:end])

                reaction_symbol = Symbol(char_1 * "_" * char_2 * "_" * char_3 * "_to_" * char_4)
                
                elem_1 = [Symbol(char_1), Symbol(char_2), Symbol(char_3)];
                elem_2 = [Symbol(char_4)];

                Q_value = parse(Float64, dataset[(n + 53): (n+ 64)])                       
                
                reaction_info = JinaReactionRate(reaction_symbol, elem_1, elem_2, Q_value, a, set_label, res_rate, char_flag)
                add_to_references(dictionary, reference_dictionary, reaction_symbol, reaction_info)
                
            end
            
        else

            reaction = false
            chap += 1

        end
    end
        
end

##

# p to H1
# n stays small
# Deuterium / Tritium --> H2 & H3

References = Dict()
Jina_Rates = Dict()
read_set(file_contents, Jina_Rates, References)

## 

References = Dict()
Jina_Rates = Dict()

##

using BenchmarkTools
@benchmark read_set(file_contents, Jina_Rates, References)

##

print(Jina_Rates)

##

# print(References)

##

for (key, value) in References
    println("$key: $value")
end


##

# lijst niet leeg laten: ook de enkele rates erin zetten
# definitie van reaction rate bekijken



# Reaction Rates


function reaction_rate_Jina(a_0, a, eos00::EOSResults{TT})::TT where{TT}
    T_9 = (eos00.T / 1e9)
    sum_term = sum(i -> a[i] * T_9^((2i - 5)/3), 1:5)
    x = exp(a_0 + sum_term + a[6] * log(T_9))
    return x
end

##

a_0 = 1.0
a = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
# T_9 = 3.0

r = EOSResults{Float64}()
eos = EOS.IdealEOS(false)

set_EOS_resultsTρ!(eos, r, log(3.0e9), log(1), [1.0, 0.0], [:H1, :He4])

result = reaction_rate_Jina(a_0, a, r);
println(result)

##

function get_reaction_rate(reaction::JinaReactionRate, eos00::EOSResults{TT}, xa::AbstractVector{TT}, xa_index::Dict{Symbol,Int})::TT where{TT}

    T_9 = (eos00.T / 1e9)
    a = reaction.coeff

    # sum_term = sum(i -> a[i] * T_9^((2i - 5)/3), 2:6)
    λ = a[1] +  a[7] * log(T_9)

    for i in 2:6
        λ += a[i] * T_9^((2(i-1) - 5)/3)
        
    end

    # λ = exp(a[1] + sum_term + a[7] * log(T_9))

    return exp(λ)    

end


##

r = EOSResults{Float64}()
eos = EOS.IdealEOS(false)
logT = LinRange(7.0,10.0,100)
rates = zeros(100)

for i in eachindex(logT)
    # set_EOS_resultsTρ!(eos, r, logT[i]*log(10), log(1), [1.0, 0.0], [:H1, :He4])
    r.T = 10^(logT[i])
    rates[i] = get_reaction_rate(Jina_Rates[:P_P_to_D], r,[1.0, 0.0], Dict{Symbol, Int}() )
 
end

##

rates
# Jina_Rates[:P_D_to_He3].coeff

##

using CairoMakie

f = Figure();
ax = Axis(f[1, 1]; xlabel=L"\log_{10}(T_\mathrm{eff}/[K])", ylabel=L"\log_{10}(L/L_\odot)", xreversed=false)
# history = StellarModels.get_history_dataframe_from_hdf5("history.hdf5")
lines!(ax, logT, log10.(rates))
# ylims!(ax, -3,5)
f

##
# extra dict with all references (how many rates of that type)
# add returns to kipp_rates



