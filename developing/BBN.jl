using Jems.Evolution
using Jems.StellarModels
using Jems.NuclearNetworks
using Jems.ReactionRates
using Jems.Constants
using Jems.Chem
using Jems.DualSupport
using Jems.Plotting

##
bbrates = [(:jina_rates, :n_to_H1_wc12_w_x_0),
           (:jina_rates, :n_H1_to_D2_an06_n_x_0),
           (:jina_rates, :n_H1_to_D2_an06_n_x_1),
           (:jina_rates, :n_H1_to_D2_an06_n_x_2),
           (:jina_rates, :D2_to_n_H1_an06_n_v_0),
           (:jina_rates, :D2_to_n_H1_an06_n_v_1),
           (:jina_rates, :D2_to_n_H1_an06_n_v_2),
           (:jina_rates, :H1_D2_to_He3_de04_n_x_0),
           (:jina_rates, :H1_D2_to_He3_de04_x_x_0),
           (:jina_rates, :He3_to_H1_D2_de04_n_v_0),
           (:jina_rates, :He3_to_H1_D2_de04_x_v_0),
           (:jina_rates, :D2_D2_to_H1_T3_go17_n_x_0),
           (:jina_rates, :H1_T3_to_D2_D2_go17_n_v_0),
           (:jina_rates, :D2_D2_to_n_He3_gi17_n_x_0),
           (:jina_rates, :n_He3_to_D2_D2_gi17_n_v_0),
           (:jina_rates, :D2_T3_to_n_He4_de04_x_x_0),
           (:jina_rates, :D2_T3_to_n_He4_de04_x_x_1),
           (:jina_rates, :n_He4_to_D2_T3_de04_x_v_0),
           (:jina_rates, :n_He4_to_D2_T3_de04_x_v_1),
           (:jina_rates, :D2_He3_to_H1_He4_de04_x_x_0),
           (:jina_rates, :D2_He3_to_H1_He4_de04_x_x_1),
           (:jina_rates, :H1_He4_to_D2_He3_de04_x_v_0),
           (:jina_rates, :H1_He4_to_D2_He3_de04_x_v_1),
           (:jina_rates, :He4_T3_to_Li7_de04_x_x_0),
           (:jina_rates, :Li7_to_He4_T3_de04_x_v_0),
           (:jina_rates, :He4_He3_to_Be7_cd08_n_x_0),
           (:jina_rates, :He4_He3_to_Be7_cd08_n_x_1),
           (:jina_rates, :Be7_to_He4_He3_cd08_n_v_0),
           (:jina_rates, :Be7_to_He4_He3_cd08_n_v_1),
           (:jina_rates, :n_He3_to_H1_T3_de04_x_x_0),
           (:jina_rates, :n_He3_to_H1_T3_de04_x_x_1),
           (:jina_rates, :H1_T3_to_n_He3_de04_x_v_0),
           (:jina_rates, :H1_T3_to_n_He3_de04_x_v_1),
           (:jina_rates, :n_Be7_to_H1_Li7_db18_x_x_0),
           (:jina_rates, :H1_Li7_to_n_Be7_db18_x_v_0),
           (:jina_rates, :H1_Li7_to_He4_He4_de04_x_x_0),
           (:jina_rates, :H1_Li7_to_He4_He4_de04_r_x_0),
           (:jina_rates, :H1_Li7_to_He4_He4_de04_x_x_1),
           (:jina_rates, :H1_Li7_to_He4_He4_de04_r_x_1),
           (:jina_rates, :He4_He4_to_H1_Li7_de04_x_v_0),
           (:jina_rates, :He4_He4_to_H1_Li7_de04_r_v_0),
           (:jina_rates, :He4_He4_to_H1_Li7_de04_x_v_1),
           (:jina_rates, :He4_He4_to_H1_Li7_de04_r_v_1)]
net = NuclearNetwork([:n, :H1, :D2, :T3, :He3, :He4, :Be7, :Li7],
                     bbrates)

function equation_composition(oz::OneZone, k::Int, iso_name::Symbol)  # needs this signature for TypeStableEquations
    # Get mass fraction for this iso
    X00 = get_00_dual(oz.props.xa[oz.network.xa_index[iso_name]])
    dXdt_nuc::typeof(X00) = 0
    reactions_in = oz.network.species_reactions_in[oz.network.xa_index[iso_name]]
    for reaction_in in reactions_in
        rate = get_00_dual(oz.props.rates[reaction_in[1]])
        dXdt_nuc -= rate * reaction_in[2] * Chem.isotope_list[iso_name].A * AMU
    end
    reactions_out = oz.network.species_reactions_out[oz.network.xa_index[iso_name]]
    for reaction_out in reactions_out
        rate = get_00_dual(oz.props.rates[reaction_out[1]])
        dXdt_nuc += rate * reaction_out[2] * Chem.isotope_list[iso_name].A * AMU
    end
    Xi = get_value(oz.prv_step_props.xa[oz.network.xa_index[iso_name]])  # is never a dual!!
    return ((X00 - Xi) / oz.props.dt - dXdt_nuc)
end

##
function T(oz::OneZone)
    return 0.1 * Constants.MEV_TO_ERGS / (KERG * sqrt(oz.props.time / 132))
end

t₀ = 13.797e9 * SECYEAR
T_CMB = 2.72548  # K
a_eq = 1.0 / (1 + 3387)
t_eq = t₀ * (a_eq)^1.5

function T_new(oz::OneZone)
    return T_CMB * sqrt(t_eq / oz.props.time) / a_eq
end

Ω_b = 0.02237
ρ_today = Ω_b * 1.8788e-29  # g/cm^3
function ρ_new(oz::OneZone)
    return ρ_today * (t_eq / oz.props.time)^(1.5) / (a_eq^3)
end

conversion = KERG^4 / (Constants.HBAR^3 * CLIGHT^5)
function ρ(oz::OneZone)
    photondensity = pi^2 / 15 * (oz.props.T)^4 * conversion
    return photondensity * 6.23e-10  # baryon density
end

oz = OneZone(equation_composition, net);

## IC
oz.props.time = 5.0  # s
oz.props.T = T_new(oz)  # K
oz.props.ρ = ρ_new(oz)  # g cm^{-3}
oz.props.ind_vars = zeros(oz.network.nspecies)
oz.props.ind_vars[oz.network.xa_index[:H1]] = 0.85
oz.props.ind_vars[oz.network.xa_index[:n]] = 0.15
# oz.props.ind_vars[oz.network.xa_index[:n]] = 1.0
oz.props.dt_next = 1.0
oz.props.model_number = 0

# ## first eval setup
# StellarModels.evaluate_model_properties!(oz, oz.props);
# StellarModels.cycle_props!(oz);
# oz.props = deepcopy(oz.prv_step_props);
# oz.props.dt = oz.prv_step_props.dt_next

# ##
# get_value.(oz.props.xa)

# ##
# Evolution.eval_jacobian_eqs!(oz)

# ##
# oz.solver_data.jacobian_D[1]
# oz.solver_data.eqs_numbers

# ##
# Evolution.thomas_algorithm!(oz)
# oz.solver_data.solver_corr

## evolution
function do_one_zone_burn!(oz::OneZone)
    # before loop actions
    StellarModels.create_output_files!(oz)
    StellarModels.evaluate_model_properties!(oz, oz.props)  # set the initial condition as the result of a previous phantom step
    retry_count = 0

    # evolution loop, be sure to have sensible termination conditions or this will go on forever!
    while true
        StellarModels.cycle_props!(oz)  # move props of previous step to prv_step_props of current step

        # step loop
        oz.solver_data.newton_iters = 0
        max_steps = oz.opt.solver.newton_max_iter
        if (oz.props.model_number == 0)
            max_steps = oz.opt.solver.newton_max_iter_first_step
        end

        exit_evolution = false
        retry_step = false
        oz.props = deepcopy(oz.prv_step_props)
        oz.props.dt = oz.prv_step_props.dt_next  # dt of this step becomes dt_next of previous
        oz.props.T = T_new(oz)
        oz.props.ρ = ρ_new(oz)

        corr = oz.solver_data.solver_corr
        equs = oz.solver_data.eqs_numbers

        # evaluate the equations for the first step
        Evolution.eval_jacobian_eqs!(oz)  # heavy lifting happens here!
        for i = 1:max_steps
            Evolution.thomas_algorithm!(oz)  # here as well

            (abs_max_corr, i_corr) = findmax(abs, corr)
            signed_max_corr = corr[i_corr]
            corr_nz = i_corr ÷ oz.nvars + 1
            corr_equ = i_corr % oz.nvars
            rel_corr = abs_max_corr / eps(oz.props.ind_vars[i_corr])

            # scale correction
            if oz.props.model_number == 0
                correction_multiplier = min(1.0, oz.opt.solver.initial_model_scale_max_correction / abs_max_corr)
            else
                correction_multiplier = min(1.0, oz.opt.solver.scale_max_correction / abs_max_corr)
            end
            if correction_multiplier < 1
                corr .*= correction_multiplier
            end

            # first try applying correction and see if it would give negative luminosity
            oz.props.ind_vars .+= corr
            oz.solver_data.newton_iters = i

            # evaluate the equations after correction and get residuals
            try
                StellarModels.evaluate_model_properties!(oz, oz.props)
                Evolution.eval_jacobian_eqs!(oz)  # heavy lifting happens here!

                (max_res, i_res) = findmax(abs, equs)
                res_nz = i_res ÷ oz.nvars + 1
                res_equ = i_res % oz.nvars

                #reporting
                if oz.opt.solver.report_solver_progress &&
                   i % oz.opt.solver.solver_progress_iter == 0
                    @show oz.props.model_number, i, rel_corr, signed_max_corr, corr_nz, corr_equ, max_res, res_nz,
                          res_equ
                end
                #check if tolerances are satisfied
                if rel_corr < oz.opt.solver.relative_correction_tolerance &&
                   max_res < oz.opt.solver.maximum_residual_tolerance
                    if oz.props.model_number == 0
                        println("Found first model")
                    end
                    break  # successful, break the step loop
                end
            catch e
                if isa(e, InterruptException)
                    throw(e)
                end
                println("Error while evaluating equations")
                showerror(stdout, e)
                retry_step = true
            end
            # if not, determine if we give up or retry
            if i == max_steps
                if retry_count > 10
                    exit_evolution = true
                    println("Too many retries, ending simulation")
                else
                    retry_count = retry_count + 1
                    retry_step = true
                    println("Failed to converge step $(oz.props.model_number) with timestep $(oz.props.dt/SECYEAR), retrying")
                end
            end
        end

        if retry_step
            StellarModels.uncycle_props!(oz)  # reset props to what prv_step_props contains, ie mimic state at end of previous step
            oz.props.dt_next *= oz.opt.timestep.dt_retry_decrease # adapt dt
            continue  # go back to top of evolution loop
        end

        if (exit_evolution)
            println("Terminating evolution")
            break
        end

        # step must be successful at this point
        retry_count = 0

        # increment age and model number since we accept the step.
        oz.props.time += oz.props.dt
        oz.props.model_number += 1

        # write state in oz.props and potential history/profiles.
        StellarModels.evaluate_model_properties!(oz, oz.props)
        StellarModels.write_data(oz)
        StellarModels.write_terminal_info(oz)

        if oz.opt.plotting.do_plotting && oz.props.model_number == 1
            Plotting.init_plots!(oz)
        elseif oz.opt.plotting.do_plotting && oz.props.model_number % oz.opt.plotting.plotting_interval == 0
            Plotting.update_plotting!(oz)
        end

        # check termination conditions
        if (oz.props.model_number >= oz.opt.termination.max_model_number)
            # StellarModels.write_terminal_info(oz; now=true)
            println("Reached maximum model number")
            break
        elseif oz.props.time >= oz.opt.termination.max_time * SECYEAR
            println("max time reached, ending simulation")
            break
        end

        # get dt for coming step
        oz.props.dt_next = Evolution.get_dt_next(oz)
    end
    if oz.opt.plotting.do_plotting
        Plotting.end_of_evolution(oz)
    end
    StellarModels.shut_down_IO!(oz)
end

open("example_options.toml", "w") do file
    write(file,
          """
          [solver]
          newton_max_iter_first_step = 1000
          initial_model_scale_max_correction = 0.2
          newton_max_iter = 50
          scale_max_correction = 0.1
          report_solver_progress = false
          solver_progress_iter = 1

          [timestep]
          dt_max_increase = 2.0
          delta_Xc_limit = 0.005
          max_dt = 1e-7

          [termination]
          max_model_number = 2000
          max_time = 3e-3

          [plotting]
          do_plotting = false
          wait_at_termination = false

          plotting_interval = 1

          window_specs = ["history"]
          window_layout = [[1, 1]]
          yaxes_log = [true]

          history_xaxis = "T"
          history_yaxes = ["n", "H1", "D2", "He4", "Li7"]

          [io]
          profile_interval = 50
          terminal_header_interval = 100
          terminal_info_interval = 10
          history_values = ["age", "dt", "model_number", "T", "ρ",
                            "n", "H1", "D2", "T3", "He3", "He4", "Li7", "Be7"]

          """)
end
StellarModels.set_options!(oz.opt, "./example_options.toml")
@time do_one_zone_burn!(oz)


##
using CairoMakie, LaTeXStrings, MathTeXEngine
basic_theme = Theme(fonts=(regular=texfont(:text), bold=texfont(:bold),
                           italic=texfont(:italic), bold_italic=texfont(:bolditalic)),
                    fontsize=30, size=(1000, 750), linewidth=7,
                    Axis=(xlabelsize=40, ylabelsize=40, titlesize=40, xgridvisible=false, ygridvisible=false,
                          spinewidth=2.5, xminorticksvisible=true, yminorticksvisible=true, xtickalign=1, ytickalign=1,
                          xminortickalign=1, yminortickalign=1, xticksize=14, xtickwidth=2.5, yticksize=14,
                          ytickwidth=2.5, xminorticksize=7, xminortickwidth=2.5, yminorticksize=7, yminortickwidth=2.5,
                          xticklabelsize=35, yticklabelsize=35, xticksmirrored=true, yticksmirrored=true),
                    Legend=(patchsize=(70, 10), framevisible=false, patchlabelgap=20, rowgap=10))
set_theme!(basic_theme)
# GLMakie.activate!()
##
### Plot the history
f = Figure();
ax = Axis(f[1, 1]; xlabel=L"k_\textrm{B}T\,\textrm{(MeV)}", ylabel=L"\log_{10}(X)", xreversed=true, xscale=log10, yticksmirrored=false)
history = StellarModels.get_history_dataframe_from_hdf5("history.hdf5");
temperature = history[!, "T"] .* KERG ./ Constants.MEV_TO_ERGS;
lines!(ax, temperature, log10.(history[!, "n"]), label=L"n")
lines!(ax, temperature, log10.(history[!, "H1"]), label=L"^1H")
lines!(ax, temperature, log10.(history[!, "D2"]), label=L"^2H")
lines!(ax, temperature, log10.(history[!, "He3"]), label=L"^3He")
lines!(ax, temperature, log10.(history[!, "He4"]), label=L"^4He")
lines!(ax, temperature, log10.(history[!, "Li7"]), label=L"^7Li")
lines!(ax, temperature, log10.(history[!, "Be7"]), label=L"^7Be")
axislegend(position=:lc)
# hidespines!(ax, :r)

# ax2 = Axis(f[1, 1]; yaxisposition=:right, ylabel=L"\log_{10}(T, ρ)",xreversed = true, xscale=log10, yticksmirrored=false)
# hidespines!(ax2, :l, :t, :b)
# hidexdecorations!(ax2, ticks=false)
# lines!(ax2, temperature, log10.(history[!, "T"]), color=:red, label=L"T")
# lines!(ax2, temperature, log10.(history[!, "ρ"]), color=:blue, label=L"ρ")
# axislegend(position=:rb)
ylims!(ax, -11, 1)
f

##  number fractions
println("BBN results")
println("He4 ", history[end, "He4"])
println("number fractions")
for el in ["D2", "He3", "Li7"]
    massfraction = history[end, el]
    println(el, " ", massfraction / (Chem.isotope_list[Symbol(el)].mass * history[end, "H1"]))
end