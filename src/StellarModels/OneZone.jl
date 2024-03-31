
using ForwardDiff

@kwdef mutable struct OneZone{TNUMBER<:Real,TDUALFULL<:ForwardDiff.Dual,
                              TPROPS<:StellarModels.AbstractModelProperties,
                              TNET<:NuclearNetworks.AbstractNuclearNetwork,
                              TSOLVER<:StellarModels.AbstractSolverData} <: AbstractModel
    nvars::Int  # This is the sum of hydro vars and species
    var_names::Vector{Symbol}  # List of variable names
    vari::Dict{Symbol,Int}  # Maps variable names to ind_vars vector

    ## Properties related to the solver ##
    # List of equations to be solved.
    composition_equation::StellarModels.TypeStableEquation{OneZone{TNUMBER,TDUALFULL,TPROPS,TNET,TSOLVER},TDUALFULL}
    # cache to store residuals and solver matrices
    solver_data::TSOLVER

    # Microphyical models
    network::TNET

    # Properties that define the model
    prv_step_props::TPROPS  # properties of the previous step
    props::TPROPS  # properties during and after newton solving

    # Space for used defined options, defaults are in Options.jl
    opt::StellarModels.Options

    # object holding plotting things, ie figures, data to plot.
    plt::StellarModels.Plotter

    # Output files
    history_file::HDF5.File
end

"""
    OneZone(varnames::Vector{Symbol}, composition_equation::Function,
                nvars::Int, nspecies::Int)

Constructor for a `OneZone` instance, using `varnames` for the independent variables, the composition equation
to be solved, number of independent variables `nvars`, number of species in the network `nspecies`
"""
function OneZone(compostion_equation::Function, network::NuclearNetwork, use_static_arrays=true, number_type=Float64)
    nvars = network.nspecies
    var_names_full = network.species_names
    # link var_names to the correct index so you can do ind_var[vari[:lnT]] = 'some temperature'
    vari::Dict{Symbol,Int} = Dict()
    for i in eachindex(network.species_names)
        vari[var_names_full[i]] = i
    end

    solver_data = StellarModels.SolverData(nvars, 1, 0, use_static_arrays, number_type)

    # properties
    prv_step_props = OneZoneProperties(nvars, length(network.reactions), network.nspecies, number_type)
    props = OneZoneProperties(nvars, length(network.reactions), network.nspecies, number_type)

    # create type stable function objects
    dual_sample = ForwardDiff.Dual(zero(number_type), (zeros(number_type, nvars)...))
    tpe_stbl_func = StellarModels.TypeStableEquation{OneZone{number_type,typeof(dual_sample),typeof(props),
                                                             typeof(network),typeof(solver_data)},
                                                     typeof(dual_sample)}(compostion_equation)

    opt = StellarModels.Options()  # create options object
    plt = StellarModels.Plotter()

    # create the stellar model
    oz = OneZone(; nvars=nvars,
                 var_names=var_names_full, vari=vari,
                 solver_data=solver_data,
                 composition_equation=tpe_stbl_func, network=network,
                 prv_step_props=prv_step_props, props=props,
                 opt=opt, plt=plt,
                 history_file=HDF5.File(-1, ""))
    return oz
end

@kwdef mutable struct OneZoneProperties{TN,TDual,TCellDualData} <: AbstractModelProperties
    # scalar quantities
    dt::TN  # Timestep of the current evolutionary step (s)
    dt_next::TN
    time::TN  # Age of the model (s)
    model_number::Int
    nz = 1  # duh; however  solver needs to know what nz is

    # T and ρ for the zone
    T::TN
    ρ::TN

    # array of the values of the independent variables, everything should be reconstructable from this (mesh dependent)
    ind_vars::Vector{TN}

    # independent variables (duals constructed from the ind_vars array)
    xa::Vector{TCellDualData}   # dim-less
    xa_dual::Vector{TDual}      # only the cell duals wrt itself

    # rates
    rates::Vector{TCellDualData}  # g^-1 s^-1
    rates_dual::Vector{TDual}     # only cell duals wrt itself

    ϵ_nuc::TN
end

function OneZoneProperties(nvars::Int, nrates::Int, nspecies::Int, ::Type{TN}) where {TN<:Real}
    # define the types
    CDDTYPE = CellDualData{nvars + 1,3 * nvars + 1,TN}  # full dual arrays
    TD = typeof(ForwardDiff.Dual(zero(TN), (zeros(TN, nvars))...))  # only the cell duals

    # create the vector containing the independent variables
    ind_vars = zeros(TN, nvars)

    xa_dual = zeros(TD, nvars)
    xa = Vector{CDDTYPE}(undef, nspecies)
    for j = 1:nspecies
        xa[j] = CellDualData(nvars, TN; is_ind_var=true, ind_var_i=nvars - nspecies + j)
    end
    rates_dual = zeros(TD, nrates)
    rates = Vector{CDDTYPE}(undef, nrates)
    for k = 1:nrates
        rates[k] = CellDualData(nvars, TN)
    end

    return OneZoneProperties(; ind_vars=ind_vars, model_number=zero(Int), dt=zero(TN), dt_next=zero(TN), time=zero(TN),
                             xa=xa, xa_dual=xa_dual, rates=rates, rates_dual=rates_dual, T=zero(TN), ρ=zero(TN),
                             ϵ_nuc=zero(TN))
end

"""
    function evaluate_stellar_model_properties!(oz, props::StellarModelProperties{TDual, TCellDualData}) where
        {TDual <: ForwardDiff.Dual, TCellDualData}

Evaluates the stellar model properties `props` from the `ind_vars` array. The goal is to save the 'state' of the
StellarModel so we can easily get properties like rates, eos, opacity values, and retrace if a retry is called.
This does _not_ update the mesh/ind_vars arrays.
"""
function evaluate_one_zone_properties!(oz, props::OneZoneProperties{TN,TDual}) where {TN<:Real,TDual<:ForwardDiff.Dual}
    # update independent variables
    for j = 1:(oz.network.nspecies)
        update_cell_dual_data_value!(props.xa[j], props.ind_vars[oz.nvars - oz.network.nspecies + j])
        props.xa_dual[j] = get_cell_dual(props.xa[j])
    end

    # evaluate rates
    set_rates_for_network!(props.rates_dual, oz.network, props.T, props.ρ, props.xa_dual)
    for j in eachindex(props.rates_dual)
        update_cell_dual_data!(props.rates[j], props.rates_dual[j])
    end

    props.ϵ_nuc = 0.0
    for j in eachindex(props.rates_dual)
        props.ϵ_nuc += props.rates_dual[j].value * oz.network.reactions[j].Qvalue
    end
end

function cycle_props!(oz::OneZone)
    temp_props = oz.prv_step_props
    oz.prv_step_props = oz.props
    oz.props = temp_props
end

function uncycle_props!(oz::OneZone)
    temp_props = oz.props
    oz.props = oz.prv_step_props
    oz.prv_step_props = temp_props
end

function get_dt_next(oz::OneZone)
    dt_next = oz.props.dt  # this it calculated at end of step, so props.dt is the dt we used to do this step

    X = get_value(oz.props.xa[1, oz.network.xa_index[:H1]])
    Xold = get_value(oz.prv_step_props.xa[1, oz.network.xa_index[:H1]])
    ΔX = abs(X - Xold)

    dt_nextX = dt_next * oz.opt.timestep.delta_Xc_limit / ΔX

    min_dt = dt_next * oz.opt.timestep.dt_max_decrease
    dt_next = min(oz.opt.timestep.dt_max_increase * dt_next, dt_nextX)
    dt_next = max(dt_next, min_dt)
    println("new dt", dt_next)
    return dt_next
end

function do_one_zone_burn!(oz::OneZone)
    # before loop actions
    # StellarModels.create_output_files!(oz)
    evaluate_one_zone_properties!(oz, oz.props)  # set the initial condition as the result of a previous phantom step
    retry_count = 0

    # evolution loop, be sure to have sensible termination conditions or this will go on forever!
    while true
        cycle_props!(oz)  # move props of previous step to prv_step_props of current step

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

        corr = oz.solver_data.solver_corr
        equs = oz.solver_data.eqs_numbers

        # evaluate the equations for the first step
        eval_jacobian_eqs!(oz)  # heavy lifting happens here!
        for i = 1:max_steps
            thomas_algorithm!(oz)  # here as well

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
                evaluate_one_zone_properties!(oz, oz.props)
                eval_jacobian_eqs!(oz)  # heavy lifting happens here!

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
            uncycle_props!(oz)  # reset props to what prv_step_props contains, ie mimic state at end of previous step
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
        evaluate_one_zone_properties!(oz, oz.props)
        StellarModels.write_data(oz)
        StellarModels.write_terminal_info(oz)

        # if oz.opt.plotting.do_plotting && oz.props.model_number == 1
        #     Plotting.init_plots!(oz)
        # elseif oz.opt.plotting.do_plotting && oz.props.model_number % oz.opt.plotting.plotting_interval == 0
        #     Plotting.update_plotting!(oz)
        # end

        # check termination conditions
        if (oz.props.model_number > oz.opt.termination.max_model_number)
            # StellarModels.write_terminal_info(oz; now=true)
            println("Reached maximum model number")
            break
        end

        # get dt for coming step
        oz.props.dt_next = get_dt_next(oz)
    end
    # if oz.opt.plotting.do_plotting
    #     Plotting.end_of_evolution(oz)
    # end
    StellarModels.close_output_files!(oz)
end
