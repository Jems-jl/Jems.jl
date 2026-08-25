#=
# DualEvol.jl

This notebook provides a simple example of evolving a star with dual numbers. The dual numbers are used to compute the
derivative of the luminosity with respect to the mass of the star.
=#
using BenchmarkTools
using ForwardDiff
using Jems.Chem
using Jems.Constants
using Jems.EOS
using Jems.Opacity
using Jems.NuclearNetworks
using Jems.Turbulence
using Jems.StellarModels
using Jems.Evolution
using Jems.Plotting
using Jems.DualSupport



##
# set up the basic properties of the model
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16], [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
nz = 1000
nextra = 100
eos = EOS.IdealEOS(false)
opacity = Opacity.SimpleElectronScatteringOpacity()
turbulence = Turbulence.BasicMLT(1.0)

##
#= Now we define the dual numbers we'll be using.
We need to use two different types of duals, one for derivatives with respect to initial parameters, and one that JEMS
uses internally. We need to specify the internal dual type here expicitly, because we need to register both dual types
with the ForwardDiff package in order to for it to know the ordering of derivatives.
Finally, we the external type is:
    ForwardDiff.Dual{external_tag, Float64, 2}
where the tag is a symbol that we register with ForwardDiff, and the number of partials is 2, one for the initial mass,
and one of the initial hydrogen fraction X.
=#
external_tag = ForwardDiff.Tag{:external, nothing}
ForwardDiff.tagcount(external_tag)  # register tag
internal_tag = ForwardDiff.Tag{:internal, nothing}
ForwardDiff.tagcount(internal_tag)
external_number_type = ForwardDiff.Dual{external_tag, Float64, 2}

sm = StellarModel(StellarModels.DefaultStellarEquationSet(), nz, nextra, net, eos, opacity, turbulence; number_type=external_number_type, internal_dual_tag=internal_tag);

## setup the initial condition. Since we want derivatives with respect to the mass, we need to use a dual number for
# the mass. Its partial derivative, ∂logM/∂logM, is trivially 1.0, so we set the dual number to be 1.0 * MSUN + 1.0 * ε.
n = 3
logM_dual = ForwardDiff.Dual{external_tag}(log(1.0), 1.0, 0.0)  # we use the derivative to log(M) in units of Msun
mass_dual = 10.0^logM_dual*MSUN
X_dual = ForwardDiff.Dual{external_tag}(0.7154, 0.0, 1.0)  # we use the derivative to X, the initial hydrogen fraction
StellarModels.n_polytrope_initial_condition!(n, sm, nz, X_dual, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], mass_dual,
                                             100 * RSUN; initial_dt=10 * SECYEAR)
Evolution.compute_starting_model_properties!(sm)  # populate the stellar model with the initial condition

## let's see what the numbers of the stellar model look like.

#= 
The luminosity at the surface is of type LocalDualData, which is a container with parametric type
    LocalDualData{N + 1, 3N + 1, TNUMBER, DUAL_TAG}
where N is the number of independent variables, TNUMBER is the type of the external number of the calculations, and 
DUAL_TAG is the tag of the internal dual numbers used by JEMS.
=#
print(typeof(sm.props.L[sm.props.nz]))

#= getting its value is done with the `get_value` function, which returns a number of type TNUMBER, in this case
a dual number of the external type, carrying the derivative to initial mass. Note that even the time is a dual number
of the external type.
=#
print(typeof(get_value(sm.props.L[sm.props.nz])))
print(get_value(sm.props.lnρ[1]))
print(sm.props.time[1])

## define some setting for actual evolution
open("example_options.toml", "w") do file
    write(file,
          """
          [remesh]
          do_remesh = true

          [solver]
          newton_max_iter_first_step = 1000
          initial_model_scale_max_correction = 0.2
          newton_max_iter = 10
          scale_max_correction = 0.1

          [timestep]
          dt_max_increase = 1.5
          delta_R_limit = 0.01
          delta_Tc_limit = 0.005
          delta_Xc_limit = 0.001

          [termination]
          max_model_number = 5000
          min_center_X = 0.01

          [io]
          profile_interval = 50
          terminal_header_interval = 100
          terminal_info_interval = 100

          """)
end
StellarModels.set_options!(sm.opt, "./example_options.toml")
rm(sm.opt.io.hdf5_history_filename; force=true)
rm(sm.opt.io.hdf5_profile_filename; force=true)


## run model
@time Evolution.do_evolution_loop!(sm, plotter=Plotting.NullPlotter());

##
using GLMakie, LaTeXStrings
set_theme!(Plotting.basic_theme())
##
history = StellarModels.get_history_dataframe_from_hdf5("history.hdf5");
M_derivatives = StellarModels.get_ith_partial_history_dataframe_from_hdf5("history.hdf5", 1);
X_derivatives = StellarModels.get_ith_partial_history_dataframe_from_hdf5("history.hdf5", 2);
get_first_partial(x) = x.partials[1]
get_second_partial(x) = x.partials[2]
get_real_value(x) = x.value
## Let's see the HRD
# It's a very luminous 1 Msun star, but that is due to the opacities used (or rather the lack thereof).
f = Figure();
ax = Axis(f[1, 1]; xlabel=L"\log_{10}(T_\text{eff}/[K])", ylabel=L"\log_{10}(L/L_\odot)", xreversed=true)
lines!(ax, log10.(history[!, "T_surf"]), log10.(history[!, "L_surf"]))
f

##
#=
Let's look at the derivative of the luminosity wrt the mass.
We see its 3 for most of the evolution, as expected, but we also see it drops significantly just as the star
descends to settle on the main sequence.
Remember that this derivative is computed at constant model number, not at any constant physical parameter. Meaning that
if all else stays the same, e.g. at around model number 530, the derivative drops to about 1, meaning that if everything
else was the same, even the different timesteps that go from model 1 to model 530, the luminosity at model 530 of a 
slightly more massive star would scale linearly with mass.
In phases where the star is changing luminosity, the partial derivative is therefore not very meaningful, since
    dL/dm = ∂L/∂m - dL/dX * ∂X/∂m
where X is any other parameter that is changing with time and input mass. If L is not changing with time, the slope
dL/dX is zero, and the partial derivative is equal to the total derivative.
This also shows why the derivative on the main sequence is distorted: 3.1 and climbing to 3.5 vs the expected 3.0.
=#
f = Figure();
ax = Axis(f[1, 1]; ylabel=L"L(t)",
          xticklabelsvisible=false)
lines!(ax, history[!, "model_number"], log10.(history[!, "L_surf"]))
vlines!(ax, 528, color=:red, linestyle=:dash)
ax2 = Axis(f[2, 1]; xlabel="model number", ylabel=L"\frac{\partial \ln(L/L_\odot)}{\partial \ln(M/M_\odot)}")
lines!(ax2, history[!, "model_number"], log10(exp(1)) ./ history[!, "L_surf"] .* M_derivatives[!, "L_surf"])
vlines!(ax2, 528, color=:red, linestyle=:dash)
linkxaxes!(ax, ax2)
display(GLMakie.Screen(), f)


##
#=
Instead, we can compute the derivative at constant surface temprature and at constant Xᵪ, which are more physically 
meaningful as the star contracts and arrives on the main sequence, respectively.
The widely known mass-luminosity relation: L ∝ M^3 μ^4, is valid for stars that are fully radiative with constant
opacity, and the EOS being the ideal gas.
=#

#=
First, we reconstruct the dual numbers from the outputs.
=#
Lduals = ForwardDiff.Dual{}.(history[!, "L_surf"], M_derivatives[!, "L_surf"], X_derivatives[!, "L_surf"]);
X_c_duals = ForwardDiff.Dual{}.(history[!, "X_center"], M_derivatives[!, "X_center"], X_derivatives[!, "X_center"]);
T_duals = ForwardDiff.Dual{}.(history[!, "T_surf"], M_derivatives[!, "T_surf"], X_derivatives[!, "T_surf"]);

#=
Then, we interpolate the luminosity with respect to the central hydrogen and the temperature, which acts like performing
the implicit function theorem to get a total derivative at constant X_c and T, respectively.
=#
ts_for_X_c = []
Ldual_int_at_X_c = []
X_c_eval = []
for i in eachindex(Lduals)[1:end-1]
    if abs(X_c_duals[i+1].value - X_c_duals[i].value) < 1e-10
        continue
    end
    # interpolate Ldual to a range of X_c values in between X_c_duals[i] and X_c_duals[i+1]
    X_c_range = range(X_c_duals[i].value, X_c_duals[i+1].value, length=10)
    X_c_eval = vcat(X_c_eval, X_c_range)
    for X_c in X_c_range
        Ldual_interp = interpolate_to_value(X_c_duals[i], X_c_duals[i+1], Lduals[i], Lduals[i+1], X_c)
        t_interp = interpolate_to_value(X_c_duals[i].value, X_c_duals[i+1].value, history[!, "age"][i], history[!, "age"][i+1], X_c)
        push!(ts_for_X_c, t_interp)
        push!(Ldual_int_at_X_c, Ldual_interp)
    end
end
ts_for_T = []
Ldual_int_at_T = []
T_eval = []
for i in eachindex(Lduals)[1:end-1]
    if abs(T_duals[i+1].value - T_duals[i].value) < 1e-10
        continue
    end
    # interpolate Ldual to a range of T values in between T_duals[i] and T_duals[i+1]
    T_range = range(T_duals[i].value, T_duals[i+1].value, length=10)
    T_eval = vcat(T_eval, T_range)
    for T in T_range
        Ldual_interp = interpolate_to_value(T_duals[i], T_duals[i+1], Lduals[i], Lduals[i+1], T)
        t_interp = interpolate_to_value(T_duals[i].value, T_duals[i+1].value, history[!, "age"][i], history[!, "age"][i+1], T)
        push!(ts_for_T, t_interp)
        push!(Ldual_int_at_T, Ldual_interp)
    end
end


##

#=
We see that with the corrections to X_c we are a lot closer the expected value of 3.0 during the bulk of the main
sequence. Toward the end, the rise is due to the change in opacity (which is a function of molecular weight).
Also, at ZAMS, the zigzag is a result of the slope dL/dX_c being very large (L changes a lot when settling on the main
sequencewhile X_c still changes slowly).

In the panel with the derivative at constant T, we see that for the bulk of the pre-main sequence contraction (from 1e5
to 1e6 years), the (total) derivative is still close to 3.0. On the main sequence, this derivative is not very
meaningful, because T is mostly contanst, making this derivative very sensitive. Also, a different mass star will have
a much different temprature during its main sequence evolution.
=#
f = Figure();
ax1 = Axis(f[1, 1]; xscale=log10, xticklabelsvisible=false,
          xminorticksvisible=true,
          xminorticks=IntervalsBetween(10),
          xticks=(10 .^ (2:1:9)), ylabel=L"\frac{\partial \ln(L/L_\odot)}{\partial \ln(M/M_\odot)}|_{X_c}",)
scatter!(ax1, history[!, "age"], 1 ./ (log(10) * history[!, "L_surf"]) .* M_derivatives[!, "L_surf"])
lines!(ax1, history[!, "age"][1:(end - 1)],
       1 ./ (log(10) * history[!, "L_surf"][1:(end - 1)]) .*
    construct_derivative.(M_derivatives[!, "L_surf"][1:(end - 1)], M_derivatives[!, "X_center"][1:(end - 1)],
    get_slope.(Ref(history[!, "X_center"]), Ref(history[!, "L_surf"]), 1:(length(history[!, "X_center"]) - 1))), color=:orange)
scatter!(ax1, ts_for_X_c, 1 ./ (log(10) * get_real_value.(Ldual_int_at_X_c)) .* get_first_partial.(Ldual_int_at_X_c), color=:red)
ylims!(ax1, 2, 4.0)

ax2 = Axis(f[2, 1]; xlabel="Age (yr)", ylabel=L"\frac{\partial \ln(L/L_\odot)}{\partial \ln(M/M_\odot)}|_{T}",
           xscale=log10,
           xminorticksvisible=true,
           xminorticks=IntervalsBetween(10),
           xticks=(10 .^ (2:1:9)))
scatter!(ax2, history[!, "age"], 1 ./ (log(10) * history[!, "L_surf"]) .* M_derivatives[!, "L_surf"])
lines!(ax2, history[!, "age"][1:(end - 1)],
    1 ./ (log(10) *  history[!, "L_surf"][1:(end - 1)]) .*
    construct_derivative.(M_derivatives[!, "L_surf"][1:(end - 1)], M_derivatives[!, "T_surf"][1:(end - 1)],
    get_slope.(Ref(history[!, "T_surf"]), Ref(history[!, "L_surf"]), 1:(length(history[!, "T_surf"]) - 1))), color=:orange)
scatter!(ax2, ts_for_T, 1 ./ (log(10) * get_real_value.(Ldual_int_at_T)) .* get_first_partial.(Ldual_int_at_T),
         color=:red)
linkxaxes!(ax1, ax2)
ylims!(ax2, 2, 4.0)
display(GLMakie.Screen(), f)

##
#=
Next, we investigate the derivative of the luminosity with respect to the initial hydrogen fraction X. We see that it
is negative, as expected, since a higher initial hydrogen fraction means a higher opacity. Our opacity law:
    κ = 0.2 * (1 + X) cm²/g

=#
f = Figure();
ax1 = Axis(f[1, 1]; xscale=log10, xticklabelsvisible=false,
           xminorticksvisible=true,
           xminorticks=IntervalsBetween(10),
           xticks=(10 .^ (2:1:9)), ylabel=L"\frac{\partial \log L}{\partial X_{\mathrm{init}}}|_{X_c}",)
scatter!(ax1, history[!, "age"], 1 ./ (log(10) * history[!, "L_surf"]) .* X_derivatives[!, "L_surf"])
lines!(ax1, history[!, "age"][1:(end - 1)],
       1 ./ (log(10) * history[!, "L_surf"][1:(end - 1)]) .*
       construct_derivative.(X_derivatives[!, "L_surf"][1:(end - 1)], X_derivatives[!, "X_center"][1:(end - 1)],
       get_slope.(Ref(history[!, "X_center"]), Ref(history[!, "L_surf"]), 1:(length(history[!, "X_center"]) - 1))), color=:orange)
scatter!(ax1, ts_for_X_c, 1 ./ (log(10) * get_real_value.(Ldual_int_at_X_c)) .* get_second_partial.(Ldual_int_at_X_c), color=:red)
ylims!(ax1, -2, -0.5)
ax2 = Axis(f[2, 1]; xlabel="Age (yr)", ylabel=L"\frac{\partial \log L}{\partial X_{\mathrm{init}}}|_{T}",
           xscale=log10,
           xminorticksvisible=true,
           xminorticks=IntervalsBetween(10),
           xticks=(10 .^ (2:1:9)))
scatter!(ax2, history[!, "age"], 1 ./ (log(10) * history[!, "L_surf"]) .* X_derivatives[!, "L_surf"])
lines!(ax2, history[!, "age"][1:(end - 1)],
       1 ./ (log(10) * history[!, "L_surf"][1:(end - 1)]) .*
       construct_derivative.(X_derivatives[!, "L_surf"][1:(end - 1)], X_derivatives[!, "T_surf"][1:(end - 1)],
       get_slope.(Ref(history[!, "T_surf"]), Ref(history[!, "L_surf"]), 1:(length(history[!, "T_surf"]) - 1))), color=:orange)
scatter!(ax2, ts_for_T, 1 ./ (log(10) * get_real_value.(Ldual_int_at_T)) .* get_second_partial.(Ldual_int_at_T),
         color=:red)
linkxaxes!(ax1, ax2)
ylims!(ax2, -4, -0.5)
f