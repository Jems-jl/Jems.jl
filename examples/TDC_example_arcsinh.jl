using BenchmarkTools 
using Revise
using Jems.Chem
using Jems.Constants
using Jems.EOS
using Jems.Opacity
using Jems.NuclearNetworks
using Jems.Turbulence
using Jems.StellarModels 
using Jems.Evolution
using Jems.Plotting
using Jems.Interpolations
include("TDC_StellarModels2.jl")
include("TDC_InitialCondition.jl")
##
# Get MESA data if not available
if !isdir("MESA_data")
    mkdir("MESA_data")
    download("https://zenodo.org/records/19722306/files/mesa-26.04.1.zip", "MESA_data/mesa-26.04.1.zip")
    
    cd("MESA_data")
    run(`unzip -q mesa-26.04.1.zip`)
    
    mv("mesa-26.04.1/kap/kap_data.tar.xz", "./kap_data.tar.xz")
    run(`tar -xJf kap_data.tar.xz`)
    rm("kap_data.tar.xz")

    mkdir("eosFreeEOS_data_filtered")
    mv("mesa-26.04.1/eos/eosFreeEOS_data.tar.xz", "./eos_data.tar.xz")
    run(`tar -xJf eos_data.tar.xz`)
    rm("eos_data.tar.xz")

    extracted_eos_dir = "eosFreeEOS_data"
    target_pattern = r"^mesa-FreeEOS_\d+z\d+x\.data$"
    
    for file in readdir(extracted_eos_dir)
        if occursin(target_pattern, file)
            mv(joinpath(extracted_eos_dir, file), joinpath("eosFreeEOS_data_filtered", file))
        end
    end

    rm("mesa-26.04.1", recursive=true)
    rm("eosFreeEOS_data", recursive=true)
    mv("eosFreeEOS_data_filtered", "eosFreeEOS_data")
    rm("mesa-26.04.1.zip")
    
    cd("../")
end
##
EXCESS_FACTOR[] = 1e-8
α_w_FACTOR[] = 0.25
α_Λ_FACTOR[] = 1.0

##
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16], [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
# net = merge_nuclear_networks([NuclearNetworks.networks[:JINA_PPI],
#                             NuclearNetworks.networks[:JINA_PPII],
#                             NuclearNetworks.networks[:JINA_CNOI],
#                             NuclearNetworks.networks[:JINA_CNOII],])
nz = 1000
nextra = 100
#eos = EOS.IdealEOS(true)
eos_table = EOSTableCollector("MESA_data/eosFreeEOS_data", include_radiation = true)
#println("Memory used: $(Base.summarysize(eos_table) / 1024^3) GB")
low_T_collection = OpacityTableCollector("MESA_data/kap_data", "lowT_fa05_gs98") 
high_T_collection = OpacityTableCollector("MESA_data/kap_data","oplib_agss09" ) 
opacity = CompositeOpacity(low_T_collection, high_T_collection, 3.8, 4.2)
# eos = EOS.IdealEOS(true)
# opacity = Opacity.SimpleElectronScatteringOpacity()
turbulence = Turbulence.BasicMLT(1.0)
##
nz = 1000
nextra = 10000
sm = StellarModel(TDCEquationSet(), nz, nextra, net, eos_table, opacity, turbulence);
##
nz = sm_alt.props.nz
nextra = 10000
sm = StellarModel(TDCEquationSet(), nz, nextra, net, eos_table, opacity, turbulence);

##
sm.props.m[1:nz] = sm_alt.props.m[1:nz]
sm.props.dm[1:nz] = sm_alt.props.dm[1:nz]
for i in 1:nz
    for k in 1:4
        sm.props.ind_vars[(i-1)*sm.nvars + k] = sm_alt.props.ind_vars[(i-1)*sm_alt.nvars + k]
    end
    if i == sm.props.nz
        v_turb = get_value(sm_alt.props.turb_res[i-1].v_turb)
    else
        v_turb = get_value(sm_alt.props.turb_res[i].v_turb)
    end
    v_turb_floored = max(1e-5, v_turb)
    sm.props.ind_vars[(i-1)*sm.nvars + 5] = max(1e-5, asinh(abs((v_turb_floored) / sqrt(2.0))))
    for k in 1:sm.network.nspecies
        sm.props.ind_vars[(i-1)*sm.nvars + 5 + k] = sm_alt.props.ind_vars[(i-1)*sm_alt.nvars + 4 + k]
    end
end
Evolution.compute_starting_model_properties!(sm)
sm.props.dt = sm_alt.props.dt

##

n = 1.5
# Polytropic Initial condition
tdc_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], 
                                            1 * MSUN, 80 * RSUN; initial_dt=10 * SECYEAR)
Evolution.compute_starting_model_properties!(sm)

# Adiabatic Initial condition
##
using CairoMakie
radius = Float64[]
input_radii = Float64[]
density = Float64[]
temp = Float64[]
for i in 1:10:300
n = 1.5
tdc_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], 
                                            1 * MSUN, i * RSUN; initial_dt=0.1 * SECYEAR)
Evolution.compute_starting_model_properties!(sm)
r = exp(get_value(sm.props.lnr[sm.props.nz]))/RSUN
den = exp(get_value(sm.props.lnρ[1]))
t = exp(get_value(sm.props.lnT[1]))
push!(radius,r)
push!(input_radii, Float64(i))
push!(temp, t)
push!(density, den)
end 
##
f= Figure(resolution = (1200, 800));
ax1 = Axis(f[1,1];xlabel = L"input radius", ylabel = "central temp", title = "variation")
scatter!(ax1, input_radii, log10.(temp), color = :blue, markersize = 10)
f
#
##
@benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
end

##
#=
And next we benchmark the evaluation of the model equations and construction of the Jacobian:
=#
@benchmark begin
    Evolution.eval_jacobian_eqs!($sm)
end
##


##
#=

To benchmark the linear solver itself we need to perform
the jacobian evaluation as a setup for the benchmark. This is because the solver
destroys the Jacobian to perform in-place operations.
=#

@benchmark begin
    Evolution.block_tridiagonal_solver!($sm, $sm.solver_data)
end setup=(Evolution.eval_jacobian_eqs!($sm))
##

function calculate_nabla_tdc(sm:: StellarModel, k :: Int)
    L = get_00_dual(sm.props.L[k]) * LSUN
    γ₀ = abs(get_00_dual(sm.props.gamma_turb[k]))
    ω = (sinh(γ₀))^2
    m₀ = sm.props.m[k]
    r₀ = exp(get_00_dual(sm.props.lnr[k]))
    if k == sm.props.nz
        P = get_00_dual(sm.props.eos_res[k].P)
        ρ = get_00_dual(sm.props.eos_res[k].ρ)
        T = get_00_dual(sm.props.eos_res[k].T)
        κ = get_00_dual(sm.props.κ[k])
        cₚ = get_00_dual(sm.props.eos_res[k].cₚ)
        Hₚ = P / (ρ * CGRAV * m₀ / r₀^2)
        Λ =  α_Λ_FACTOR[]*Hₚ #1/(1/ (α_Λ_FACTOR[]* Hₚ) + 1/r₀)
        k_rad = 16 * SIGMA_SB * T^3 / (3 * κ * ρ)
        α₂ = ρ*cₚ*0.5*sqrt(2/3)*Λ*sqrt(ω)
        ∇ᵣ = 3 * κ * L * P / (16π * CRAD * CLIGHT * CGRAV * m₀ * T^4)
        ∇ₐ = get_00_dual(sm.props.eos_res[k].∇ₐ)
        SA = (∇ᵣ - ∇ₐ)*(1 + α₂/k_rad)^(-1)        
        ∇ = ∇ₐ + SA
        return ∇
    end 

    P_face = exp(get_00_dual(sm.props.lnP_face[k]))
    ρ_face = exp(get_00_dual(sm.props.lnρ_face[k]))
    T_face = exp(get_00_dual(sm.props.lnT_face[k]))
    cₚ =  get_00_dual(sm.props.cₚ_face[k])
    κ = get_00_dual(sm.props.κ_face[k])
    Hₚ = P_face / (ρ_face * CGRAV * m₀ / r₀^2) 
    Λ = 1/(1/( α_Λ_FACTOR[]*Hₚ) + 1/r₀)
    m₀ = sm.props.m[k]
    k_rad = 16 * SIGMA_SB * T_face^3 / (3 * κ * ρ_face)
    α₂ = ρ_face*cₚ*0.5*sqrt(2/3)*Λ*sqrt(ω)
    ∇ᵣ = 3 * κ * L * P_face / (16π * CRAD * CLIGHT * CGRAV * m₀ * T_face^4)
    ∇ₐ = get_00_dual(sm.props.∇ₐ_face[k])
    SA = (∇ᵣ - ∇ₐ)*(1 + α₂/k_rad)^(-1)
    ∇ = ∇ₐ + SA 
    return ∇
end 
function get_nabla_tdc(sm, k)
    return calculate_nabla_tdc(sm, k).value
end

function get_D_turb(sm::StellarModel, k:: Int)
    r = exp(get_00_dual(sm.props.lnr[k]))
    P = exp(get_00_dual(sm.props.lnP_face[k]))
    ρ = exp(get_00_dual(sm.props.lnρ_face[k]))
    γ_face_00 = abs(get_00_dual(sm.props.gamma_turb[k]))
    ω = (sinh(γ_face_00))^2
    m = sm.props.m[k]
    g = CGRAV * m / (r^2)
    Hp = P / ( ρ * g)
    Λ = 1/(1/(Hp) + 1/r)
    D = 1/3 * Λ * sqrt(2 * ω)
    return ForwardDiff.value(D)
end 

function get_lambda_turb(sm::StellarModel, k:: Int)
    r = exp(get_00_dual(sm.props.lnr[k]))
    P = exp(get_00_dual(sm.props.lnP_face[k]))
    ρ = exp(get_00_dual(sm.props.lnρ_face[k]))
    m = sm.props.m[k]
    g = CGRAV * m / (r^2)
    Hp = P / ( ρ * g)
    Λ = 1/(1/(Hp) + 1/r)
    return ForwardDiff.value(Λ)
end 

function interpolate_boundary_live(sm::StellarModel, i_rad::Int)
    i_conv = i_rad - 1
    logr_rad  = get_value(sm.props.lnr[i_rad])
    logr_conv = get_value(sm.props.lnr[i_conv])
    logP_rad  = get_value(sm.props.lnP_face[i_rad])
    logP_conv = get_value(sm.props.lnP_face[i_conv])
    logρ_rad  = get_value(sm.props.lnρ_face[i_rad])
    logρ_conv = get_value(sm.props.lnρ_face[i_conv])
    m_rad  = sm.props.m[i_rad] 
    m_conv = sm.props.m[i_conv]
    dnabla_rad  = get_value(sm.props.turb_res[i_rad].∇ᵣ) - get_value(sm.props.∇ₐ_face[i_rad])
    dnabla_conv = get_value(sm.props.turb_res[i_conv].∇ᵣ) - get_value(sm.props.∇ₐ_face[i_conv])
    dnabla_total = dnabla_rad - dnabla_conv
    if abs(dnabla_total) < 1e-12 
        logr_sch = (logr_rad + logr_conv) / 2.0
        logP_sch = (logP_rad + logP_conv) / 2.0
        logρ_sch = (logρ_rad + logρ_conv) / 2.0
        m_sch = (m_rad + m_conv) / 2.0
    else
        logr_sch = logr_conv - dnabla_conv * (logr_rad - logr_conv) / dnabla_total
        logP_sch = logP_conv - dnabla_conv * (logP_rad - logP_conv) / dnabla_total
        logρ_sch = logρ_conv - dnabla_conv * (logρ_rad - logρ_conv) / dnabla_total
        m_sch = m_conv - dnabla_conv * (m_rad - m_conv) / dnabla_total
    end
    return logr_sch, logP_sch, logρ_sch, m_sch
end

function interpolate_ov_boundary(sm::StellarModel, i_ov:: Int, threshold:: Float64)
    i_d = i_ov - 1 
    logr_d  = get_value(sm.props.lnr[i_d])
    logr_ov = get_value(sm.props.lnr[i_ov])
    D_ov = get_D_turb(sm, i_ov)
    D_d = get_D_turb(sm, i_d)
    D_ov_safe = max(D_ov, 1e-20)
    D_d_safe = max(D_d, 1e-20)
    target_logD = log(threshold)
    logD_ov = log(D_ov_safe)
    logD_d = log(D_d_safe)
    diff_logD = logD_ov - logD_d
    diff_r = logr_ov - logr_d
    if abs(diff_logD) < 1e-15
        return (logr_ov + logr_d) / 2.0
    end
    slope = diff_logD / diff_r
    return logr_d + (target_logD - logD_d) / slope
end

function calculate_overshoot_length(sm:: StellarModel) 
    OVERSHOOT_THRESHOLD = 1e5
    nz_interior  = sm.props.nz - 1 
    i_Sch = nothing
    for k in 2:nz_interior
        nabla_ad = get_value(sm.props.∇ₐ_face[k])
        nabla_rad = get_value(sm.props.turb_res[k].∇ᵣ)
        if nabla_ad > nabla_rad
            i_Sch = k
            break
        end   
    end 
    if isnothing(i_Sch)
         return 0.0
    end
    logr_sch, logP_sch, logρ_sch, m_sch = interpolate_boundary_live(sm, i_Sch)
    HP_Sch = exp(logP_sch) / (exp(logρ_sch) * CGRAV * m_sch / exp(logr_sch)^2)
    i_ov = nothing
    for k = i_Sch + 1:nz_interior
        if get_D_turb(sm, k) < OVERSHOOT_THRESHOLD
            i_ov = k
            break
        end
    end
    if isnothing(i_ov)
        return (0.0, ForwardDiff.value(exp(logr_sch)), ForwardDiff.value(HP_Sch), 0.0, 0.0)
    end
    r_overshoot = exp(interpolate_ov_boundary(sm, i_ov, OVERSHOOT_THRESHOLD))
    overshooting_distance_cm = abs(r_overshoot - exp(logr_sch))
    return (ForwardDiff.value(overshooting_distance_cm / HP_Sch), ForwardDiff.value(exp(logr_sch)), ForwardDiff.value(HP_Sch), ForwardDiff.value(r_overshoot), ForwardDiff.value(overshooting_distance_cm))
end

function get_T_eff(sm::StellarModel)
    i_surf = sm.props.nz
    L_surf = get_00_dual(sm.props.L[i_surf]) * LSUN
    r_surf = exp(get_00_dual(sm.props.lnr[i_surf]))
    T_eff = (L_surf / (4 * pi * r_surf^2 * SIGMA_SB))^(0.25)
    
    return ForwardDiff.value(T_eff)
end

##
function gammaTurb_arcsin_res(sm::StellarModel, k::Int)
    ###Constants###
    C_d = 8/3 * sqrt(2/3)
    γ_face_00 = abs(get_00_dual(sm.props.gamma_turb[k]))
    ω_face_00 = (sinh(γ_face_00))^2
    dm_cell_00 = sm.props.dm[k]
    m_cell_00 = sm.props.m[k]
    # variables defined at face which are same in all cases
    r_face_00 = exp(get_00_dual(sm.props.lnr[k]))
    L_face_00  = get_00_dual(sm.props.L[k]) * LSUN

    if k == 1 
        κ_inner_face_00  = get_00_dual(sm.props.κ[k])
        P_inner_face_00 = get_00_dual(sm.props.eos_res[k].P)
        ρ_inner_face_00 = get_00_dual(sm.props.eos_res[k].ρ)
        T_inner_face_00 = get_00_dual(sm.props.eos_res[k].T)
        ∇ₐ_inner_face_00 = get_00_dual(sm.props.eos_res[k].∇ₐ)
        cₚ_inner_face_00 = get_00_dual(sm.props.eos_res[k].cₚ)

        # cell centre values of the upper cell (p1)
        κ_cc_p1  = get_p1_dual(sm.props.κ[k+1])
        P_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].P)
        ρ_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].ρ)
        T_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].T)
        ∇ₐ_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].∇ₐ)
        cₚ_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].cₚ)

        # thermodynamic variables at face (00 cell)
        κ_face_00 = κ_inner_face_00 + (κ_cc_p1 - κ_inner_face_00)*(m_cell_00/(sm.props.dm[k] + 0.5* sm.props.dm[k+1]))
        P_face_00 = P_inner_face_00 + (P_cc_p1 - P_inner_face_00)*(m_cell_00/(sm.props.dm[k] + 0.5* sm.props.dm[k+1]))
        ρ_face_00 = ρ_inner_face_00 + (ρ_cc_p1 - ρ_inner_face_00)*(m_cell_00/(sm.props.dm[k] + 0.5* sm.props.dm[k+1]))
        T_face_00 = T_inner_face_00 + (T_cc_p1 - T_inner_face_00)*(m_cell_00/(sm.props.dm[k] + 0.5* sm.props.dm[k+1]))
        cₚ_face_00 = cₚ_inner_face_00 + (cₚ_cc_p1 - cₚ_inner_face_00)*(m_cell_00/(sm.props.dm[k] + 0.5* sm.props.dm[k+1]))
        ∇ₐ_face_00 = ∇ₐ_inner_face_00 + (∇ₐ_cc_p1 - ∇ₐ_inner_face_00)*(m_cell_00/(sm.props.dm[k] + 0.5* sm.props.dm[k+1]))

        # mixing term 
        F_00 = get_00_dual(sm.props.turb_flux_cc[k])
        mixing_term = (F_00/  (sm.props.dm[k]+0.5*sm.props.dm[k+1]))

    elseif k == sm.props.nz 

        # thermodynamic variables at face (00 cell)
        P_face_00  = get_00_dual(sm.props.eos_res[k].P)
        ρ_face_00  = get_00_dual(sm.props.eos_res[k].ρ)
        T_face_00  = get_00_dual(sm.props.eos_res[k].T)
        κ_face_00  = get_00_dual(sm.props.κ[k]) 
        ∇ₐ_face_00 = get_00_dual(sm.props.eos_res[k].∇ₐ)
        cₚ_face_00 = get_00_dual(sm.props.eos_res[k].cₚ)

        # mixing term 
        F_m1 = get_m1_dual(sm.props.turb_flux_cc[k-1])
        mixing_term = -(F_m1/  (0.5*sm.props.dm[k])) 
    else 
        # thermodynamic variables at face (00 cell )
        P_face_00 = exp(get_00_dual(sm.props.lnP_face[k]))
        ρ_face_00 = exp(get_00_dual(sm.props.lnρ_face[k]))
        T_face_00 = exp(get_00_dual(sm.props.lnT_face[k]))
        κ_face_00 = get_00_dual(sm.props.κ_face[k])
        ∇ₐ_face_00 = get_00_dual(sm.props.∇ₐ_face[k])
        cₚ_face_00 = get_00_dual(sm.props.cₚ_face[k])

        #mixing term 
        F_m1 = get_m1_dual(sm.props.turb_flux_cc[k-1])
        F_00 = get_00_dual(sm.props.turb_flux_cc[k])
        #mixing_term = ((F_00 - F_m1)/  (0.5*sm.props.dm[k-1] + 0.5*sm.props.dm[k]))
        mixing_term = ((F_00 - F_m1)/  (0.5*sm.props.dm[k] + 0.5*sm.props.dm[k+1]))
    end 

    # variables derived from primary variables at face (00 cell)
    Hₚ_face_00 = P_face_00 / (ρ_face_00 * m_cell_00 * CGRAV / r_face_00^2)
    Λ_face_00 = 1 / (1 / Hₚ_face_00 + 1 / r_face_00)
    ∇ᵣ_face_00 = (3 * κ_face_00 * L_face_00 * P_face_00) / (16π * CRAD * CLIGHT * CGRAV * m_cell_00 * T_face_00^4)
    τᵣ_face_00 = (cₚ_face_00 * κ_face_00 * ρ_face_00^2 * Λ_face_00^2) / (48 * SIGMA_SB * T_face_00^3)
    c_s_face_00 = sqrt(P_face_00 / ρ_face_00)
    k_rad_face_00 = (16 * SIGMA_SB * T_face_00^3) / (3 * κ_face_00 * ρ_face_00)
    α₂_face_00 = ρ_face_00 * cₚ_face_00 * 0.5 * sqrt(2 / 3) * Λ_face_00 * abs(sinh(γ_face_00))
    α₁_face_00 = ∇ₐ_face_00 * T_face_00 * Λ_face_00 * 0.5 * sqrt(2 / 3) * cₚ_face_00 / Hₚ_face_00^2
    SA_face_00 = (∇ᵣ_face_00 - ∇ₐ_face_00) * (1 + α₂_face_00 / k_rad_face_00)^(-1)
    beta =  α₂_face_00/k_rad_face_00
    # Terms for calculating the residual 
     gamma_old = abs(get_value(sm.start_step_props.gamma_turb[k]))
    omega_var_term = (sinh(γ_face_00)^2- sinh(gamma_old)^2) / sm.props.dt
    source_term = α₁_face_00 * SA_face_00 * abs(sinh(γ_face_00))
    turb_dissipation_term = C_d * (abs(sinh(γ_face_00)))^3 / Λ_face_00
    rad_dissipation_term = ω_face_00 / τᵣ_face_00
    excess_term = C_d * (c_s_face_00 * EXCESS_FACTOR[])^3 / Λ_face_00

    R5 =
    omega_var_term +
    turb_dissipation_term +
    rad_dissipation_term -
    source_term -
    excess_term - mixing_term

    term_scale =
    abs(omega_var_term) +
    abs(turb_dissipation_term) +
    abs(rad_dissipation_term) +
    abs(source_term) +
    abs(excess_term) + 
    abs(mixing_term)

    relative_residual = abs(R5) / term_scale
    relative_residual =
    term_scale > 0 ?
    abs(R5) / term_scale :
    abs(R5)

    return  ForwardDiff.value(omega_var_term), ForwardDiff.value(turb_dissipation_term), ForwardDiff.value(rad_dissipation_term), ForwardDiff.value(source_term), ForwardDiff.value(mixing_term), ForwardDiff.value(excess_term), ForwardDiff.value(relative_residual)
    
    end



    

##
using LaTeXStrings
StellarModels.add_profile_option!(sm, "gamma_turb_energy", "unitless", (sm, k) -> ((get_value(sm.props.gamma_turb[k]))))
StellarModels.add_profile_option!(sm, "turb_flux_cc", "unitless", (sm,k) -> ((get_value(sm.props.turb_flux_cc[k]))))
StellarModels.add_profile_option!(sm, "D_turb_kuff", "unitless", (sm, k) -> ((get_value(sm.props.D_turb[k]))))
StellarModels.add_profile_option!(sm, "kappa", "unitless", (sm, k) -> ((get_value(sm.props.κ[k]))))
StellarModels.add_profile_option!(sm, "nablaa_face", "unitless", (sm, k) -> get_value(sm.props.∇ₐ_face[k]), label="\nabla_\text{a,face}")
StellarModels.add_profile_option!(sm, "nablar_face", "unitless", (sm, k) -> get_value(sm.props.turb_res[k].∇ᵣ),  label="\nabla_\text{r,face}")
StellarModels.add_profile_option!(sm, "nabla_face_k", "unitless", (sm, k) -> get_nabla_tdc(sm,k), label="\nabla_\text{face}")
StellarModels.add_profile_option!(sm, "v_turb", "unitless", (sm, k) -> get_value(sm.props.turb_res[k].v_turb), label="\nabla_\text{r,face}")
StellarModels.add_profile_option!(sm,"omega_var_term", "unitless", (sm,k) -> gammaTurb_arcsin_res(sm,k)[1], label="\nabla_\text{r,face}")
StellarModels.add_profile_option!(sm,"turb_dissipation_term", "unitless", (sm,k) -> gammaTurb_arcsin_res(sm,k)[2], label="\nabla_\text{r,face}")
StellarModels.add_profile_option!(sm,"rad_dissipation_term", "unitless", (sm,k) -> gammaTurb_arcsin_res(sm,k)[3], label="\nabla_\text{r,face}")
StellarModels.add_profile_option!(sm,"source_term", "unitless", (sm,k) -> gammaTurb_arcsin_res(sm,k)[4], label="\nabla_\text{r,face}")
StellarModels.add_profile_option!(sm,"mixing_term", "unitless", (sm,k) -> gammaTurb_arcsin_res(sm,k)[5], label="\nabla_\text{r,face}")
StellarModels.add_profile_option!(sm,"excess_term", "unitless", (sm,k) -> gammaTurb_arcsin_res(sm,k)[6], label="\nabla_\text{r,face}")
StellarModels.add_profile_option!(sm,"rel_res", "unitless", (sm,k) -> gammaTurb_arcsin_res(sm,k)[7], label="\nabla_\text{r,face}")
# StellarModels.add_profile_option!(sm,"lambda", "unitless", (sm,k) -> gammaTurb_arcsin_res(sm,k)[10], label="\nabla_\text{r,face}")
# StellarModels.add_profile_option!(sm,"A_p1", "unitless", (sm,k) -> gammaTurb_arcsin_res(sm,k)[7], label="\nabla_\text{r,face}")
# StellarModels.add_profile_option!(sm,"A_00", "unitless", (sm,k) -> gammaTurb_arcsin_res(sm,k)[8], label="\nabla_\text{r,face}")
# StellarModels.add_profile_option!(sm,"F_p1", "unitless", (sm,k) -> gammaTurb_arcsin_res(sm,k)[9], label="\nabla_\text{r,face}")
# StellarModels.add_profile_option!(sm,"F_00", "unitless", (sm,k) -> gammaTurb_arcsin_res(sm,k)[10], label="\nabla_\text{r,face}")
# StellarModels.add_profile_option!(sm,"t_p1", "unitless", (sm,k) -> gammaTurb_arcsin_res(sm,k)[11], label="\nabla_\text{r,face}")
# StellarModels.add_profile_option!(sm,"t_00", "unitless", (sm,k) -> gammaTurb_arcsin_res(sm,k)[12], label="\nabla_\text{r,face}")
# StellarModels.add_profile_option!(sm,"t_m1", "unitless", (sm,k) -> gammaTurb_arcsin_res(sm,k)[13], label="\nabla_\text{r,face}")
StellarModels.add_history_option!(sm, "alpha_overshoot", "H_p", sm ->calculate_overshoot_length(sm)[1], label=L"\text{Alpha_ov}\,[\alpha_{ov}]")
# StellarModels.add_history_option!(sm, "Sch_radius", "Rsun", sm ->calculate_overshoot_length(sm)[2], label=L"\text{Alpha_ov}\,[\alpha_{ov}]")
# StellarModels.add_history_option!(sm, "pressure_scale_height", "unitless", sm ->calculate_overshoot_length(sm)[3], label=L"\text{Alpha_ov}\,[\alpha_{ov}]")
# StellarModels.add_history_option!(sm, "ov_radius", "unitless", sm ->calculate_overshoot_length(sm)[4], label=L"\text{Alpha_ov}\,[\alpha_{ov}]")
# StellarModels.add_history_option!(sm, "ov_distance", "unitless", sm ->calculate_overshoot_length(sm)[5], label=L"\text{Alpha_ov}\,[\alpha_{ov}]")
StellarModels.add_history_option!(sm, "T_eff", "unitless", sm ->get_T_eff(sm), label=L"\text{Alpha_ov}\,[\alpha_{ov}]")
##
# "tau", "rad_diss", "vis_diss", "Lambda", "excess_term"
open("example_options.toml", "w") do file
    write(file,
          """
          [remesh]
          do_remesh = false
          [solver]
          newton_max_iter_first_step = 5000
          initial_model_scale_max_correction = 0.2
          newton_max_iter = 50
          scale_max_correction = 1.0
          solver_progress_iter = 1
          relative_correction_tolerance= 1e14
          maximum_residual_tolerance = 1e-2
          use_preconditioning = true
          report_solver_progress = true
          [timestep]
          dt_max_increase = 1.5
          delta_R_limit = 0.01
          delta_Tc_limit = 0.01
          delta_Xc_limit = 0.005

          [termination]
          max_model_number = 20000
          max_center_T = 1e12
          min_center_H1 = 0.5

          [io]
          profile_interval = 20
          hdf5_history_filename = "history_models_5MSUN_woking.hdf5"
          hdf5_profile_filename = "profiles_models_5MSUN_working.hdf5"
          terminal_header_interval = 100
          terminal_info_interval = 100
          profile_values = ["zone", "mass", "dm", "log10_rho", "log10_r", "log10_P", "log10_T", "luminosity", "X", "Y","D_face", "nablaa_face", "nablar_face","nabla_face", "gamma_turb_energy", "D_turb_kuff", "nabla_face_k", "kappa", "v_turb", "omega_var_term", "turb_dissipation_term", "rad_dissipation_term", "source_term", "excess_term", "rel_res"]
          history_values = ["model_number", "age", "dt", "star_mass","X_center", "T_surf", "T_eff", "L_surf","alpha_overshoot"]

    
          """)
end
StellarModels.set_options!(sm.opt, "./example_options.toml")
# rm(sm.opt.io.hdf5_history_filename; force=true)
# rm(sm.opt.io.hdf5_profile_filename; force=true)


##
#Configure live plots. To turn off one can use `plotter = Plotting.NullPlotter()`
using GLMakie
GLMakie.activate!()
set_theme!(Plotting.basic_theme())
f = Figure(size=(1400,750))
plots = [Plotting.HRPlot(f[1,1]),
         Plotting.TRhoProfile(f[1,2]),
         Plotting.KippenLine(f[2,1], xaxis=:time, time_units=:Gyr),
         Plotting.AbundancePlot(f[2,2],net,log_yscale=true, ymin=1e-3),
         Plotting.HistoryPlot(f[1,3], sm, x_name="age", y_name="X_center", othery_name="Y_center", link_yaxes=true),
         Plotting.ProfilePlot(f[2,3], sm, x_name="mass", y_name="log10_rho", othery_name="log10_T")]
plotter = Plotting.Plotter(fig=f,plots=plots)

##
plotter = Plotting.NullPlotter()
##
#set initial condition and run model
n = 1.5
# StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], 
#                                             1 * MSUN, 100 * RSUN; initial_dt=10 * SECYEAR)
Z = 0.0142 *2
Y = 0.249 + 1.67 * Z
X = 1 - Y - Z
tdc_initial_condition!(n, sm, nz, X, Z, 0.0, Chem.abundance_lists[:ASG_09], 
                                        5.0 * MSUN, 80 *  RSUN; initial_dt=0.000001 * SECYEAR)

##         
@time Evolution.do_evolution_loop!(sm, plotter=plotter); 


##
open("gamma_target_no_dual.log", "w") do io
    redirect_stdout(io) do 
        redirect_stderr(io) do
            @time Evolution.do_evolution_loop!(sm, plotter=plotter)

        end
    end
end

##
using Serialization

sm_run = deserialize(
    "/Users/rdbnath/Documents/Jems.jl/sm_checkpoint.jls"
)

@show sm_run.props.model_number
##
nz_restart = sm_run.props.nz
nextra_restart = 10000

sm = StellarModel(
    TDCEquationSet(),
    nz_restart,
    nextra_restart,
    net,
    eos_table,
    opacity,
    turbulence,
);

##
using Serialization

# The termination check occurs before dt_next is refreshed.
sm.props.dt_next = Evolution.get_dt_next(sm)

checkpoint = "/Users/rdbnath/Documents/Jems.jl/sm_checkpoint.jls"
serialize(checkpoint, sm)

##
open("example_options.toml", "w") do file
    write(file,
          """
          [remesh]
          do_remesh = true
          [solver]
          newton_max_iter_first_step = 5000
          initial_model_scale_max_correction = 0.2
          newton_max_iter = 50
          scale_max_correction = 1.0
          solver_progress_iter = 1
          relative_correction_tolerance= 1e14
          maximum_residual_tolerance = 1e-2
          use_preconditioning = true
          report_solver_progress = true
          [timestep]
          dt_max_increase = 1.5
          delta_R_limit = 0.01
          delta_Tc_limit = 0.01
          delta_Xc_limit = 0.005

          [termination]
          max_model_number = 20000
          max_center_T = 1e12
          min_center_H1 = 1e-4

          [io]
          profile_interval = 1
          hdf5_history_filename = "history_models_final.hdf5"
          hdf5_profile_filename = "profiles_models_final.hdf5"
          terminal_header_interval = 100
          terminal_info_interval = 100
          profile_values = ["zone", "mass", "dm", "log10_rho", "log10_r", "log10_P", "log10_T", "luminosity", "X", "Y","D_face", "nablaa_face", "nablar_face","nabla_face", "gamma_turb_energy", "D_turb_kuff", "nabla_face_k", "kappa"]
          history_values = ["model_number", "age", "dt", "star_mass","X_center", "T_surf", "T_eff", "L_surf","alpha_overshoot","Sch_radius", "pressure_scale_height", "ov_radius", "ov_distance"]

    
          
          

          """)
end
sm_run.opt.io.hdf5_history_keep_open = false
sm_run.opt.io.hdf5_profile_keep_open = false
sm_run.output_files_created = false
StellarModels.set_options!(sm_run.opt, "./example_options.toml")
##
sm_run.props.dt_next = Jems.Evolution.get_dt_next(sm_run)
##
using GLMakie
GLMakie.activate!()
set_theme!(Plotting.basic_theme())
f = Figure(size=(1400,750))
plots = [Plotting.HRPlot(f[1,1]),
         Plotting.TRhoProfile(f[1,2]),
         Plotting.KippenLine(f[2,1], xaxis=:time, time_units=:Gyr),
         Plotting.AbundancePlot(f[2,2],net,log_yscale=true, ymin=1e-3),
         Plotting.HistoryPlot(f[1,3], sm, x_name="age", y_name="X_center", othery_name="Y_center", link_yaxes=true),
         Plotting.ProfilePlot(f[2,3], sm, x_name="mass", y_name="log10_rho", othery_name="log10_T")]
plotter = Plotting.Plotter(fig=f,plots=plots)

@time Jems.Evolution.do_evolution_loop!(
    sm_run;
    plotter = plotter
)

##
using Printf

function check_tdc_jacobian(
    sm;
    zones = [89],
    equations = [2],
    variables = [:lnρ, :lnT, :lnr, :lum, :gamma_turb],
    steps = (1e-4, 1e-5, 1e-6),
)
    @assert sm.props !== sm.start_step_props
    @assert sm.props.nz == sm.start_step_props.nz
    @assert sm.props.dt > 0

    nz, nv = sm.props.nz, sm.nvars
    n = nz * nv

    @assert all(k -> 1 <= k <= nz, zones)
    @assert all(eq -> 1 <= eq <= nv, equations)
    @assert all(s -> s > 0, steps)

    x = @view sm.props.ind_vars[1:n]
    saved_x = copy(x)
    results = NamedTuple[]

    function evaluate!()
        Jems.StellarModels.evaluate_stellar_model_properties!(
            sm, sm.props
        )
        Jems.Evolution.eval_jacobian_eqs!(sm)
        return nothing
    end

    try
        evaluate!()

        # Preserve the Jacobian before perturbing anything.
        baseline = Dict(
            k => (
                L = copy(sm.solver_data.jacobian_L[k]),
                D = copy(sm.solver_data.jacobian_D[k]),
                U = copy(sm.solver_data.jacobian_U[k]),
            )
            for k in zones
        )

        f0 = copy(@view sm.solver_data.eqs_numbers[1:n])

        for k in zones
            @printf(
                "\nEquation zone %d: m/M=%.8e\n",
                k, sm.props.m[k] / sm.props.m[nz],
            )

            for eq in equations
                @printf(
                    "Baseline F[%d,%d] = % .12e\n",
                    k, eq, f0[(k - 1) * nv + eq],
                )
            end

            for cell in max(1, k - 1):min(nz, k + 1)
                block = cell < k ? baseline[k].L :
                        cell > k ? baseline[k].U :
                                   baseline[k].D

                for variable in variables
                    j = sm.vari[variable]
                    index = (cell - 1) * nv + j
                    x0 = saved_x[index]

                    # Log variables use absolute log increments.
                    # Luminosity is scaled to its local magnitude.
                    scale = variable == :lum ?
                            max(abs(x0), 1e-8) : 1.0

                    for step in steps
                        h = step * scale

                        # Stay on the same side of gamma=0.
                        # Relative gamma increments also separate
                        # the three tests when gamma is very small.
                        if variable == :gamma_turb
                            if x0 == 0
                                println(
                                    "Skipping gamma at cell ", cell,
                                    ": derivative at zero may be nonsmooth.",
                                )
                                break
                            end

                            h = min(
                                step * min(abs(x0), 1.0),
                                0.1 * abs(x0),
                            )
                        end

                        xp, xm = x0 + h, x0 - h
                        if !(xm < x0 < xp)
                            println(
                                "Skipping unrepresentable step: cell=",
                                cell, " variable=", variable, " h=", h,
                            )
                            continue
                        end

                        copyto!(x, saved_x)
                        x[index] = xp
                        evaluate!()
                        fp = [
                            sm.solver_data.eqs_numbers[
                                (k - 1) * nv + eq
                            ]
                            for eq in equations
                        ]

                        copyto!(x, saved_x)
                        x[index] = xm
                        evaluate!()
                        fm = [
                            sm.solver_data.eqs_numbers[
                                (k - 1) * nv + eq
                            ]
                            for eq in equations
                        ]

                        for (a, eq) in enumerate(equations)
                            fbase = f0[(k - 1) * nv + eq]
                            ad = block[eq, j]
                            fd = (fp[a] - fm[a]) / (xp - xm)

                            forward = (fp[a] - fbase) / (xp - x0)
                            backward = (fbase - fm[a]) / (x0 - xm)

                            absolute_error = abs(ad - fd)
                            relative_error = absolute_error /
                                max(abs(ad), abs(fd), 1e-30)

                            push!(results, (
                                zone = k,
                                equation = eq,
                                cell = cell,
                                variable = variable,
                                h = h,
                                AD = ad,
                                FD = fd,
                                forward = forward,
                                backward = backward,
                                abs_error = absolute_error,
                                rel_error = relative_error,
                            ))

                            @printf(
                                "eq=%d cell=%d %-10s h=%.3e AD=% .8e FD=% .8e rel=%.3e FD+=% .8e FD-=% .8e\n",
                                eq, cell, string(variable), h,
                                ad, fd, relative_error,
                                forward, backward,
                            )
                        end
                    end
                end
            end
        end
    finally
        copyto!(x, saved_x)
        evaluate!()
    end

    return results
end



checks_gamma = check_tdc_jacobian(
    sm_run;
    zones = [89],
    equations = [5],
)
##
iz = 1014
idx = (iz-1)*sm.nvars + 5

γ_before = get_value(sm.props.gamma_turb[iz])
Δγ_raw = sm.solver_data.solver_corr[idx]
Δγ_applied = correction_multiplier * Δγ_raw

println(
    "iter=$i",
    " gamma=$γ_before",
    " signgamma=$(sign(γ_before))",
    " smooth_deriv=$(γ_before / sqrt(γ_before^2 + 1e-16))",
    " raw_dgamma=$Δγ_raw",
    " applied_dgamma=$Δγ_applied",
    " predicted_gamma=$(γ_before + Δγ_applied)"
)