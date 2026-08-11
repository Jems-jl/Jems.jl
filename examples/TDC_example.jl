using BenchmarkTools
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
include("TDC_StellarModels.jl")
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
EXCESS_FACTOR[] = 1e-6
α_w_FACTOR[] = 0.25

##
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16], [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
nz = 1000
nextra = 100
#eos = EOS.IdealEOS(true)
eos_table = EOSTableCollector("MESA_data/eosFreeEOS_data", include_radiation = true)
println("Memory used: $(Base.summarysize(eos_table) / 1024^3) GB")
low_T_collection = OpacityTableCollector("MESA_data/kap_data", "lowT_fa05_gs98") 
high_T_collection = OpacityTableCollector("MESA_data/kap_data","oplib_agss09" ) 
opacity = CompositeOpacity(low_T_collection, high_T_collection, 3.8, 4.2)
# eos = EOS.IdealEOS(true)
# opacity = Opacity.SimpleElectronScatteringOpacity()
turbulence = Turbulence.BasicMLT(2.0)
##
nz = 1000
nextra = 700
sm = StellarModel(TDCEquationSet(), nz, nextra, net, eos_table, opacity, turbulence);



##

n = 1.5
# Polytropic Initial condition
tdc_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], 
                                            100 * MSUN, 3000 * RSUN; initial_dt=10 * SECYEAR)
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
#=

To benchmark the linear solver itself we need to perform
the jacobian evaluation as a setup for the benchmark. This is because the solver
destroys the Jacobian to perform in-place operations.
=#

@benchmark begin
    Evolution.block_tridiagonal_solver!($sm, $sm.solver_data)
end setup=(Evolution.eval_jacobian_eqs!($sm))
##

function get_D_turb(sm::StellarModel, k:: Int)
    α_w   = 0.25
    r = exp(get_00_dual(sm.props.lnr[k]))
    P = exp(get_00_dual(sm.props.lnP_face[k]))
    ρ = exp(get_00_dual(sm.props.lnρ_face[k]))
    
    γ_face_00 = get_00_dual(sm.props.gamma_turb[k])
    ω = exp(γ_face_00)
    m = sm.props.m[k]

    g = CGRAV * m / (r^2)
    Hp = P / ( ρ * g)
    Λ = 1/(1/Hp + 1/r)
    v = sqrt(2 * ω)
    D = α_w * Λ * sqrt(2 * ω)

    return ForwardDiff.value(D)

end 

function get_lambda_turb(sm::StellarModel, k:: Int)
    α_w   = 0.25
    r = exp(get_00_dual(sm.props.lnr[k]))
    P = exp(get_00_dual(sm.props.lnP_face[k]))
    ρ = exp(get_00_dual(sm.props.lnρ_face[k]))
    
    γ_face_00 = get_00_dual(sm.props.gamma_turb[k])
    ω = exp(γ_face_00)
    m = sm.props.m[k]

    g = CGRAV * m / (r^2)
    Hp = P / ( ρ * g)
    Λ = 1/(1/Hp + 1/r)
    v = sqrt(2 * ω)
    D = α_w * Λ * sqrt(2 * ω)

    return ForwardDiff.value(Λ)

end 

"""
Function for calculating nabla_tdc for profile output
"""

function calculate_nabla_tdc(sm:: StellarModel, k :: Int)
    L = get_00_dual(sm.props.L[k]) * LSUN
    γ₀ = get_00_dual(sm.props.gamma_turb[k])
    ω = exp(γ₀)
    m₀ = sm.props.m[k]
    r₀ = exp(get_00_dual(sm.props.lnr[k]))
    if k == sm.props.nz
        P = get_00_dual(sm.props.eos_res[k].P)
        ρ = get_00_dual(sm.props.eos_res[k].ρ)
        T = get_00_dual(sm.props.eos_res[k].T)
        κ = get_00_dual(sm.props.κ[k])
        cₚ = get_00_dual(sm.props.eos_res[k].cₚ)
        Hₚ = P / (ρ * CGRAV * m₀ / r₀^2)
        Λ = 1/(1/Hₚ + 1/r₀)
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
    Λ = 1/(1/Hₚ + 1/r₀)
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
    # Calls the calculation function and extracts the pure numerical value.
    return calculate_nabla_tdc(sm, k).value
end
"""
End of function
"""

"""
Beginning of function for setting up alpha_overshoot history functions
"""

function interpolate_boundary_live(sm::StellarModel, i_rad::Int)
    i_conv = i_rad - 1
    
    #Radius 
    logr_rad  = get_value(sm.props.lnr[i_rad])
    logr_conv = get_value(sm.props.lnr[i_conv])
    #Pressure
    logP_rad  = get_value(sm.props.lnP_face[i_rad])
    logP_conv = get_value(sm.props.lnP_face[i_conv])
    #density 
    logρ_rad  = get_value(sm.props.lnρ_face[i_rad])
    logρ_conv = get_value(sm.props.lnρ_face[i_conv])
    #mass
    m_rad  = sm.props.m[i_rad] 
    m_conv = sm.props.m[i_conv]

    # Delta-Nabla (Δ∇ = ∇_rad - ∇_ad)
    dnabla_rad  = get_value(sm.props.turb_res[i_rad].∇ᵣ) - get_value(sm.props.∇ₐ_face[i_rad])
    dnabla_conv = get_value(sm.props.turb_res[i_conv].∇ᵣ) - get_value(sm.props.∇ₐ_face[i_conv])
    
    # Linear Interpolation Calculation
    dnabla_total = dnabla_rad - dnabla_conv
    
    if abs(dnabla_total) < 1e-12 
        # Avoid division by zero: return the midpoint if profiles are flat
        logr_sch = (logr_rad + logr_conv) / 2.0
        logP_sch = (logP_rad + logP_conv) / 2.0
        logρ_sch = (logρ_rad + logρ_conv) / 2.0
        m_sch = (m_rad + m_conv) / 2.0
    else
        # Interpolation: logr_sch = logr_conv - Δ∇_conv * (Δlogr / Δ(Δ∇))
        logr_sch = logr_conv - dnabla_conv * (logr_rad - logr_conv) / dnabla_total
        logP_sch = logP_conv - dnabla_conv * (logP_rad - logP_conv) / dnabla_total
        logρ_sch = logρ_conv - dnabla_conv * (logρ_rad - logρ_conv) / dnabla_total
        m_sch = m_conv - dnabla_conv * (m_rad - m_conv) / dnabla_total
    end
    
    # Return the interpolated natural log radius (ln r_sch)
    return logr_sch, logP_sch, logρ_sch, m_sch
end

function interpolate_ov_boundary(sm::StellarModel, i_ov:: Int)
    i_d = i_ov - 1 
    logr_d  = get_value(sm.props.lnr[i_d])
    logr_ov = get_value(sm.props.lnr[i_ov])

    # D_ov = get_value(sm.props.D_turb[i_ov])
    # D_d = get_value(sm.props.D_turb[i_d])

    D_ov = get_D_turb(sm, i_ov)
    D_d = get_D_turb(sm, i_d)
    #Logarithmic interpolation 
    #Putting a safer lower limit 
    D_ov_safe = max(D_ov, 1e-20)
    D_d_safe = max(D_d, 1e-20)
    target_logD = log(1e5)
    logD_ov = log(D_ov_safe)
    logD_d = log(D_d_safe)

    # diff_D = D_ov - D_d
    diff_logD = logD_ov - logD_d
    diff_r = logr_ov - logr_d

    # Avoid div by zero
    if abs(diff_logD) < 1e-15
        return (logr_ov + logr_d) / 2.0
    end

    slope = diff_logD / diff_r

    logr_boundary = logr_d + (target_logD - logD_d) / slope
    return logr_boundary
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
    
    # Interpolate to find boundary properties
    logr_sch, logP_sch, logρ_sch, m_sch = interpolate_boundary_live(sm, i_Sch)
    P_boundary = exp(logP_sch)
    ρ_boundary = exp(logρ_sch)
    r_boundary = exp(logr_sch) # Radius in cm
    m_boundary = m_sch
    
    HP_Sch = P_boundary / (ρ_boundary * CGRAV * m_boundary / r_boundary^2)

    
    i_ov = nothing

    #calculate overshoot length using interpolation 
    for k = i_Sch + 1:nz_interior
        D_turb = get_D_turb(sm,k)
        if D_turb < OVERSHOOT_THRESHOLD
            i_ov = k
            break
        end
    end
    
    if isnothing(i_ov)
        # If mixing doesn't fall below threshold before the surface
        return (0.0,  ForwardDiff.value(r_boundary), ForwardDiff.value(HP_Sch),0.0, 0.0)
    end
    logr_ov = interpolate_ov_boundary(sm, i_ov) #Redundant calculations 
    r_overshoot = exp(logr_ov) # Radius in cm 
   
    overshooting_distance_cm = abs(r_overshoot - r_boundary)
    alpha_ov = overshooting_distance_cm / HP_Sch
    
    return (ForwardDiff.value(alpha_ov), ForwardDiff.value(r_boundary), ForwardDiff.value(HP_Sch), ForwardDiff.value(r_overshoot), ForwardDiff.value(overshooting_distance_cm))
end

##

function calculate_tau_CD_by_lambda(sm :: StellarModel, k :: Int)

###Constants###
C_d = 8/3 * sqrt(2/3)

if k == 1
    γ_face_00 = get_00_dual(sm.props.gamma_turb[k]) #defined at the outer face 
    ω_face_00 = exp(γ_face_00)   #defined at the outer face 
    dm_cell_00 = sm.props.dm[k]  
    m_cell_00 = sm.props.m[k]
    r_face_00 = exp(get_00_dual(sm.props.lnr[k])) #defined at the outer face 
    κ_face_00  = get_00_dual(sm.props.κ[k])  #defined at outer face 
    L_face_00  = get_00_dual(sm.props.L[k]) * LSUN #defined at outer face 

    P_inner_face_00 = get_00_dual(sm.props.eos_res[k].P) #defined at the inner face 
    P_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].P)
    ρ_inner_face_00 = get_00_dual(sm.props.eos_res[k].ρ) #defined at the inner face 
    ρ_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].ρ)
    T_inner_face_00 = get_00_dual(sm.props.eos_res[k].T) #defined at the inner face 
    T_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].T) 
    ∇ₐ_inner_face_00 = get_00_dual(sm.props.eos_res[k].∇ₐ) #defined at the inner face 
    ∇ₐ_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].∇ₐ)
    cₚ_inner_face_00 = get_00_dual(sm.props.eos_res[k].cₚ) #defined at the inner face 
    cₚ_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].cₚ)

    ## Interpolating the values to get the outer face values ##
    P_face_00 = P_inner_face_00 + (P_cc_p1 - P_inner_face_00)*(m_cell_00/(sm.props.m[k] + 0.5* sm.props.m[k+1]))
    ρ_face_00 = ρ_inner_face_00 + (ρ_cc_p1 - ρ_inner_face_00)*(m_cell_00/(sm.props.m[k] + 0.5* sm.props.m[k+1]))
    T_face_00 = T_inner_face_00 + (T_cc_p1 - T_inner_face_00)*(m_cell_00/(sm.props.m[k] + 0.5* sm.props.m[k+1]))
    cₚ_face_00 = cₚ_inner_face_00 + (cₚ_cc_p1 - cₚ_inner_face_00)*(m_cell_00/(sm.props.m[k] + 0.5* sm.props.m[k+1]))
    ∇ₐ_face_00 = ∇ₐ_inner_face_00 + (∇ₐ_cc_p1 - ∇ₐ_inner_face_00)*(m_cell_00/(sm.props.m[k] + 0.5* sm.props.m[k+1]))
    
    Hₚ_face_00 = P_face_00 / (ρ_face_00 * m_cell_00 * CGRAV / r_face_00^2)
    Λ_face_00 = 1 / (1 / Hₚ_face_00 + 1 / r_face_00)
    ∇ᵣ_face_00 = (3 * κ_face_00 * L_face_00 * P_face_00) / (16π * CRAD * CLIGHT * CGRAV * m_cell_00 * T_face_00^4)
    τᵣ_face_00 = (cₚ_face_00 * κ_face_00 * ρ_face_00^2 * Λ_face_00^2) / (48 * SIGMA_SB * T_face_00^3)
    c_s_face_00 = sqrt(P_face_00 / ρ_face_00)
    k_rad_face_00 = (16 * SIGMA_SB * T_face_00^3) / (3 * κ_face_00 * ρ_face_00)  ##
    α₂_face_00 = ρ_face_00 * cₚ_face_00 * 0.5 * sqrt(2 / 3) * Λ_face_00 * sqrt(ω_face_00)
    α₁_face_00 = ∇ₐ_face_00 * T_face_00 * Λ_face_00 * 0.5 * sqrt(2 / 3) * cₚ_face_00 / Hₚ_face_00^2
    SA_face_00 = (∇ᵣ_face_00 - ∇ₐ_face_00) * (1 + α₂_face_00 / k_rad_face_00)^(-1)
    dgammadt_face_00 = (ω_face_00- exp(get_value(sm.start_step_props.gamma_turb[k]))) / sm.props.dt
    
    """
    mixing term = 1/ dm(i,face) [ A(2) * (ω(2)- ω(1))/ dm(2,cell) ]

    """

    ### p1 terms needed as F_1 = 0 and F_2 = 0 is the next cell so everything is wrt to i+1 so the p1 cell, also mixing term calcualted at cell centre ###
    P_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].P)
    ρ_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].ρ)
    T_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].T)
    dm_face_p1 = 0.5*(sm.props.dm[k] + sm.props.dm[k+1])##
    dm_cell_p1 = sm.props.dm[k+1] 
    m_face_p1 = 0.5*(sm.props.m[k] + sm.props.m[k+1])
    L_cc_p1 = 0.5*(get_p1_dual(sm.props.L[k+1]) + get_00_dual(sm.props.L[k])) * LSUN
    r_cc_p1 = 0.5*(exp(get_p1_dual(sm.props.lnr[k+1])) + exp(get_00_dual(sm.props.lnr[k])))
    Hₚ_cc_p1 = P_cc_p1 / (ρ_cc_p1 * m_face_p1 * CGRAV / r_cc_p1^2)
    ω_cc_p1 = 0.5*(exp(get_p1_dual(sm.props.gamma_turb[k+1])) + exp(get_00_dual(sm.props.gamma_turb[k])))

    Λ_cc_p1 = 1 / (1 / Hₚ_cc_p1 + 1 / r_cc_p1)
    A_p1 = (4π  *ρ_cc_p1* r_cc_p1^2)^2 * Λ_cc_p1 * α_w_FACTOR[] * sqrt(ω_cc_p1)

    F_p1 = (A_p1 / dm_cell_p1) * (exp(get_p1_dual(sm.props.gamma_turb[k+1])) - exp(get_00_dual(sm.props.gamma_turb[k]))) 

    # Different terms for residual at k = 1
    mixing_term =  (F_p1/  dm_face_p1)  
    omega_var_term = ((γ_face_00- get_value(sm.start_step_props.gamma_turb[k])) / sm.props.dt) * ω_face_00
    source_term = α₁_face_00 * SA_face_00 *sqrt(ω_face_00)
    turb_dissipation_term = C_d * (ω_face_00)^(3/2) / Λ_face_00
    rad_dissipation_term = ω_face_00 / τᵣ_face_00
    excess_term = (c_s_face_00)^2 * EXCESS_FACTOR[]/ τᵣ_face_00

    return  ForwardDiff.value(τᵣ_face_00), ForwardDiff.value(rad_dissipation_term), ForwardDiff.value(turb_dissipation_term),  ForwardDiff.value(Λ_face_00), ForwardDiff.value(excess_term)
end

if k == sm.props.nz

    # ==============================================================================
    # THERMODYNAMICS & GEOMETRY (From EOS Results = Face Values)
    # ==============================================================================
    # As per instruction: EOS results here are defined at the face
    γ_face_00 = get_00_dual(sm.props.gamma_turb[k])
    ω_face_00 = exp(γ_face_00)
    ω_cc_00 = 0.5*(exp(get_00_dual(sm.props.gamma_turb[k])) + exp(get_m1_dual(sm.props.gamma_turb[k-1])))
    dm_cell_00 = sm.props.dm[k]
    m_face_00  = sm.props.m[k] # Mass at the outer boundary
    m_cc_00 = 0.5*(sm.props.m[k] + sm.props.m[k-1])
    r_face_00  = exp(get_00_dual(sm.props.lnr[k]))
    r_cc_00 = 0.5*(exp(get_00_dual(sm.props.lnr[k])) + exp(get_m1_dual(sm.props.lnr[k-1])))
    L_face_00  = get_00_dual(sm.props.L[k]) * LSUN
    
    # Thermodynamics directly from EOS (treated as Face values)
    P_face_00  = get_00_dual(sm.props.eos_res[k].P)
    ρ_face_00  = get_00_dual(sm.props.eos_res[k].ρ)
    T_face_00  = get_00_dual(sm.props.eos_res[k].T)
    
    # Opacity and gradients (Assuming these are available in props or eos_res)
    # If kappa/nabla are in eos_res, use those. If in props arrays, fetch index k.
    κ_face_00  = get_00_dual(sm.props.κ[k]) 
    ∇ₐ_face_00 = get_00_dual(sm.props.eos_res[k].∇ₐ)
    cₚ_face_00 = get_00_dual(sm.props.eos_res[k].cₚ)

    Hₚ_face_00 = P_face_00 / (ρ_face_00 * m_cc_00 * CGRAV / r_cc_00^2)
    Λ_face_00  = 1 / (1 / Hₚ_face_00 + 1 / r_cc_00)
    
    ∇ᵣ_face_00 = (3 * κ_face_00 * L_face_00 * P_face_00) / (16π * CRAD * CLIGHT * CGRAV * m_face_00 * T_face_00^4)
    τᵣ_face_00 = (cₚ_face_00 * κ_face_00 * ρ_face_00^2 * Λ_face_00^2) / (48 * SIGMA_SB * T_face_00^3)
    c_s_face_00 = sqrt(P_face_00 / ρ_face_00)
    k_rad_face_00 = (16 * SIGMA_SB * T_face_00^3) / (3 * κ_face_00 * ρ_face_00)
    
    α₂_face_00 = ρ_face_00 * cₚ_face_00 * 0.5 * sqrt(2 / 3) * Λ_face_00 * sqrt(ω_face_00)
    α₁_face_00 = ∇ₐ_face_00 * T_face_00 * Λ_face_00 * 0.5 * sqrt(2 / 3) * cₚ_face_00 / Hₚ_face_00^2
    SA_face_00 = (∇ᵣ_face_00 - ∇ₐ_face_00) * (1 + α₂_face_00 / k_rad_face_00)^(-1)
    
    #dgammadt_face_00 = (ω_face_00 - exp(get_value(sm.start_step_props.gamma_turb[k]))) / sm.props.dt
    dgammadt_face_00 = (γ_face_00- get_value(sm.start_step_props.gamma_turb[k])) / sm.props.dt
    # ==============================================================================
    #  FLUX CALCULATION (A_00)
    # ==============================================================================
    # We use the same Face values for A_00 as they are the definitive properties at k
    A_00 = (4π * ρ_face_00 * r_cc_00^2)^2 * Λ_face_00 * α_w_FACTOR[] * sqrt(ω_cc_00)
    
    # Flux entering from k-1
    F_00 = (A_00 / sm.props.dm[k]) * (exp(get_00_dual(sm.props.gamma_turb[k])) - exp(get_m1_dual(sm.props.gamma_turb[k-1])))

    # ==============================================================================
    # RESIDUAL
    # ==============================================================================
    mixing_term = -(F_00 / sm.props.dm[k]) # Flux out (F_p1) is zero at surface
    omega_var_term = dgammadt_face_00 * ω_face_00
    
    source_term = α₁_face_00 * SA_face_00 * sqrt(ω_face_00)
    turb_dissipation_term = C_d * (ω_face_00)^(3/2) / Λ_face_00
    rad_dissipation_term = ω_face_00 / τᵣ_face_00
    excess_term = (c_s_face_00)^2 * EXCESS_FACTOR[]/ τᵣ_face_00
    return ForwardDiff.value(τᵣ_face_00), ForwardDiff.value(rad_dissipation_term), ForwardDiff.value(turb_dissipation_term),  ForwardDiff.value(Λ_face_00), ForwardDiff.value(excess_term)
end



### Other Calculations : 1 < k < nz ###
begin 

    ### face Values are required for all other terms except the mixing term ###
    γ_face_00 = get_00_dual(sm.props.gamma_turb[k])
    ω_face_00 = exp(γ_face_00) ###
    dm_cell_00 = sm.props.dm[k] 
    m_cell_00 = sm.props.m[k]
    r_face_00 = exp(get_00_dual(sm.props.lnr[k]))
    P_face_00 = exp(get_00_dual(sm.props.lnP_face[k]))
    ρ_face_00 = exp(get_00_dual(sm.props.lnρ_face[k]))
    T_face_00 = exp(get_00_dual(sm.props.lnT_face[k]))
    κ_face_00 = get_00_dual(sm.props.κ_face[k])
    ∇ₐ_face_00 = get_00_dual(sm.props.∇ₐ_face[k])
    cₚ_face_00 = get_00_dual(sm.props.cₚ_face[k])
    L_face_00 = get_00_dual(sm.props.L[k]) * LSUN

    Hₚ_face_00 = P_face_00 / (ρ_face_00 * m_cell_00 * CGRAV / r_face_00^2)
    Λ_face_00 = 1 / (1 / Hₚ_face_00 + 1 / r_face_00)
    ∇ᵣ_face_00 = (3 * κ_face_00 * L_face_00 * P_face_00) / (16π * CRAD * CLIGHT * CGRAV * m_cell_00 * T_face_00^4)
    τᵣ_face_00 = (cₚ_face_00 * κ_face_00 * ρ_face_00^2 * Λ_face_00^2) / (48 * SIGMA_SB * T_face_00^3)
    c_s_face_00 = sqrt(P_face_00 / ρ_face_00)
    k_rad_face_00 = (16 * SIGMA_SB * T_face_00^3) / (3 * κ_face_00 * ρ_face_00)  ##
    α₂_face_00 = ρ_face_00 * cₚ_face_00 * 0.5 * sqrt(2 / 3) * Λ_face_00 * sqrt(ω_face_00)
    α₁_face_00 = ∇ₐ_face_00 * T_face_00 * Λ_face_00 * 0.5 * sqrt(2 / 3) * cₚ_face_00 / Hₚ_face_00^2
    SA_face_00 = (∇ᵣ_face_00 - ∇ₐ_face_00) * (1 + α₂_face_00 / k_rad_face_00)^(-1)
    #dgammadt_face_00 = (ω_face_00- exp(get_value(sm.start_step_props.gamma_turb[k]))) / sm.props.dt
    dgammadt_face_00 = (γ_face_00- get_value(sm.start_step_props.gamma_turb[k])) / sm.props.dt


    ## cell centre values required for mixing term ##
    ## Both 00 and p1 terms are needed for F_i and F_i+1

    ## 00 values ##
    # γ_cc_00 = 0.5*(get_00_dual(sm.props.gamma_turb[k]) + get_m1_dual(sm.props.gamma_turb[k-1])) 
    # ω_cc_00 = exp(γ_cc_00) ###
    ω_cc_00 = 0.5*(exp(get_00_dual(sm.props.gamma_turb[k])) + exp(get_m1_dual(sm.props.gamma_turb[k-1])))#exp(γ_cc_00) ###
    dm_face_00 = 0.5*(sm.props.dm[k] + sm.props.dm[k-1]) ##
    m_face_00 = 0.5*(sm.props.m[k] + sm.props.m[k-1])
    r_cc_00 = 0.5*(exp(get_00_dual(sm.props.lnr[k])) + exp(get_m1_dual(sm.props.lnr[k-1])))
    P_cc_00 = get_00_dual(sm.props.eos_res[k].P)
    ρ_cc_00 = get_00_dual(sm.props.eos_res[k].ρ)
    L_cc_00 = 0.5*(get_00_dual(sm.props.L[k]) + get_m1_dual(sm.props.L[k-1])) * LSUN
    Hₚ_cc_00 = P_cc_00 / (ρ_cc_00 * m_face_00 * CGRAV / r_cc_00^2)
    Λ_cc_00 = 1 / (1 / Hₚ_cc_00 + 1 / r_cc_00)
    A_00 = (4π  *ρ_cc_00 * r_cc_00^2)^2 * Λ_cc_00 * α_w_FACTOR[] * sqrt(ω_cc_00)
    F_00 = (A_00 / sm.props.dm[k]) * (exp(get_00_dual(sm.props.gamma_turb[k])) - exp(get_m1_dual(sm.props.gamma_turb[k-1])))

    ## P1 values ##
    P_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].P)
    ρ_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].ρ)
    T_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].T)
    dm_face_p1 = 0.5*(sm.props.dm[k] + sm.props.dm[k+1])
    m_face_p1 = 0.5*(sm.props.m[k] + sm.props.m[k+1]) ##
    L_cc_p1 = 0.5*(get_p1_dual(sm.props.L[k+1]) + get_00_dual(sm.props.L[k])) * LSUN
    r_cc_p1 = 0.5*(exp(get_p1_dual(sm.props.lnr[k+1])) + exp(get_00_dual(sm.props.lnr[k])))
    Hₚ_cc_p1 = P_cc_p1 / (ρ_cc_p1 * m_face_p1 * CGRAV / r_cc_p1^2)
    Λ_cc_p1 = 1 / (1 / Hₚ_cc_p1 + 1 / r_cc_p1)
    ω_cc_p1 = 0.5*(exp(get_p1_dual(sm.props.gamma_turb[k+1])) + exp(get_00_dual(sm.props.gamma_turb[k])))
    A_p1 = (4π *ρ_cc_p1* r_cc_p1^2)^2 * Λ_cc_p1 * α_w_FACTOR[] * sqrt(ω_cc_p1)
    F_p1 = (A_p1 / sm.props.dm[k+1]) * (exp(get_p1_dual(sm.props.gamma_turb[k+1])) - exp(get_00_dual(sm.props.gamma_turb[k])))

    # Calculation of all terms for residual 
    mixing_term = (F_p1 - F_00) / dm_face_p1
    omega_var_term = dgammadt_face_00 * ω_face_00 
    source_term = α₁_face_00 * SA_face_00 * sqrt(ω_face_00)
    turb_dissipation_term = C_d * (ω_face_00)^(3/2) / Λ_face_00
    rad_dissipation_term = ω_face_00 / τᵣ_face_00
    excess_term = (c_s_face_00)^2 * EXCESS_FACTOR[]/ τᵣ_face_00

     return ForwardDiff.value(τᵣ_face_00), ForwardDiff.value(rad_dissipation_term), ForwardDiff.value(turb_dissipation_term),  ForwardDiff.value(Λ_face_00), ForwardDiff.value(excess_term)
end 
end 
    

##
StellarModels.add_profile_option!(sm, "gamma_turb_energy", "unitless", (sm, k) -> ((get_value(sm.props.gamma_turb[k]))))
StellarModels.add_profile_option!(sm, "D_turb_kuff", "unitless", (sm, k) -> ((get_value(sm.props.D_turb[k]))))
StellarModels.add_profile_option!(sm, "kappa", "unitless", (sm, k) -> ((get_value(sm.props.κ[k]))))
StellarModels.add_profile_option!(sm, "nablaa_face", "unitless", (sm, k) -> get_value(sm.props.∇ₐ_face[k]), label="\nabla_\text{a,face}")
StellarModels.add_profile_option!(sm, "nablar_face", "unitless", (sm, k) -> get_value(sm.props.turb_res[k].∇ᵣ),  label="\nabla_\text{r,face}")
StellarModels.add_profile_option!(sm, "nabla_face_k", "unitless", (sm, k) -> get_value(sm.props.turb_res[k].∇), label="\nabla_\text{face}")
StellarModels.add_profile_option!(sm, "v_turb", "unitless", (sm, k) -> get_value(sm.props.turb_res[k].v_turb), label="\nabla_\text{r,face}")
StellarModels.add_profile_option!(sm,"tau", "unitless", (sm,k) -> calculate_tau_CD_by_lambda(sm,k)[1], label="\nabla_\text{r,face}")
StellarModels.add_profile_option!(sm,"rad_diss", "unitless", (sm,k) -> calculate_tau_CD_by_lambda(sm,k)[2], label="\nabla_\text{r,face}")
StellarModels.add_profile_option!(sm,"vis_diss", "unitless", (sm,k) -> calculate_tau_CD_by_lambda(sm,k)[3], label="\nabla_\text{r,face}")
StellarModels.add_profile_option!(sm,"Lambda", "unitless", (sm,k) -> calculate_tau_CD_by_lambda(sm,k)[4], label="\nabla_\text{r,face}")
StellarModels.add_profile_option!(sm,"excess_term", "unitless", (sm,k) -> calculate_tau_CD_by_lambda(sm,k)[5], label="\nabla_\text{r,face}")
StellarModels.add_history_option!(sm, "alpha_overshoot", "H_p", sm ->calculate_overshoot_length(sm)[1], label=L"\text{Alpha_ov}\,[\alpha_{ov}]")
StellarModels.add_history_option!(sm, "Sch_radius", "Rsun", sm ->calculate_overshoot_length(sm)[2], label=L"\text{Alpha_ov}\,[\alpha_{ov}]")
StellarModels.add_history_option!(sm, "pressure_scale_height", "unitless", sm ->calculate_overshoot_length(sm)[3], label=L"\text{Alpha_ov}\,[\alpha_{ov}]")
StellarModels.add_history_option!(sm, "ov_radius", "unitless", sm ->calculate_overshoot_length(sm)[4], label=L"\text{Alpha_ov}\,[\alpha_{ov}]")
StellarModels.add_history_option!(sm, "ov_distance", "unitless", sm ->calculate_overshoot_length(sm)[5], label=L"\text{Alpha_ov}\,[\alpha_{ov}]")
# "tau", "rad_diss", "vis_diss", "Lambda", "excess_term"
open("example_options.toml", "w") do file
    write(file,
          """
          [remesh]
          do_remesh = true

          [solver]
          newton_max_iter_first_step = 1000
          initial_model_scale_max_correction = 0.2
          newton_max_iter = 400
          scale_max_correction = 1.0
          solver_progress_iter = 1
          relative_correction_tolerance = 1e15
          maximum_residual_tolerance = 1e-2
          use_preconditioning = true

          [timestep]
          dt_max_increase = 1.5
          delta_R_limit = 0.01
          delta_Tc_limit = 0.01
          delta_Xc_limit = 0.005

          [termination]
          max_model_number = 20000
          max_center_T = 1e12

          [io]
          profile_interval = 1
          terminal_header_interval = 100
          terminal_info_interval = 100
          profile_values = ["zone", "mass", "dm", "log10_rho", "log10_r", "log10_P", "log10_T", "luminosity",
                                      "X", "Y", "kappa", "nablaa_face", "nablar_face","v_turb", "gamma_turb_energy","D_turb_kuff","nabla_face_k"]
          history_values = ["model_number", "age", "dt", "star_mass", "X_center", "alpha_overshoot","pressure_scale_height", "Sch_radius", "ov_radius", "ov_distance"]
          """)
end
StellarModels.set_options!(sm.opt, "./example_options.toml")
rm(sm.opt.io.hdf5_history_filename; force=true)
rm(sm.opt.io.hdf5_profile_filename; force=true)

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
#set initial condition and run model
n = 1.5
# StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], 
#                                             1 * MSUN, 100 * RSUN; initial_dt=10 * SECYEAR)
tdc_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], 
                                           5 * MSUN, 80 *  RSUN; initial_dt=0.01 * SECYEAR)

           
@time Evolution.do_evolution_loop!(sm, plotter=plotter); 
  
##

# check for NaNs/Infs in the Jacobian 
for k in 1:1000
    if !all(isfinite, sm.solver_data.jacobian_D[k]) 
        println("CRITICAL: Non-finite value (NaN/Inf) found at zone k = $k")
        println("A[k] block: ", sm.solver_data.jacobian_D[k])
end

end