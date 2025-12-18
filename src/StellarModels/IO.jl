using HDF5
using DataFrames
using Printf

export add_history_option, add_profile_option,
       history_output_units, history_output_functions, history_output_labels,
       profile_output_units, profile_output_functions, profile_output_labels

ioinit = false

const width = 9
const decimals = 4
const floatstr = "%#$width.$decimals" * "g "
const intstr = "%$width" * "i "

mutable struct TerminalHeader
    header::String
    linefmts::Vector{Printf.Format}
end

const terminal_header = TerminalHeader("", [])

function setup_header(sm::StellarModel)
    terminal_header.linefmts = Vector{String}(undef, 2)
    terminal_header.linefmts[1] = Printf.Format(intstr * floatstr^6 * intstr * "\n")
    terminal_header.linefmts[2] = Printf.Format(floatstr^7 * intstr * "\n")

    terminal_header.header = """
        model     logdt      logL   logTeff     logPs     logρs    H_cntr     iters
         mass       age      logR     logTc     logPc     logρc   He_cntr     zones
    -------------------------------------------------------------------------------
    """
end

function setup_header(oz::OneZone)
    lines = oz.network.nspecies ÷ 6 + 1
    lastline = oz.network.nspecies % 6
    if lastline == 0
        terminal_header.linefmts = Vector{Printf.Format}(undef, lines)
    else
        terminal_header.linefmts = Vector{Printf.Format}(undef, lines + 1)
    end
    terminal_header.linefmts[1] = Printf.Format(intstr * floatstr^4 * intstr * "\n")

    j = oz.network.nspecies
    i = 2
    while j > 6
        terminal_header.linefmts[i] = Printf.Format(floatstr^6 * "\n")
        j -= 6
        i += 1
    end
    terminal_header.linefmts[end] = Printf.Format(floatstr^j * "\n")

    terminal_header.header = """
        model     logdt       age      logT      logρ     iters
    """

    for (j, species) in enumerate(oz.network.species_names)
        if j % 6 != 1
            speciesstr = lpad(String(species), width+1)
        else
            speciesstr = lpad(String(species), width)
        end
        terminal_header.header *= speciesstr
        if j == oz.network.nspecies || j % 6 == 0
            terminal_header.header *= "\n"
        end
    end

    terminal_header.header *= """
    -----------------------------------------------------------
    """
end

function add_history_option!(m, name, unit, func; label::Union{LaTeXStrings.LaTeXString, String}="")
    if haskey(m.history_output_units, name)
        throw(ArgumentError("Key $name is already part of the history output options"))
    end
    m.history_output_units[name] = unit
    m.history_output_functions[name] = func
    if label == ""
        m.history_output_labels[name] = name
    else
        m.history_output_labels[name] = label
    end
end

#### ====== Equation for profiles of turbulent energy equation terms ======= #####

function Energy_cal(sm::StellarModel)

###Constants###
C_d = 8/3 * sqrt(2/3)
α_w = 0.25 


###Intiating summation terms###
total_omega_var_term = 0.0 
total_mixing_term = 0.0
total_source_term = 0.0
total_turb_dissipation_term = 0.0 
total_rad_dissipation_term = 0.0
total_excess_term = 0.0 
total_normalized_mixing = 0.0

total_E_turb_new = 0.0
total_E_turb_old = 0.0
for k in 1:sm.props.nz
   
## Inner boundary condition (k = 1)
if k == 1
    ### face Values required for calculating all the other terms except mixing term ###
    γ_face_00 = get_00_dual(sm.props.gamma_turb[k])
    ω_face_00 = exp(γ_face_00)  ###
    dm_cell_00 = sm.props.dm[k] 
    m_cell_00 = sm.props.m[k]
    r_face_00 = exp(get_00_dual(sm.props.lnr[k]))
    P_face_00 = get_00_dual(sm.props.eos_res[k].P)
    ρ_face_00 = get_00_dual(sm.props.eos_res[k].ρ)
    T_face_00 = get_00_dual(sm.props.eos_res[k].T)
    κ_face_00  = get_00_dual(sm.props.κ[k]) 
    ∇ₐ_face_00 = get_00_dual(sm.props.eos_res[k].∇ₐ)
    cₚ_face_00 = get_00_dual(sm.props.eos_res[k].cₚ)
    L_face_00  = get_00_dual(sm.props.L[k]) * LSUN

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
    A_p1 = (4π  *ρ_cc_p1* r_cc_p1^2)^2 * Λ_cc_p1 * α_w * sqrt(ω_cc_p1)

    F_p1 = (A_p1 / dm_cell_p1) * (exp(get_p1_dual(sm.props.gamma_turb[k+1])) - exp(get_00_dual(sm.props.gamma_turb[k]))) 

    # Different terms for residual at k = 1
    mixing_term =  (F_p1/  dm_face_p1)  
    omega_var_term = dgammadt_face_00 
    source_term = α₁_face_00 * SA_face_00 *sqrt(ω_face_00)
    turb_dissipation_term = C_d * (ω_face_00)^(3/2) / Λ_face_00
    rad_dissipation_term = ω_face_00 / τᵣ_face_00
    excess_term = C_d * (c_s_face_00 * 1e-7)^3 / Λ_face_00
    normalized_mixing = mixing_term/ (abs(source_term) + abs(excess_term) + abs(turb_dissipation_term) + abs(rad_dissipation_term))
    


## Outer boundary condition (k = sm.props.nz)
elseif k == sm.props.nz

    # ==============================================================================
    # 1. THERMODYNAMICS & GEOMETRY (From EOS Results = Face Values)
    # ==============================================================================
    # As per instruction: EOS results here are defined at the face
    γ_face_00 = get_00_dual(sm.props.gamma_turb[k])
    ω_face_00 = exp(γ_face_00)
    ω_cc_00 = 0.5*(exp(get_00_dual(sm.props.gamma_turb[k])) + exp(get_m1_dual(sm.props.gamma_turb[k-1])))
    # Geometry
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
    
    dgammadt_face_00 = (ω_face_00 - exp(get_value(sm.start_step_props.gamma_turb[k]))) / sm.props.dt
    
    # ==============================================================================
    # 3. FLUX CALCULATION (A_00)
    # ==============================================================================
    # We use the same Face values for A_00 as they are the definitive properties at k
    A_00 = (4π * ρ_face_00 * r_cc_00^2)^2 * Λ_face_00 * α_w * sqrt(ω_cc_00)
    
    # Flux entering from k-1
    F_00 = (A_00 / sm.props.dm[k]) * (exp(get_00_dual(sm.props.gamma_turb[k])) - exp(get_m1_dual(sm.props.gamma_turb[k-1])))

    # ==============================================================================
    # 4. RESIDUAL
    # ==============================================================================
    mixing_term = -(F_00 / sm.props.dm[k]) # Flux out (F_p1) is zero at surface
    omega_var_term = dgammadt_face_00 
    
    source_term = α₁_face_00 * SA_face_00 * sqrt(ω_face_00)
    turb_dissipation_term = C_d * (ω_face_00)^(3/2) / Λ_face_00
    rad_dissipation_term = ω_face_00 / τᵣ_face_00
    excess_term = C_d * (c_s_face_00 * 1e-7)^3 / Λ_face_00
    normalized_mixing = mixing_term/ (abs(source_term) + abs(excess_term) + abs(turb_dissipation_term) + abs(rad_dissipation_term))


### Other Calculations : 1 < k < nz ###
else

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
    dgammadt_face_00 = (ω_face_00- exp(get_value(sm.start_step_props.gamma_turb[k]))) / sm.props.dt
    


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
    A_00 = (4π  *ρ_cc_00 * r_cc_00^2)^2 * Λ_cc_00 * α_w * sqrt(ω_cc_00)
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
    A_p1 = (4π *ρ_cc_p1* r_cc_p1^2)^2 * Λ_cc_p1 * α_w * sqrt(ω_cc_p1)
    F_p1 = (A_p1 / sm.props.dm[k+1]) * (exp(get_p1_dual(sm.props.gamma_turb[k+1])) - exp(get_00_dual(sm.props.gamma_turb[k])))

    # Calculation of all terms for residual 
    mixing_term = (F_p1 - F_00) / dm_face_p1
    omega_var_term = dgammadt_face_00  
    source_term = α₁_face_00 * SA_face_00 * sqrt(ω_face_00)
    turb_dissipation_term = C_d * (ω_face_00)^(3/2) / Λ_face_00
    rad_dissipation_term = ω_face_00 / τᵣ_face_00
    excess_term = C_d * (c_s_face_00 * 1e-7)^3 / Λ_face_00  
    normalized_mixing = mixing_term/ (abs(source_term) + abs(excess_term) + abs(turb_dissipation_term) + abs(rad_dissipation_term))
    end 
integration_dm = 0.0
        if k == sm.props.nz
            integration_dm = sm.props.dm[k]
        else
            integration_dm = 0.5 * (sm.props.dm[k] + sm.props.dm[k+1])
        end  

total_source_term += source_term * integration_dm
total_mixing_term += mixing_term * integration_dm
total_omega_var_term += omega_var_term * integration_dm
total_turb_dissipation_term += turb_dissipation_term * integration_dm
total_rad_dissipation_term += rad_dissipation_term * integration_dm
total_excess_term += excess_term * integration_dm 
# total_normalized_mixing += normalized_mixing * integration_dm

ω_new = ω_face_00 
ω_old = exp(get_value(sm.start_step_props.gamma_turb[k]))
        
total_E_turb_new += ω_new * integration_dm
total_E_turb_old += ω_old * integration_dm
end 

delta_E_actual = total_E_turb_new - total_E_turb_old
delta_E_actual = total_omega_var_term * sm.props.dt
total_physics_rate = total_source_term - total_turb_dissipation_term - total_rad_dissipation_term + total_excess_term + total_mixing_term
delta_E_expected = total_physics_rate * sm.props.dt

abs_error = delta_E_actual - delta_E_expected
relative_error = abs_error / (total_E_turb_new + 1e-50)
total_normalized_mixing = total_mixing_term * sm.props.dt/total_E_turb_new
    return return (
    ForwardDiff.value(total_source_term),
    ForwardDiff.value(total_mixing_term),
    ForwardDiff.value(total_omega_var_term),
    ForwardDiff.value(total_turb_dissipation_term),
    ForwardDiff.value(total_rad_dissipation_term),
    ForwardDiff.value(total_excess_term),
    ForwardDiff.value(total_normalized_mixing), 
    ForwardDiff.value(abs_error),
    ForwardDiff.value(relative_error)
)
end 



#### ====== End of equation ======= #####
"""
Beginning of function for setting up alpha_overshoot history functions
"""

function interpolate_boundary_live(sm::StellarModel, i_rad::Int)
    i_conv = i_rad - 1
    
    #Radius 
    logr_rad  = get_value(sm.props.lnr[i_rad])
    logr_conv = get_value(sm.props.lnr[i_conv])
    #Pressure
    logP_rad  = get_value(sm.props.eos_res[i_rad].P)
    logP_conv = get_value(sm.props.eos_res[i_conv].P)
    #density 
    logρ_rad  = get_value(sm.props.lnρ[i_rad])
    logρ_conv = get_value(sm.props.lnρ[i_conv])
    #mass
    m_rad  = sm.props.m[i_rad] 
    m_conv = sm.props.m[i_conv]

    # Delta-Nabla (Δ∇ = ∇_rad - ∇_ad)
    dnabla_rad  = get_value(sm.props.turb_res[i_rad].∇ᵣ) - get_value(sm.props.∇ₐ_face[i_rad])
    dnabla_conv = get_value(sm.props.turb_res[i_conv].∇ᵣ) - get_value(sm.props.∇ₐ_face[i_conv])
    
    # 2. Linear Interpolation Calculation
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

    D_ov = get_value(sm.props.D_turb[i_ov])
    D_d = get_value(sm.props.D_turb[i_d])

    diff_D = D_ov - D_d
    diff_r = logr_ov - logr_d
    slope = diff_D / diff_r

    logr_boundary = logr_d + (1e5 - D_d) / slope
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
    P_boundary = logP_sch
    ρ_boundary = exp(logρ_sch)
    r_boundary = exp(logr_sch) # Radius in cm
    m_boundary = m_sch
    
    HP_Sch = P_boundary / (ρ_boundary * CGRAV * m_boundary / r_boundary^2)

    
    i_ov = nothing

    #todo : calculate overshoot length using interpolation 
    for k = i_Sch + 1:nz_interior
        D_turb = get_value(sm.props.D_turb[k])
        if D_turb < OVERSHOOT_THRESHOLD
            i_ov = k
            break
        end
    end
    
    if isnothing(i_ov)
        # If mixing doesn't fall below threshold before the surface
        return 0.0 
    end
    logr_ov = interpolate_ov_boundary(sm, i_ov) #Redundant calculations 
    r_overshoot = exp(logr_ov) # Radius in cm 
   
    overshooting_distance_cm = abs(r_overshoot - r_boundary)
    alpha_ov = overshooting_distance_cm / HP_Sch
    
    return alpha_ov 
end

function get_overshoot(sm)
    # Calls the calculation function and extracts the pure numerical value.
    return calculate_overshoot_length(sm)
end
"""
End of function 
"""

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
function setup_model_history_functions!(sm::StellarModel)
    # general properties
    add_history_option!(sm, "age", "year", sm -> sm.props.time / SECYEAR, label=L"\text{age}\,[\text{yr}]")
    add_history_option!(sm, "dt", "year", sm -> sm.props.dt / SECYEAR, label=L"\Delta t\,[\text{yr}]")
    add_history_option!(sm, "model_number", "unitless", sm -> sm.props.model_number, label=L"\text{Model Number}")
    add_history_option!(sm, "star_mass", "Msun", sm -> sm.props.mstar / MSUN, label=L"\text{Mass}\,[M_\odot]")
    add_history_option!(sm, "alpha_overshoot", "H_p", get_overshoot, label=L"\text{Alpha_ov}\,[\alpha_{ov}]")
    # surface properties
    add_history_option!(sm, "R_surf", "Rsun", sm -> exp(get_value(sm.props.lnr[sm.props.nz])) / RSUN, label=L"R_\text{surf}\,[R_\odot]")
    add_history_option!(sm, "L_surf", "Lsun", sm -> get_value(sm.props.L[sm.props.nz]), label=L"\text{surf}\,[L_\odot]")
    add_history_option!(sm, "T_surf", "K", sm -> exp(get_value(sm.props.lnT[sm.props.nz])), label=L"T_\text{surf}\,[\text{K}]")
    add_history_option!(sm, "rho_surf", "g*cm^-3", sm -> exp(get_value(sm.props.lnρ[sm.props.nz])), label=L"\rho_\text{surf}\,[\text{g\,cm^{-3}}]")
    add_history_option!(sm, "P_surf", "dyne", sm -> exp(get_value(sm.props.eos_res[sm.props.nz].P)), label=L"P_\text{surf}\,[\text{dyne}]")
    add_history_option!(sm, "X_surf", "unitless", sm -> get_value(sm.props.xa[sm.props.nz, sm.network.xa_index[:H1]]), label=L"X_\text{surf}")
    add_history_option!(sm, "Y_surf", "unitless", sm -> get_value(sm.props.xa[sm.props.nz, sm.network.xa_index[:He4]]), label=L"Y_\text{surf}")

    # central properties
    add_history_option!(sm, "T_center", "K", sm -> exp(get_value(sm.props.lnT[1])), label=L"T_\text{c}\,[\text{K}]")
    add_history_option!(sm, "rho_center", "g*cm^-3", sm -> exp(get_value(sm.props.lnρ[1])), label=L"\rho_\text{c}\,[\text{g\,cm^{-3}}]")
    add_history_option!(sm, "P_center", "dyne", sm -> get_value(sm.props.eos_res[1].P), label=L"P_\text{c}\,[\text{dyne}]")
    add_history_option!(sm, "X_center", "unitless", sm -> get_value(sm.props.xa[1, sm.network.xa_index[:H1]]), label=L"X_\text{c}")
    add_history_option!(sm, "Y_center", "unitless", sm -> get_value(sm.props.xa[1, sm.network.xa_index[:He4]]), label=L"Y_\text{c}")

    #Turb History options 
    add_history_option!(sm, "budget_source", "erg/s", sm -> Energy_cal(sm)[1], label=L"E_\text{src}")
    add_history_option!(sm, "budget_mixing", "erg/s", sm -> Energy_cal(sm)[2], label=L"E_\text{mix}")
    add_history_option!(sm, "budget_omega_var", "erg/s", sm -> Energy_cal(sm)[3], label=L"dE/dt")
    add_history_option!(sm, "budget_diss_turb", "erg/s", sm -> Energy_cal(sm)[4], label=L"E_\text{diss,turb}")
    add_history_option!(sm, "budget_diss_rad", "erg/s", sm -> Energy_cal(sm)[5], label=L"E_\text{diss,rad}")
    add_history_option!(sm, "budget_diss_excess", "erg/s", sm -> Energy_cal(sm)[6], label=L"E_\text{diss,ex}")
    add_history_option!(sm, "budget_abs_error", "erg", sm -> Energy_cal(sm)[8], label=L"E_\text{abs_err,ex}")
    add_history_option!(sm, "budget_rel_error", "unitless", sm -> Energy_cal(sm)[9], label=L"E_\text{rek_err,ex}")
    add_history_option!(sm,"budget_norm_mix", "unitless", sm -> Energy_cal(sm)[7], label = L"E_\text{Norm_mix}" )
    # add species
    for j in eachindex(sm.network.species_names)
        species = sm.network.species_names[j]
        add_history_option!(sm, String(species)*"_center", "unitless",
            sm -> get_value(sm.props.xa[1,sm.network.xa_index[species]]), label="$(String(species))_\\text{c}")
        add_history_option!(sm, String(species)*"_surf", "unitless",
            sm -> get_value(sm.props.xa[sm.props.nz,sm.network.xa_index[species]]), label="$(String(species))_\\text{s}")
    end
end


function setup_model_history_functions!(oz::OneZone)
    # general properties
    add_history_option!(oz, "age", "year", oz -> oz.props.time / SECYEAR, label=L"\text{Age}\,[\text{yr}]")
    add_history_option!(oz, "dt", "year", oz -> oz.props.dt / SECYEAR, label=L"\Delta t\,[\text{yr}]")
    add_history_option!(oz, "model_number", "unitless", oz -> oz.props.model_number, label="\text{Model Number}")

    add_history_option!(oz, "T", "K", oz -> oz.props.T, label=L"T\,[\text{K}]")
    add_history_option!(oz, "rho", "g*cm^-3", oz -> oz.props.ρ, label=L"\rho\,[\text{g\,cm^{-3}}]")
    for j in eachindex(oz.network.species_names)
        species = oz.network.species_names[j]
        add_history_option!(oz, String(species), "unitless", oz -> get_value(oz.props.xa[oz.network.xa_index[species]]))
    end
end

function add_profile_option!(m, name, unit, func; label::Union{LaTeXStrings.LaTeXString, String}="")
    if haskey(m.profile_output_units, name)
        throw(ArgumentError("Key $name is already part of the history output options"))
    end
    m.profile_output_units[name] = unit
    m.profile_output_functions[name] = func
    if label == ""
        m.profile_output_labels[name] = name
    else
        m.profile_output_labels[name] = label
    end
end

function setup_model_profile_functions!(sm::StellarModel)
    # general properties
    add_profile_option!(sm, "zone", "unitless", (sm, k) -> k, label=L"\text{Zone}")
    add_profile_option!(sm, "mass", "Msun", (sm, k) -> sm.props.m[k] / MSUN, label=L"\text{Mass}\,[M_\odot]")
    add_profile_option!(sm, "dm", "Msun", (sm, k) -> sm.props.dm[k] / MSUN, label=L"\Delta m\,[M_\odot]")

    # thermodynamic properties
    add_profile_option!(sm, "log10_r", "log10(Rsun)", (sm, k) -> get_value(sm.props.lnr[k]) * log10(ℯ) - log10(RSUN), label=L"\log_{10}(r/R_\odot)")
    add_profile_option!(sm, "log10_P", "log10(dyne)", (sm, k) -> log10(get_value(sm.props.eos_res[k].P)), label=L"\log_{10}(P/[\text{dyne}])")
    add_profile_option!(sm, "log10_T", "log10(K)", (sm, k) -> get_value(sm.props.lnT[k]) * log10(ℯ), label=L"\log_{10}(T/[\text{K}])")
    add_profile_option!(sm, "log10_rho", "log10_(g*cm^-3)", (sm, k) -> get_value(sm.props.lnρ[k]) * log10(ℯ), label=L"\log_{10}(\rho/[\text{g\,cm^{-3}}])")
    add_profile_option!(sm, "luminosity", "Lsun", (sm, k) -> get_value(sm.props.L[k]) / LSUN, label=L"L/L_\odot")

    # abundances
    add_profile_option!(sm, "X", "unitless", (sm, k) -> get_value(sm.props.xa[k, sm.network.xa_index[:H1]]), label=L"X")
    add_profile_option!(sm, "Y", "unitless", (sm, k) -> get_value(sm.props.xa[k, sm.network.xa_index[:He4]]), label=L"Y")

    # temperature gradients
    add_profile_option!(sm, "nabla_a_face", "unitless", (sm, k) -> get_value(sm.props.∇ₐ_face[k]), label=L"\nabla_\text{a,face}")
    add_profile_option!(sm, "nabla_r_face", "unitless", (sm, k) -> get_value(sm.props.turb_res[k].∇ᵣ),  label=L"\nabla_\text{r,face}")
    add_profile_option!(sm, "nabla_face", "unitless", (sm, k) -> get_value(sm.props.turb_res[k].∇), label=L"\nabla_\text{face}")
    add_profile_option!(sm, "D_face", "cm^2*s^{-1}", (sm, k) -> get_value(sm.props.turb_res[k].D_turb), label=L"D_\text{face}\,[\text{cm^2\,s^{-1}}]")

    #extras 
    add_profile_option!(sm, "velocity_turb", "unitless", (sm, k) -> sqrt(2*exp(get_value(sm.props.gamma_turb[k]))))
    add_profile_option!(sm, "turb_energy", "unitless", (sm, k) -> (exp(get_value(sm.props.gamma_turb[k]))))
    add_profile_option!(sm, "nabla_tdc", "unitless", get_nabla_tdc)
    add_profile_option!(sm, "D_face_kuhfuss", "unitless", (sm, k) -> get_value(sm.props.D_turb[k]))

    
end

function init_IO(m::AbstractModel)
    setup_header(m)
    setup_model_history_functions!(m)
    if isa(m, OneZone)
        global ioinit = true
        return
    end
    setup_model_profile_functions!(m)
    global ioinit = true
end

function clear_IO()
    terminal_header.header = ""
    terminal_header.linefmts = []
    global ioinit = false
end

"""
    create_output_files(sm::StellarModel)

Creates output files for history and profile data
"""
function create_output_files!(m::AbstractModel)
    # Create history file
    m.history_file = h5open(m.opt.io.hdf5_history_filename, "w")
    data_cols = m.opt.io.history_values
    ncols = length(data_cols)

    # verify validity of column names
    for i in eachindex(data_cols)
        if data_cols[i] ∉ keys(m.history_output_functions)
            throw(ArgumentError("Invalid name for history data column, :$(data_cols[i])"))
        end
    end

    # Create history dataset in HDF5 file
    # Dataset is created with size (0, ncols), we will add rows by using the HDF5.set_extent_dims function
    # the (-1, ncols) is used to define the maximum extent of the dataset, -1 indicates that it is unbound
    # in number of rows. The chunk size is used for compression. Smaller chunk sizes will result in worse
    # compression but faster writes.
    # The compression level can be anywhere between 0 and 9, 0 being no compression 9 being the highest.
    # Compression is lossless.
    history = create_dataset(m.history_file, "history", Float64, ((0, ncols), (-1, ncols)),
                             chunk=(m.opt.io.hdf5_history_chunk_size, ncols),
                             compress=m.opt.io.hdf5_history_compression_level)

    # next up, include the units for all quantities. No need to recheck columns.
    attrs(history)["column_units"] = [m.history_output_units[data_cols[i]] for i in eachindex(data_cols)]
    # Finally, place column names
    attrs(history)["column_names"] = [data_cols[i] for i in eachindex(data_cols)]
    if (!m.opt.io.hdf5_history_keep_open)
        close(m.history_file)
    end

    if isa(m, OneZone)
        return
    end

    # Create profile file
    m.profiles_file = h5open(m.opt.io.hdf5_profile_filename, "w")
    data_cols = m.opt.io.profile_values
    # verify validity of column names
    for i in eachindex(data_cols)
        if data_cols[i] ∉ keys(m.profile_output_functions)
            throw(ArgumentError("Invalid name for profile data column, :$(data_cols[i])"))
        end
    end
    if (!m.opt.io.hdf5_profile_keep_open)
        close(m.profiles_file)
    end
end

function shut_down_IO!(m)
    if (m.opt.io.hdf5_history_keep_open)
        close(m.history_file)
    end
    if (m.opt.io.hdf5_profile_keep_open)
        close(m.profiles_file)
    end
    if ioinit
        clear_IO()
    end
end

"""
    write_data(sm::StellarModel)

Saves data (history/profile) for the current model, as required by the settings in `sm.opt.io`.
"""
function write_data(m::AbstractModel)
    # do history
    if (m.opt.io.history_interval > 0)
        file_exists = isfile(m.opt.io.hdf5_history_filename)
        if !file_exists
            throw(ErrorException("History file does not exist at $(m.opt.io.hdf5_history_filename)"))
        end
        if (m.props.model_number % m.opt.io.history_interval == 0)
            if (!m.opt.io.hdf5_history_keep_open)
                m.history_file = h5open(m.opt.io.hdf5_history_filename, "r+")
            end
            data_cols = m.opt.io.history_values
            ncols = length(data_cols)

            # after being sure the header is there, print the data
            history = m.history_file["history"]
            HDF5.set_extent_dims(history, (size(history)[1] + 1, ncols))
            for i in eachindex(data_cols)
                history[end, i] = m.history_output_functions[data_cols[i]](m)
            end
            if (!m.opt.io.hdf5_history_keep_open)
                close(m.history_file)
            end
        end
    end

    if isa(m, OneZone)
        return
    end

    # do profile
    if (m.opt.io.profile_interval > 0)
        file_exists = isfile(m.opt.io.hdf5_profile_filename)
        if !file_exists  # create file if it doesn't exist yet
            throw(ErrorException("Profile file does not exist at $(m.opt.io.hdf5_profile_filename)"))
        end
        if (m.props.model_number % m.opt.io.profile_interval == 0)
            if (!m.opt.io.hdf5_profile_keep_open)
                m.profiles_file = h5open(m.opt.io.hdf5_profile_filename, "r+")
            end
            data_cols = m.opt.io.profile_values
            ncols = length(data_cols)
            # Save current profile
            profile = create_dataset(m.profiles_file,
                                     "$(lpad(m.props.model_number,m.opt.io.hdf5_profile_dataset_name_zero_padding,"0"))",
                                     Float64, ((m.props.nz, ncols), (m.props.nz, ncols));
                                     chunk=(m.opt.io.hdf5_profile_chunk_size, ncols),
                                     compress=m.opt.io.hdf5_profile_compression_level)

            # next up, include the units for all quantities. No need to recheck columns.
            attrs(profile)["column_units"] = [m.profile_output_units[data_cols[i]] for i in eachindex(data_cols)]
            # Place column names
            attrs(profile)["column_names"] = [data_cols[i] for i in eachindex(data_cols)]

            # store data
            for i in eachindex(data_cols), k = 1:(m.props.nz)
                profile[k, i] = m.profile_output_functions[data_cols[i]](m, k)
            end
            if (!m.opt.io.hdf5_profile_keep_open)
                close(m.profiles_file)
            end
        end
    end
end

function write_terminal_info(sm::StellarModel; now::Bool=false)
    if sm.props.model_number == 1 || sm.props.model_number % sm.opt.io.terminal_header_interval == 0 || now
        print(terminal_header.header)
    end
    if sm.props.model_number == 1 || sm.props.model_number % sm.opt.io.terminal_info_interval == 0 || now
        Printf.format(stdout, terminal_header.linefmts[1],
                      sm.props.model_number,
                      log10(sm.props.dt / SECYEAR),
                      log10(get_value(sm.props.L[sm.props.nz])),
                      log10(ℯ) * get_value(sm.props.lnT[sm.props.nz]),
                      log10(get_value(sm.props.eos_res[sm.props.nz].P)),
                      log10(ℯ) * get_value(sm.props.lnρ[sm.props.nz]),
                      get_value(sm.props.xa[1, sm.network.xa_index[:H1]]),
                      sm.solver_data.newton_iters)
        Printf.format(stdout, terminal_header.linefmts[2],
                      sm.props.mstar / MSUN,
                      sm.props.time / SECYEAR,
                      log10(ℯ) * get_value(sm.props.lnr[sm.props.nz]) - log10(RSUN),
                      log10(ℯ) * get_value(sm.props.lnT[1]),
                      log10(get_value(sm.props.eos_res[1].P)),
                      log10(ℯ) * get_value(sm.props.lnρ[1]),
                      get_value(sm.props.xa[1, sm.network.xa_index[:He4]]),
                      sm.props.nz)
        println()
    end
end

function write_terminal_info(oz::OneZone; now::Bool=false)
    if oz.props.model_number == 1 || oz.props.model_number % oz.opt.io.terminal_header_interval == 0 || now
        print(terminal_header.header)
    end
    if oz.props.model_number == 1 || oz.props.model_number % oz.opt.io.terminal_info_interval == 0 || now
        Printf.format(stdout, terminal_header.linefmts[1],
                      oz.props.model_number,
                      log10(oz.props.dt / SECYEAR),
                      oz.props.time / SECYEAR,
                      log10(oz.props.T),
                      log10(oz.props.ρ),
                      oz.solver_data.newton_iters)
        j = oz.network.nspecies
        i = 1
        while j > 6
            Printf.format(stdout, terminal_header.linefmts[i + 1], (get_value.(oz.props.xa[((i - 1) * 6 + 1):(6i)]))...)
            i += 1
            j -= 6
        end
        Printf.format(stdout, terminal_header.linefmts[end], (get_value.(oz.props.xa[((i - 1) * 6 + 1):end]))...)
        println()
    end
end

"""
    get_history_dataframe_from_hdf5(hdf5_filename)

Returns a DataFrame object built from an hdf5 file, named `hdf5_filename`.
"""
function get_history_dataframe_from_hdf5(hdf5_filename)
    h5open(hdf5_filename) do history_file
        return DataFrame(history_file["history"][:, :], attrs(history_file["history"])["column_names"])
    end
end

"""
    get_profile_names_from_hdf5(hdf5_filename)

Retruns the column names of the profile data contained in the hdf5 file `hdf5_filename`.
"""
function get_profile_names_from_hdf5(hdf5_filename)
    h5open(hdf5_filename) do profiles_file
        return keys(profiles_file)
    end
end

"""
    get_profile_dataframe_from_hdf5(hdf5_filename, profile_name)

Returns a DataFrame object built from an hdf5 file, named `hdf5_filename`, considering the column named `profile_name`
"""
function get_profile_dataframe_from_hdf5(hdf5_filename, profile_name)
    h5open(hdf5_filename) do profiles_file
        return DataFrame(profiles_file[profile_name][:, :], attrs(profiles_file[profile_name])["column_names"])
    end
end
