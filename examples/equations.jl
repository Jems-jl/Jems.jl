using Jems
import Jems.StellarModels: CellDualData, FaceDualData
import Jems.EOS: EOSResults
import Jems.Turbulence: TurbResults
import Jems.Evolution
import Jems.StellarModels: update_cell_dual_data_value!, 
                           update_cell_dual_data!, 
                           update_face_dual_data!, 
                           get_cell_dual, 
                           get_face_dual, 
                           get_face_00_dual, 
                           get_face_p1_dual,
                           get_00_dual,
                           get_m1_dual,
                           get_p1_dual,
                           eval_face_property!,
                           eval_face_property_log!,
                           update_struct_cell_dual_data,
                           update_struct_face_dual_data,get_value
using ForwardDiff

struct TDCEquationSet<:Jems.StellarModels.AbstractEquationSet
end 


"""
Custom StellarModelProperties
"""

@kwdef mutable struct TDCStellarModelProperties{TN, TDual, TDualFace, TCellDualData, TFaceDualData} <: Jems.StellarModels.AbstractModelProperties
    # scalar quantities
    dt::TN  # Timestep of the current evolutionary step (s)
    dt_next::TN
    time::TN  # Age of the model (s)
    model_number::Int
    mstar::TN  # Total model mass (g)

    # array of the values of the independent variables, everything should be reconstructable from this (mesh dependent)
    ind_vars::Vector{TN}
    nz::Int
    m::Vector{TN}
    dm::Vector{TN}

    eos_res_dual::Vector{EOSResults{TDual}}
    eos_res::Vector{EOSResults{TCellDualData}}

    # independent variables (duals constructed from the ind_vars array)
    # represents a staggered mesh: T, ρ and abundances are defined in the center of each cell, L and r on the outer face
    lnT::Vector{TCellDualData}  # [K]
    lnρ::Vector{TCellDualData}  # [g cm^-3]
    lnr::Vector{TCellDualData}  # [cm]
    L::Vector{TCellDualData}    # Lsun
    gamma_turb::Vector{TCellDualData}  # turbulent energy
    xa::Matrix{TCellDualData}   # dim-less
    xa_dual::Matrix{TDual}      # only the cell duals wrt itself

    # opacity (cell centered)
    κ::Vector{TCellDualData}  # cm^2 g^-1

    # rates (cell centered)
    rates::Matrix{TCellDualData}  # g^-1 s^-1
    rates_dual::Matrix{TDual}     # only cell duals wrt itself

    # face values
    lnP_face::Vector{TFaceDualData}  # [dyne]
    lnT_face::Vector{TFaceDualData}  # [K]
    lnρ_face::Vector{TFaceDualData}  # [g cm^-3]
    κ_face::Vector{TFaceDualData}    # cm^2 g^-1
    ∇ₐ_face::Vector{TFaceDualData}   # dim-less
    ∇ᵣ_face::Vector{TFaceDualData}   # dim-less 
    δ_face::Vector{TFaceDualData}   # dim-less
    cₚ_face::Vector{TFaceDualData}   # erg K^-1 g^-1

    #D_turb 
    D_turb::Vector{TFaceDualData}  # cm^2 s^-1
    gamma_turb_face::Vector{TFaceDualData}
    # turbulence (i.e. convection, face valued)
    turb_res_dual::Vector{TurbResults{TDualFace}}
    turb_res::Vector{TurbResults{TFaceDualData}}

    # flux term for mixing equations (4πr^2ρ)^2 D / dm
    flux_term::Vector{TFaceDualData}

    ϵ_nuc::Vector{TN}

    mixing_type::Vector{Symbol}
end

function TDCStellarModelProperties(nvars::Int, nz::Int, nextra::Int, nrates::Int, nspecies::Int, vari::Dict{Symbol,Int},
                                ::Type{TN}) where {TN<:Real}

    # define the types
    CDDTYPE = CellDualData{nvars + 1,3 * nvars + 1,TN}  # full dual arrays
    FDDTYPE = FaceDualData{2 * nvars + 1,3 * nvars + 1,TN}
    TD = typeof(ForwardDiff.Dual(zero(TN), (zeros(TN, nvars))...))  # only the cell duals
    TDF = typeof(ForwardDiff.Dual(zero(TN), (zeros(TN, 2 * nvars))...))  # only the face duals

    # create the vector containing the independent variables
    ind_vars = zeros(TN, nvars * (nz + nextra))

    # result containers
    eos_res_dual = [EOSResults{TD}() for i = 1:(nz + nextra)]
    eos_res = [EOSResults{CDDTYPE}() for i = 1:(nz + nextra)]

    turb_res_dual = [TurbResults{TDF}() for i = 1:(nz + nextra)]
    turb_res = [TurbResults{FDDTYPE}() for i = 1:(nz + nextra)]

    # unpacked ind_vars
    lnT = [CellDualData(nvars, TN; is_ind_var=true, ind_var_i=vari[:lnT]) for i in 1:(nz+nextra)]
    lnρ = [CellDualData(nvars, TN; is_ind_var=true, ind_var_i=vari[:lnρ]) for i in 1:(nz+nextra)]
    lnr = [CellDualData(nvars, TN; is_ind_var=true, ind_var_i=vari[:lnr]) for i in 1:(nz+nextra)]
    L = [CellDualData(nvars, TN; is_ind_var=true, ind_var_i=vari[:lum]) for i in 1:(nz+nextra)]
    gamma_turb = [CellDualData(nvars, TN; is_ind_var=true, ind_var_i=vari[:gamma_turb]) for i in 1:(nz+nextra)]
    xa = Matrix{CDDTYPE}(undef,nz+nextra, nspecies)
    for k in 1:(nz+nextra)
        for i in 1:nspecies
            xa[k,i] = CellDualData(nvars, TN;
                        is_ind_var=true, ind_var_i=nvars-nspecies+i) # 4 in here is the number of non-composition variables being solved
        end
    end

    xa_dual = zeros(TD, nz + nextra, nspecies)
    rates_dual = zeros(TD, nz + nextra, nrates)

    # mesh
    m = zeros(TN, nz + nextra)
    dm = zeros(TN, nz + nextra)

    # for some reason using zeros just creates a bunch of instances of the same object
    # so we just initialize a vector of undef
    lnP_face = Vector{FDDTYPE}(undef, nz+nextra)#zeros(FDDTYPE, nz+nextra)
    lnρ_face = Vector{FDDTYPE}(undef, nz+nextra)#zeros(FDDTYPE, nz+nextra)
    lnT_face = Vector{FDDTYPE}(undef, nz+nextra)#zeros(FDDTYPE, nz+nextra)
    κ_face = Vector{FDDTYPE}(undef, nz+nextra)#zeros(FDDTYPE, nz+nextra)
    ∇ₐ_face = Vector{FDDTYPE}(undef, nz+nextra)#zeros(FDDTYPE, nz+nextra)
    ∇ᵣ_face = Vector{FDDTYPE}(undef, nz+nextra)#zeros(FDDTYPE, nz+nextra)
    δ_face = Vector{FDDTYPE}(undef, nz+nextra)#zeros(FDDTYPE, nz+nextra)
    cₚ_face = Vector{FDDTYPE}(undef, nz+nextra)#zeros(FDDTYPE, nz+nextra)
    κ = Vector{CDDTYPE}(undef, nz+nextra)  # zeros(CDDTYPE, nz+nextra)
    D_turb = Vector{FDDTYPE}(undef, nz+nextra) #zeros(FDDTYPE, nz+nextra)
    flux_term = Vector{FDDTYPE}(undef, nz+nextra)#zeros(FDDTYPE, nz+nextra)
    mixing_type::Vector{Symbol} = repeat([:no_mixing], nz+nextra)
    gamma_turb_face = Vector{FDDTYPE}(undef, nz+nextra)
    for k in 1:(nz+nextra)
        lnP_face[k] = FaceDualData(nvars, TN)
        lnρ_face[k] = FaceDualData(nvars, TN)
        lnT_face[k] = FaceDualData(nvars, TN)
        κ_face[k] = FaceDualData(nvars, TN)
        ∇ₐ_face[k] = FaceDualData(nvars, TN)
        ∇ᵣ_face[k] = FaceDualData(nvars, TN)
        δ_face[k] = FaceDualData(nvars, TN)
        cₚ_face[k] = FaceDualData(nvars, TN)
        κ[k] = CellDualData(nvars, TN)
        D_turb[k] = FaceDualData(nvars, TN)
        flux_term[k] = FaceDualData(nvars, TN)
        gamma_turb_face[k] = FaceDualData(nvars, TN)
    end

    rates_dual = zeros(TD, nz + nextra, nrates)
    rates = Matrix{CDDTYPE}(undef, nz + nextra, nrates)
    for k = 1:(nz + nextra)
        for i = 1:nrates
            rates[k, i] = CellDualData(nvars, TN)
        end
    end

    return TDCStellarModelProperties(; ind_vars=ind_vars, model_number=zero(Int),
                                  nz=nz, m=m, dm=dm, mstar=zero(TN),
                                  dt=zero(TN), dt_next=zero(TN), time=zero(TN),
                                  eos_res_dual=eos_res_dual,
                                  eos_res=eos_res,
                                  turb_res_dual=turb_res_dual,
                                  turb_res=turb_res,
                                  flux_term=flux_term,
                                  lnT=lnT,
                                  lnρ=lnρ,
                                  lnr=lnr,
                                  L=L,
                                  gamma_turb=gamma_turb,
                                  xa=xa,
                                  xa_dual=xa_dual,
                                  lnP_face=lnP_face,
                                  lnρ_face=lnρ_face,
                                  lnT_face=lnT_face,
                                  κ_face=κ_face,
                                  ∇ₐ_face=∇ₐ_face,
                                  ∇ᵣ_face=∇ᵣ_face,
                                  δ_face=δ_face,
                                  cₚ_face=cₚ_face,
                                  D_turb=D_turb,
                                  κ=κ,
                                  gamma_turb_face=gamma_turb_face,
                                  rates=rates,
                                  ϵ_nuc=zeros(nz + nextra),
                                  rates_dual=rates_dual,
                                  mixing_type=mixing_type)
end


function Jems.StellarModels.evaluate_stellar_model_properties!(sm, props::TDCStellarModelProperties)

    if sm.opacity isa BlendedOpacity
        ramp_steps = 1000.0  # Number of steps to reach 100% table opacity
        sm.opacity.λ = min(1.0, props.model_number / ramp_steps)

        println("Current Step: $(props.model_number) | Blending λ: $(sm.opacity.λ)")
    end

    lnT_i = sm.vari[:lnT]
    lnρ_i = sm.vari[:lnρ]
    lnr_i = sm.vari[:lnr]
    L_i = sm.vari[:lum]
    gamma_turb_i = sm.vari[:gamma_turb]

    Threads.@threads for i = 1:(props.nz)
        # update independent variables
        update_cell_dual_data_value!(props.lnT[i], props.ind_vars[(i-1)*(sm.nvars)+lnT_i])
        update_cell_dual_data_value!(props.lnρ[i], props.ind_vars[(i-1)*(sm.nvars)+lnρ_i])
        update_cell_dual_data_value!(props.lnr[i], props.ind_vars[(i-1)*(sm.nvars)+lnr_i])
        update_cell_dual_data_value!(props.L[i], props.ind_vars[(i-1)*(sm.nvars)+L_i])
        update_cell_dual_data_value!(props.gamma_turb[i], props.ind_vars[(i-1)*(sm.nvars)+gamma_turb_i])
        for j in 1:sm.network.nspecies
            update_cell_dual_data_value!(props.xa[i,j],
                            props.ind_vars[(i-1)*(sm.nvars)+(sm.nvars - sm.network.nspecies + j)])
            props.xa_dual[i,j] = get_cell_dual(props.xa[i,j])
        end

        lnT = get_cell_dual(props.lnT[i])
        lnρ = get_cell_dual(props.lnρ[i])
        xa = @view props.xa_dual[i, :]

        # evaluate EOS
        set_EOS_resultsTρ!(sm.eos, props.eos_res_dual[i], lnT, lnρ, xa, sm.network.species_names)
        update_struct_cell_dual_data(props.eos_res[i], props.eos_res_dual[i])

        # evaluate opacity
        κ_dual = get_opacity_resultsTρ(sm.opacity, lnT, lnρ, xa, sm.network.species_names)
        update_cell_dual_data!(props.κ[i], κ_dual)

        # evaluate rates
        rates = @view props.rates_dual[i, :]
        set_rates_for_network!(rates, sm.network, exp(lnT), exp(lnρ), xa)
        for j in eachindex(rates)
            update_cell_dual_data!(props.rates[i, j], rates[j])
        end

        # compute eps_nuc
        props.ϵ_nuc[i] = 0.0
        for j in eachindex(rates)
            props.ϵ_nuc[i] += rates[j].value * sm.network.reactions[j].Qvalue
        end
    end

    # do face values next
    Threads.@threads for i = 1:(props.nz - 1)
        Jems.StellarModels.eval_face_property_log!(props.κ[i], props.κ[i + 1], props.dm[i], props.dm[i+1], props.κ_face[i])
    Jems.StellarModels.eval_face_property!(props.eos_res[i].lnP, props.eos_res[i+1].lnP, props.dm[i], props.dm[i+1], props.lnP_face[i])
    Jems.StellarModels.eval_face_property!(props.eos_res[i].lnρ, props.eos_res[i+1].lnρ, props.dm[i], props.dm[i+1], props.lnρ_face[i])
    Jems.StellarModels.eval_face_property!(props.eos_res[i].lnT, props.eos_res[i+1].lnT, props.dm[i], props.dm[i+1], props.lnT_face[i])
    Jems.StellarModels.eval_face_property!(props.eos_res[i].∇ₐ, props.eos_res[i+1].∇ₐ, props.dm[i], props.dm[i+1], props.∇ₐ_face[i])
    Jems.StellarModels.eval_face_property!(props.eos_res[i].δ, props.eos_res[i+1].δ, props.dm[i], props.dm[i+1], props.δ_face[i])
    Jems.StellarModels.eval_face_property!(props.eos_res[i].cₚ, props.eos_res[i+1].cₚ, props.dm[i], props.dm[i+1], props.cₚ_face[i])
    Jems.StellarModels.eval_face_property!(props.gamma_turb[i], props.gamma_turb[i+1], props.dm[i], props.dm[i+1], props.gamma_turb_face[i])
        gamma_turb_dual = get_face_dual(props.gamma_turb_face[i])
        
        κ_face_dual = get_face_dual(props.κ_face[i])
        ρ_face_dual = exp(get_face_dual(props.lnρ_face[i]))
        T_face_dual = exp(get_face_dual(props.lnT_face[i]))
        P_face_dual = exp(get_face_dual(props.lnP_face[i]))
        ∇ₐ_face_dual = get_face_dual(props.∇ₐ_face[i])
        δ_face_dual = get_face_dual(props.δ_face[i])
        cₚ_face_dual = get_face_dual(props.cₚ_face[i])
        
        L₀_dual = get_face_00_dual(props.L[i]) * LSUN
        r_dual = exp(get_face_00_dual(props.lnr[i]))
        set_turb_results!(sm.turbulence, props.turb_res_dual[i],
                    κ_face_dual, L₀_dual, ρ_face_dual, P_face_dual, T_face_dual, r_dual,
                    δ_face_dual, cₚ_face_dual, ∇ₐ_face_dual, props.m[i])
        update_struct_face_dual_data(props.turb_res[i], props.turb_res_dual[i])

        D_turb_dual = (1/3) * sqrt(2 * exp(gamma_turb_dual)) *  1 / (1/(P_face_dual / (ρ_face_dual * CGRAV * sm.props.m[i]/ r_dual^2)) + 1/r_dual)
        update_face_dual_data!(props.D_turb[i], D_turb_dual)
        flux_term_dual = (4π*r_dual^2*ρ_face_dual)^2*D_turb_dual/
                            (0.5*(props.dm[i]+props.dm[i+1]))
        update_face_dual_data!(props.flux_term[i], flux_term_dual)

        if get_value(props.turb_res[i].∇) < get_value(props.turb_res[i].∇ᵣ)
            props.mixing_type[i] = :convection
        else
            props.mixing_type[i] = :no_mixing
        end
    end
end

function Jems.StellarModels.hydro_vars(equation_set::TDCEquationSet)pairs
    return [:lnρ, :lnT, :lnr, :lum, :gamma_turb]
end 

function Jems.StellarModels.hydro_vars_scaling(equation_set::TDCEquationSet)
    return [:log, :log, :log, :maxval,:log]
end
function split_omega(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew_low, varnew_up)
    # use same omega in both cells to preserve energy
     varnew_low[sm.vari[:gamma_turb]] = var_00[sm.vari[:gamma_turb]]
     varnew_up[sm.vari[:gamma_turb]] = var_00[sm.vari[:gamma_turb]]
end
function Jems.StellarModels.remesh_splitting(equation_set::TDCEquationSet, sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew_low, varnew_up)
    StellarModels.split_lnr_lnρ(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew_low, varnew_up)
    StellarModels.split_lum(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew_low, varnew_up)
    StellarModels.split_lnT(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew_low, varnew_up)
    StellarModels.split_xa(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew_low, varnew_up)
    split_omega(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew_low, varnew_up)
end




function Jems.StellarModels.build_properties_for_equation_set(equation_set::TDCEquationSet, nvars, nz, nextra, network, vari, number_type)
    TDCStellarModelProperties(nvars, nz, nextra,
                                   length(network.reactions), network.nspecies, vari, number_type)
end


"""
Gamma Turb equation
"""


function gammaTurb(sm::StellarModel, k::Int)

###Constants###
C_d = 8/3 * sqrt(2/3)
α_w = 0.25 


"""
mixing term = 1/ dm(i,face) [ A(i+1) * (ω(i+1)- ω(i))/ dm(i+1,cell) - A(i) * (ω(i)- ω(i-1))/ dm(i,cell)]
A(i) = (4πr^2(c,i))^2 * ρ(i,c) * Λ(i,c)^2 * √ω(i,c) { X(i,c) = 0.5(X(i,f) + X(i-1,f))}
A(1) = A(nz+1)= 0
"""

## Outer boundary condition (k = 1)
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
    A_p1 = (4π  *ρ_cc_p1* r_cc_p1^2)^2 * Λ_cc_p1 * α_w * sqrt(ω_cc_p1)

    F_p1 = (A_p1 / dm_cell_p1) * (exp(get_p1_dual(sm.props.gamma_turb[k+1])) - exp(get_00_dual(sm.props.gamma_turb[k]))) 

    # Different terms for residual at k = 1
    mixing_term =  (F_p1/  dm_face_p1)  
    omega_var_term = (γ_face_00- get_value(sm.start_step_props.gamma_turb[k])) / sm.props.dt
    source_term = α₁_face_00 * SA_face_00 *sqrt(ω_face_00)
    turb_dissipation_term = C_d * (ω_face_00)^(3/2) / Λ_face_00
    rad_dissipation_term = ω_face_00 / τᵣ_face_00
    excess_term = C_d * (c_s_face_00 * 1e-4)^3 / Λ_face_00

    return  omega_var_term  + turb_dissipation_term + rad_dissipation_term - excess_term -  source_term
end



## Outer boundary condition (k = sm.props.nz)
if k == sm.props.nz

    # ==============================================================================
    # 1. THERMODYNAMICS & GEOMETRY (From EOS Results = Face Values)
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
    omega_var_term = dgammadt_face_00 * ω_face_00
    
    source_term = α₁_face_00 * SA_face_00 * sqrt(ω_face_00)
    turb_dissipation_term = C_d * (ω_face_00)^(3/2) / Λ_face_00
    rad_dissipation_term = ω_face_00 / τᵣ_face_00
    excess_term = C_d * (c_s_face_00 * 1e-4)^3 / Λ_face_00
    return  omega_var_term  + turb_dissipation_term + rad_dissipation_term - excess_term - source_term
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
    omega_var_term = dgammadt_face_00 * ω_face_00 
    source_term = α₁_face_00 * SA_face_00 * sqrt(ω_face_00)
    turb_dissipation_term = C_d * (ω_face_00)^(3/2) / Λ_face_00
    rad_dissipation_term = ω_face_00 / τᵣ_face_00
    excess_term = C_d * (c_s_face_00 * 1e-4)^3 / Λ_face_00  

    return omega_var_term  + turb_dissipation_term + rad_dissipation_term - excess_term - source_term
end 
end

function equationTDC_temp(sm::StellarModel, k::Int)
    lnT₀ = get_00_dual(sm.props.eos_res[k].lnT)
    if k == sm.props.nz  # atmosphere boundary condition
        L₀ = get_00_dual(sm.props.L[k]) * LSUN
        r₀ = exp(get_00_dual(sm.props.lnr[k]))
        return lnT₀ - log(L₀ / (SIGMA_SB * 4π * r₀^2)) / 4  # Eddington gray, ignoring radiation pressure term
    end
    r₀ = exp(get_00_dual(sm.props.lnr[k]))
    lnT₀ = get_00_dual(sm.props.lnT[k])
    lnT₊ = get_p1_dual(sm.props.lnT[k+1])

    Pface = exp(get_00_dual(sm.props.lnP_face[k]))

    #∇ = get_00_dual(sm.props.turb_res[k].∇)

   # Calculating ∇
    L = get_00_dual(sm.props.L[k]) * LSUN
    γ₀ = get_00_dual(sm.props.gamma_turb[k])
    ω = exp(γ₀)
    # ω_face = 0.5 * (ω_c_k + exp(get_p1_dual(sm.props.gamma_turb[k+1])))
    ρ_face = exp(get_00_dual(sm.props.lnρ_face[k]))
    P_face = exp(get_00_dual(sm.props.lnP_face[k]))
    r₀ = exp(get_00_dual(sm.props.lnr[k]))
    T_face = exp(get_00_dual(sm.props.lnT_face[k]))
    ∇ₐ = get_00_dual(sm.props.∇ₐ_face[k])
    # ∇ᵣ = get_00_dual(sm.props.∇ᵣ_face[k])
    cₚ =  get_00_dual(sm.props.cₚ_face[k])
    κ = get_00_dual(sm.props.κ_face[k])
    m₀ = sm.props.m[k]
    Hₚ = P_face / (ρ_face * CGRAV * m₀ / r₀^2) #defined at face 
    Λ = 1/(1/Hₚ + 1/r₀) 
    k_rad = 16 * SIGMA_SB * T_face^3 / (3 * κ * ρ_face)
    α₂ = ρ_face*cₚ*0.5*sqrt(2/3)*Λ*sqrt(ω)
    ∇ᵣ = 3 * κ * L * P_face / (16π * CRAD * CLIGHT * CGRAV * m₀ * T_face^4)
    SA = (∇ᵣ - ∇ₐ)*(1 + α₂/k_rad)^(-1)
    ∇ = ∇ₐ + SA 
       
    dm = 0.5*(sm.props.dm[k + 1] + sm.props.dm[k])
    
    return ((lnT₊ - lnT₀) / dm + CGRAV * sm.props.m[k] / (4π * r₀^4 * Pface) * ∇) /
        (CGRAV * sm.props.m[k] / (4π * r₀^4 * Pface))

end

function Jems.Evolution.eval_cell_eqs!(sm::StellarModel, ::TDCEquationSet, k::Int)
    sm.solver_data.eqs_duals[k, 1] = Evolution.equationHSE(sm, k)
    sm.solver_data.eqs_duals[k, 2] = equationTDC_temp(sm, k)
    sm.solver_data.eqs_duals[k, 3] = Evolution.equationContinuity(sm, k)
    sm.solver_data.eqs_duals[k, 4] = Evolution.equationLuminosity(sm, k)
    sm.solver_data.eqs_duals[k, 5] = gammaTurb(sm,k)
    # evaluate all composition equations
    for i = 1:(sm.network.nspecies)
        sm.solver_data.eqs_duals[k, sm.nvars - sm.network.nspecies + i] = Evolution.equation_composition(sm, k, sm.network.species_names[i])
    end
end
