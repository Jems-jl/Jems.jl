using ForwardDiff
using Jems
import Jems.StellarModels
import Jems.EOS: EOSResults
import Jems.Turbulence: TurbResults
import Jems.Evolution

# Only import valid, existing functions from DualSupport:
using Jems.DualSupport: DualData, DualSupport, get_value, update_dual_data_local!, update_dual_data_mixed_full!, update_dual_data_mixed_m1!, update_dual_data_mixed_p1!, update_dual_data_value!, get_local_dual, get_p1_dual, get_00_dual, get_m1_dual

struct TDCEquationSet<:Jems.StellarModels.AbstractEquationSet
end 


"""
Custom StellarModelProperties
"""


@inline function eval_mixed_p1_property_log!(prop_00, prop_p1, dm_00, dm_p1, mixed_prop)
    val00 = get_00_dual(prop_00)
    valp1 = get_p1_dual(prop_p1)
    valmixed_dual = exp((dm_p1 * log(val00) + dm_00 * log(valp1)) / (dm_00 + dm_p1))
    update_dual_data_mixed_p1!(mixed_prop, valmixed_dual)
end

@inline function eval_mixed_p1_property!(prop_00, prop_p1, dm_00, dm_p1, mixed_prop)
    val00 = get_00_dual(prop_00)
    valp1 = get_p1_dual(prop_p1)
    valmixed_dual = (dm_p1 * val00 + dm_00 * valp1) / (dm_00 + dm_p1)
    update_dual_data_mixed_p1!(mixed_prop, valmixed_dual)
end
@inline function eval_mixed_p1_property_log_surface!(prop_00, prop_p1, dm_00, dm_p1, mixed_prop)
    val00 = get_00_dual(prop_00)
    valp1 = get_p1_dual(prop_p1)
    valmixed_dual = exp((dm_p1 * log(val00) + 1/2 * dm_00 * log(valp1)) / (1/2 * dm_00 + dm_p1))
    update_dual_data_mixed_p1!(mixed_prop, valmixed_dual)
end
@inline function eval_mixed_p1_property_surface!(prop_00, prop_p1, dm_00, dm_p1, mixed_prop)
    val00 = get_00_dual(prop_00)
    valp1 = get_p1_dual(prop_p1)
    valmixed_dual = (dm_p1 * val00 + 1/2 * dm_00 * valp1) / (1/2 * dm_00 + dm_p1)
    update_dual_data_mixed_p1!(mixed_prop, valmixed_dual)
end
@inline function eval_mixed_p1_property_log_core!(prop_00, prop_p1, dm_00, dm_p1, mixed_prop)
    val00 = get_00_dual(prop_00)
    valp1 = get_p1_dual(prop_p1)
    valmixed_dual = exp((1/2 * dm_p1 * log(val00) +  dm_00 * log(valp1)) / (dm_00 + 1/2 * dm_p1))
    update_dual_data_mixed_p1!(mixed_prop, valmixed_dual)
end
@inline function eval_mixed_p1_property_core!(prop_00, prop_p1, dm_00, dm_p1, mixed_prop)
    val00 = get_00_dual(prop_00)
    valp1 = get_p1_dual(prop_p1)
    valmixed_dual = (1/2 * dm_p1 * val00 +  dm_00 * valp1) / ( dm_00 + 1/2 * dm_p1)
    update_dual_data_mixed_p1!(mixed_prop, valmixed_dual)
end

@generated function update_struct_dual_data_local(obj::T1,obj_dual::T2) where{T1,T2}
    names = fieldnames(T1)
    lines::Vector{Expr} = [:(update_dual_data_local!(obj.$name, obj_dual.$name)) for name in names]
    return quote
        $(lines...)
    end
end

@generated function update_struct_dual_data_mixed_p1(obj::T1,obj_dual::T2) where{T1,T2}
    names = fieldnames(T1)
    lines::Vector{Expr} = [:(update_dual_data_mixed_p1!(obj.$name, obj_dual.$name)) for name in names]
    return quote
        $(lines...)
    end
end


@kwdef mutable struct TDCStellarModelProperties{TN,TLocalDual,TFullDual,TDualData} <: AbstractModelProperties
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

    eos_res_dual::Vector{EOSResults{TLocalDual}}
    eos_res::Vector{EOSResults{TDualData}}

    # independent variables (duals constructed from the ind_vars array)
    # represents a staggered mesh: T, ρ and abundances are defined in the center of each cell, L and r on the outer face
    lnT::Vector{TDualData}  # [K]
    lnρ::Vector{TDualData}  # [g cm^-3]
    lnr::Vector{TDualData}  # [cm]
    L::Vector{TDualData}    # Lsun
    gamma_turb::Vector{TDualData}
    xa::Matrix{TDualData}   # dim-less
    xa_dual::Matrix{TLocalDual}      # only the cell duals wrt itself

    # opacity (cell centered)
    κ::Vector{TDualData}  # cm^2 g^-1

    # rates (cell centered)
    rates::Matrix{TDualData}  # g^-1 s^-1
    rates_dual::Matrix{TLocalDual}     # only cell duals wrt itself

    # face values
    lnP_face::Vector{TDualData}  # [dyne]
    lnT_face::Vector{TDualData}  # [K]
    lnρ_face::Vector{TDualData}  # [g cm^-3]
    κ_face::Vector{TDualData}    # cm^2 g^-1
    ∇ₐ_face::Vector{TDualData}   # dim-less
    ∇ᵣ_face::Vector{TDualData}   # dim-less
    δ_face::Vector{TDualData}   # dim-less
    cₚ_face::Vector{TDualData}   # erg K^-1 g^-1
    turb_flux_cc::Vector{TDualData}
    #D_turb 
    D_turb::Vector{TDualData} 
    gamma_turb_face::Vector{TDualData}

    # turbulence (i.e. convection, face valued)
    turb_res_dual::Vector{TurbResults{TFullDual}}
    turb_res::Vector{TurbResults{TDualData}}

    # flux term for mixing equations (4πr^2ρ)^2 D / dm
    flux_term::Vector{TDualData}

    ϵ_nuc::Vector{TN}

    mixing_type::Vector{Symbol}
end


function TDCStellarModelProperties(nvars::Int, nz::Int, nextra::Int, nrates::Int, nspecies::Int, vari::Dict{Symbol,Int},
                                ::Type{TN}, ::Type{internal_tag}) where {TN<:Real,internal_tag<:ForwardDiff.Tag}

    # define the types
    DDTYPE = DualData{nvars + 1,3 * nvars + 1,TN,internal_tag}  # full dual arrays
    TDL = typeof(ForwardDiff.Dual{internal_tag}(zero(TN), (zeros(TN, nvars))...))  # only the local duals
    TDF = typeof(ForwardDiff.Dual{internal_tag}(zero(TN), (zeros(TN, 3 * nvars))...))  # only the full duals

    # create the vector containing the independent variables
    ind_vars = zeros(TN, nvars * (nz + nextra))

    # result containers
    eos_res_dual = [EOSResults{TDL}() for i = 1:(nz + nextra)]
    eos_res = [EOSResults{DDTYPE}() for i = 1:(nz + nextra)]

    turb_res_dual = [TurbResults{TDF}() for i = 1:(nz + nextra)]
    turb_res = [TurbResults{DDTYPE}() for i = 1:(nz + nextra)]

    # unpacked ind_vars
    lnT = [DualData(nvars, TN, internal_tag; is_ind_var=true, ind_var_i=vari[:lnT]) for i in 1:(nz+nextra)]
    lnρ = [DualData(nvars, TN, internal_tag; is_ind_var=true, ind_var_i=vari[:lnρ]) for i in 1:(nz+nextra)]
    lnr = [DualData(nvars, TN, internal_tag; is_ind_var=true, ind_var_i=vari[:lnr]) for i in 1:(nz+nextra)]
    L = [DualData(nvars, TN, internal_tag; is_ind_var=true, ind_var_i=vari[:lum]) for i in 1:(nz+nextra)]
    gamma_turb = [DualData(nvars, TN, internal_tag; is_ind_var=true, ind_var_i=vari[:gamma_turb]) for i in 1:(nz+nextra)]
    xa = Matrix{DDTYPE}(undef,nz+nextra, nspecies)
    for k in 1:(nz+nextra)
        for i in 1:nspecies
            xa[k,i] = DualData(nvars, TN, internal_tag;
                        is_ind_var=true, ind_var_i=nvars-nspecies+i) # nvars-nspecies in here is the number of non-composition variables being solved
        end
    end

    xa_dual = zeros(TDL, nz + nextra, nspecies)
    rates_dual = zeros(TDL, nz + nextra, nrates)

    # mesh
    m = zeros(TN, nz + nextra)
    dm = zeros(TN, nz + nextra)

    # for some reason using zeros just creates a bunch of instances of the same object
    # so we just initialize a vector of undef
    lnP_face = Vector{DDTYPE}(undef, nz+nextra)#zeros(DDTYPE, nz+nextra)
    lnρ_face = Vector{DDTYPE}(undef, nz+nextra)#zeros(DDTYPE, nz+nextra)
    lnT_face = Vector{DDTYPE}(undef, nz+nextra)#zeros(DDTYPE, nz+nextra)
    κ_face = Vector{DDTYPE}(undef, nz+nextra)#zeros(DDTYPE, nz+nextra)
    ∇ₐ_face = Vector{DDTYPE}(undef, nz+nextra)#zeros(DDTYPE, nz+nextra)
    ∇ᵣ_face = Vector{DDTYPE}(undef, nz+nextra)#zeros(DDTYPE, nz+nextra)
    δ_face = Vector{DDTYPE}(undef, nz+nextra)#zeros(DDTYPE, nz+nextra)
    cₚ_face = Vector{DDTYPE}(undef, nz+nextra)#zeros(DDTYPE, nz+nextra)
    turb_flux_cc = Vector{DDTYPE}(undef, nz+nextra)#zeros(DDTYPE, nz+nextra)
    κ = Vector{DDTYPE}(undef, nz+nextra)  # zeros(DDTYPE, nz+nextra)
    D_turb = Vector{DDTYPE}(undef, nz+nextra)  # zeros(DDTYPE, nz+nextra)
    flux_term = Vector{DDTYPE}(undef, nz+nextra)#zeros(DDTYPE, nz+nextra)
    mixing_type::Vector{Symbol} = repeat([:no_mixing], nz+nextra)
    gamma_turb_face = Vector{DDTYPE}(undef, nz+nextra)#zeros(DDTYPE, nz+nextra)
    for k in 1:(nz+nextra)
        lnP_face[k] = DualData(nvars, TN, internal_tag)
        lnρ_face[k] = DualData(nvars, TN, internal_tag)
        lnT_face[k] = DualData(nvars, TN, internal_tag)
        κ_face[k] = DualData(nvars, TN, internal_tag)
        ∇ₐ_face[k] = DualData(nvars, TN, internal_tag)
        ∇ᵣ_face[k] = DualData(nvars, TN, internal_tag)
        δ_face[k] = DualData(nvars, TN, internal_tag)
        cₚ_face[k] = DualData(nvars, TN, internal_tag)
        turb_flux_cc[k] = DualData(nvars, TN, internal_tag)
        κ[k] = DualData(nvars, TN, internal_tag)
        D_turb[k] = DualData(nvars, TN, internal_tag)
        flux_term[k] = DualData(nvars, TN, internal_tag)
        gamma_turb_face[k] = DualData(nvars, TN, internal_tag)
    end

    rates_dual = zeros(TDL, nz + nextra, nrates)
    rates = Matrix{DDTYPE}(undef, nz + nextra, nrates)
    for k = 1:(nz + nextra)
        for i = 1:nrates
            rates[k, i] = DualData(nvars, TN, internal_tag)
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
                                  turb_flux_cc = turb_flux_cc,
                                  D_turb=D_turb,
                                  κ=κ,
                                  gamma_turb_face=gamma_turb_face,
                                  rates=rates,
                                  ϵ_nuc=zeros(TN, nz + nextra),
                                  rates_dual=rates_dual,
                                  mixing_type=mixing_type)
end

@inline function smooth_abs(x, ϵ)
    s = sqrt(x*x + ϵ*ϵ)
    return (x*x) / (s + ϵ)
end

function Jems.StellarModels.evaluate_stellar_model_properties!(sm, props::TDCStellarModelProperties)
    lnT_i = sm.vari[:lnT]
    lnρ_i = sm.vari[:lnρ]
    lnr_i = sm.vari[:lnr]
    L_i = sm.vari[:lum]
    gamma_turb_i = sm.vari[:gamma_turb]
    Threads.@threads for i = 1:(props.nz)
        # update independent variables
        update_dual_data_value!(props.lnT[i], props.ind_vars[(i-1)*(sm.nvars)+lnT_i])
        update_dual_data_value!(props.lnρ[i], props.ind_vars[(i-1)*(sm.nvars)+lnρ_i])
        update_dual_data_value!(props.lnr[i], props.ind_vars[(i-1)*(sm.nvars)+lnr_i])
        update_dual_data_value!(props.L[i], props.ind_vars[(i-1)*(sm.nvars)+L_i])
        update_dual_data_value!(props.gamma_turb[i], props.ind_vars[(i-1)*(sm.nvars)+gamma_turb_i])
        for j in 1:sm.network.nspecies
            update_dual_data_value!(props.xa[i,j],
                            props.ind_vars[(i-1)*(sm.nvars)+(sm.nvars - sm.network.nspecies + j)])
            props.xa_dual[i,j] = get_local_dual(props.xa[i,j])
        end

        lnT = get_local_dual(props.lnT[i])
        lnρ = get_local_dual(props.lnρ[i])
        xa = @view props.xa_dual[i, :]

        # evaluate EOS
        set_EOS_resultsTρ!(sm.eos, props.eos_res_dual[i], lnT, lnρ, xa, sm.network.species_names)
        update_struct_dual_data_local(props.eos_res[i], props.eos_res_dual[i])

        # evaluate opacity
        κ_dual = get_opacity_resultsTρ(sm.opacity, lnT, lnρ, xa, sm.network.species_names)
        update_dual_data_local!(props.κ[i], κ_dual)

        # evaluate rates
        rates = @view props.rates_dual[i, :]
        set_rates_for_network!(rates, sm.network, exp(lnT), exp(lnρ), xa)
        for j in eachindex(rates)
            update_dual_data_local!(props.rates[i, j], rates[j])
        end

        # compute eps_nuc
        props.ϵ_nuc[i] = 0.0
        for j in eachindex(rates)
            props.ϵ_nuc[i] += rates[j].value * sm.network.reactions[j].Qvalue
        end
    end

    # do face values next
    Threads.@threads for i = 1:(props.nz - 1)
        if i == 1 
            eval_mixed_p1_property_log_core!(props.κ[i], props.κ[i + 1], props.dm[i], props.dm[i+1], props.κ_face[i])
            eval_mixed_p1_property_core!(props.eos_res[i].lnP, props.eos_res[i+1].lnP, props.dm[i], props.dm[i+1], props.lnP_face[i])
            eval_mixed_p1_property_core!(props.eos_res[i].lnρ, props.eos_res[i+1].lnρ, props.dm[i], props.dm[i+1], props.lnρ_face[i])
            eval_mixed_p1_property_core!(props.eos_res[i].lnT, props.eos_res[i+1].lnT, props.dm[i], props.dm[i+1], props.lnT_face[i])
            eval_mixed_p1_property_core!(props.eos_res[i].∇ₐ, props.eos_res[i+1].∇ₐ, props.dm[i], props.dm[i+1], props.∇ₐ_face[i])
            eval_mixed_p1_property_core!(props.eos_res[i].δ, props.eos_res[i+1].δ, props.dm[i], props.dm[i+1], props.δ_face[i])
            eval_mixed_p1_property_core!(props.eos_res[i].cₚ, props.eos_res[i+1].cₚ, props.dm[i], props.dm[i+1], props.cₚ_face[i])
        elseif i == props.nz - 1
            eval_mixed_p1_property_log_surface!(props.κ[i], props.κ[i + 1], props.dm[i], props.dm[i+1], props.κ_face[i])
            eval_mixed_p1_property_surface!(props.eos_res[i].lnP, props.eos_res[i+1].lnP, props.dm[i], props.dm[i+1], props.lnP_face[i])
            eval_mixed_p1_property_surface!(props.eos_res[i].lnρ, props.eos_res[i+1].lnρ, props.dm[i], props.dm[i+1], props.lnρ_face[i])
            eval_mixed_p1_property_surface!(props.eos_res[i].lnT, props.eos_res[i+1].lnT, props.dm[i], props.dm[i+1], props.lnT_face[i])
            eval_mixed_p1_property_surface!(props.eos_res[i].∇ₐ, props.eos_res[i+1].∇ₐ, props.dm[i], props.dm[i+1], props.∇ₐ_face[i])
            eval_mixed_p1_property_surface!(props.eos_res[i].δ, props.eos_res[i+1].δ, props.dm[i], props.dm[i+1], props.δ_face[i])
            eval_mixed_p1_property_surface!(props.eos_res[i].cₚ, props.eos_res[i+1].cₚ, props.dm[i], props.dm[i+1], props.cₚ_face[i])
        else 
            eval_mixed_p1_property_log!(props.κ[i], props.κ[i + 1], props.dm[i], props.dm[i+1], props.κ_face[i])
            eval_mixed_p1_property!(props.eos_res[i].lnP, props.eos_res[i+1].lnP, props.dm[i], props.dm[i+1], props.lnP_face[i])
            eval_mixed_p1_property!(props.eos_res[i].lnρ, props.eos_res[i+1].lnρ, props.dm[i], props.dm[i+1], props.lnρ_face[i])
            eval_mixed_p1_property!(props.eos_res[i].lnT, props.eos_res[i+1].lnT, props.dm[i], props.dm[i+1], props.lnT_face[i])
            eval_mixed_p1_property!(props.eos_res[i].∇ₐ, props.eos_res[i+1].∇ₐ, props.dm[i], props.dm[i+1], props.∇ₐ_face[i])
            eval_mixed_p1_property!(props.eos_res[i].δ, props.eos_res[i+1].δ, props.dm[i], props.dm[i+1], props.δ_face[i])
            eval_mixed_p1_property!(props.eos_res[i].cₚ, props.eos_res[i+1].cₚ, props.dm[i], props.dm[i+1], props.cₚ_face[i])
        end 
        gamma_00_dual   = abs(get_00_dual(props.gamma_turb[i]))
        gamma_p1_dual = abs(get_p1_dual(props.gamma_turb[i+1]))
        gamma_turb_dual = (props.dm[i+1] * gamma_00_dual + props.dm[i] * gamma_p1_dual) / (props.dm[i] + props.dm[i+1])
        update_dual_data_mixed_p1!(props.gamma_turb_face[i], gamma_turb_dual)
       
        κ_face_dual = get_00_dual(props.κ_face[i])
        ρ_face_dual = exp(get_00_dual(props.lnρ_face[i]))
        T_face_dual = exp(get_00_dual(props.lnT_face[i]))
        P_face_dual = exp(get_00_dual(props.lnP_face[i]))
        ∇ₐ_face_dual = get_00_dual(props.∇ₐ_face[i])
        δ_face_dual = get_00_dual(props.δ_face[i])
        cₚ_face_dual = get_00_dual(props.cₚ_face[i])

        L₀_dual = get_00_dual(props.L[i]) * LSUN
        r_dual = exp(get_00_dual(props.lnr[i]))

        set_turb_results!(sm.turbulence, props.turb_res_dual[i],
                    κ_face_dual, L₀_dual, ρ_face_dual, P_face_dual, T_face_dual, r_dual,
                    δ_face_dual, cₚ_face_dual, ∇ₐ_face_dual, props.m[i])
        update_struct_dual_data_mixed_p1(props.turb_res[i], props.turb_res_dual[i])

        D_turb_dual = (1/3) * sqrt(2) *  abs(sinh(get_00_dual(props.gamma_turb[i]))) * 1 / (1/(P_face_dual / (ρ_face_dual * CGRAV * props.m[i]/ r_dual^2)) + 1/r_dual)
        #D_turb_dual = (1/3) * sqrt(2 * exp(gamma_turb_dual)) *  1 / (1/(P_face_dual / (ρ_face_dual * CGRAV * sm.props.m[i]/ r_dual^2)) + 1/r_dual)
        update_dual_data_mixed_p1!(props.D_turb[i], D_turb_dual)
        D_turb_fixed = get_value(sm.start_step_props.D_turb[i])
        flux_term_dual = (4π*r_dual^2*ρ_face_dual)^2*D_turb_fixed/ 
                           (0.5*(props.dm[i]+props.dm[i+1]))
        update_dual_data_mixed_p1!(props.flux_term[i], flux_term_dual)

        if get_value(props.turb_res[i].∇) < get_value(props.turb_res[i].∇ᵣ)
            props.mixing_type[i] = :convection
        else
            props.mixing_type[i] = :no_mixing
        end
    end 
    Threads.@threads for i in 1:props.nz-1
        ## P1 values ##
        if i != props.nz-1
            P_cc_p1 = get_p1_dual(props.eos_res[i+1].P)
            ρ_cc_p1 = get_p1_dual(props.eos_res[i+1].ρ)
        else
            P_face_p1 = get_p1_dual(props.eos_res[i+1].P)
            ρ_face_p1 = get_p1_dual(props.eos_res[i+1].ρ)
            P_cc_00 = get_00_dual(props.eos_res[i].P)
            ρ_cc_00 = get_00_dual(props.eos_res[i].ρ)

            P_cc_p1 = P_cc_00 + (P_face_p1 - P_cc_00)*(0.5*props.dm[i]+0.5*props.dm[i+1])/(0.5*props.dm[i]+props.dm[i+1])
            ρ_cc_p1 = ρ_cc_00 + (ρ_face_p1 - ρ_cc_00)*(0.5*props.dm[i]+0.5*props.dm[i+1])/(0.5*props.dm[i]+props.dm[i+1])
        end
        m_face_p1 = 0.5*(props.m[i] + props.m[i+1])
        r_cc_p1 = 0.5*(exp(get_p1_dual(props.lnr[i+1])) + exp(get_00_dual(props.lnr[i])))

        Hₚ_cc_p1 = P_cc_p1 / (ρ_cc_p1 * m_face_p1 * CGRAV / r_cc_p1^2)
        Λ_cc_p1 = 1 / (1 / (α_Λ_FACTOR[]*Hₚ_cc_p1) + 1 / r_cc_p1)
        
    
        γ_cc_p1 = 0.5 * (abs(sinh(get_p1_dual(props.gamma_turb[i+1]))) + abs(sinh(get_00_dual(props.gamma_turb[i]))))
        A_p1 = (4π *ρ_cc_p1* r_cc_p1^2)^2 * Λ_cc_p1 * α_w_FACTOR[] * ( γ_cc_p1)  #sqrt(ω_cc_p1)    
        
        turb_flux_term = (A_p1 / props.dm[i+1]) * ((sinh(get_p1_dual(props.gamma_turb[i+1])))^2 - (sinh(get_00_dual(props.gamma_turb[i])))^2)
        update_dual_data_mixed_p1!(props.turb_flux_cc[i], turb_flux_term)
    end
end


function Jems.StellarModels.hydro_vars(equation_set::TDCEquationSet)
    return [:lnρ, :lnT, :lnr, :lum, :gamma_turb]
end 

function Jems.StellarModels.hydro_vars_scaling(equation_set::TDCEquationSet)
    return [:log, :log, :log, :maxval, :log]
end
# function split_omega(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew_low, varnew_up)
#     # use same omega in both cells to preserve energy
#      varnew_low[sm.vari[:gamma_turb]] = var_00[sm.vari[:gamma_turb]]
#      varnew_up[sm.vari[:gamma_turb]] = var_00[sm.vari[:gamma_turb]]
# end

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




function Jems.StellarModels.build_properties_for_equation_set(equation_set::TDCEquationSet, nvars, nz, nextra, network, vari, number_type, internal_dual_tag)
    TDCStellarModelProperties(nvars, nz, nextra,
                                   length(network.reactions), network.nspecies, vari, number_type, internal_dual_tag)
end


"""
Gamma Turb equation (using the radiative term for the excess source term)
"""
const α_w_FACTOR = Ref{Float64}()
const α_Λ_FACTOR = Ref{Float64}()
const EXCESS_FACTOR = Ref{Float64}()
function gammaTurb(sm::StellarModel, k::Int)
###Constants###
C_d = 8/3 * sqrt(2/3)
# α_w = 0.25 


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
    A_p1 = (4π  *ρ_cc_p1* r_cc_p1^2)^2 * Λ_cc_p1 * α_w_FACTOR[] * sqrt(ω_cc_p1)

    F_p1 = (A_p1 / dm_cell_p1) * (exp(get_p1_dual(sm.props.gamma_turb[k+1])) - exp(get_00_dual(sm.props.gamma_turb[k]))) 

    # Different terms for residual at k = 1
    mixing_term =  (F_p1/  dm_face_p1)  
    omega_var_term = ((γ_face_00- get_value(sm.start_step_props.gamma_turb[k])) / sm.props.dt) * ω_face_00
    source_term = α₁_face_00 * SA_face_00 *sqrt(ω_face_00)
    turb_dissipation_term = C_d * (ω_face_00)^(3/2) / Λ_face_00
    rad_dissipation_term = ω_face_00 / τᵣ_face_00
    excess_term = (c_s_face_00)^2 * EXCESS_FACTOR[] / τᵣ_face_00

    return  omega_var_term  + turb_dissipation_term + rad_dissipation_term - excess_term -  source_term - mixing_term
end



## Outer boundary condition (k = sm.props.nz)
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
    excess_term = (c_s_face_00)^2 * EXCESS_FACTOR[] / τᵣ_face_00
    return  omega_var_term  + turb_dissipation_term + rad_dissipation_term - excess_term - source_term - mixing_term
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
    #excess_term = (c_s_face_00)^2 * EXCESS_FACTOR[]/ τᵣ_face_00
    excess_term = C_d * (c_s_face_00 * EXCESS_FACTOR[])^3 / Λ_face_00

    return omega_var_term  + turb_dissipation_term + rad_dissipation_term - excess_term - source_term - mixing_term
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

function equationTDC_temp_arcsin(sm::StellarModel, k::Int)
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
    γ₀ = abs(get_00_dual(sm.props.gamma_turb[k]))
    ω = (sinh(γ₀))^2  # CHANGED: exp -> sin for arcsin mapping
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
    α₂ = ρ_face*cₚ*0.5*sqrt(2/3)*Λ* abs(sinh(γ₀))
    ∇ᵣ = 3 * κ * L * P_face / (16π * CRAD * CLIGHT * CGRAV * m₀ * T_face^4)
    SA = (∇ᵣ - ∇ₐ)*(1 + α₂/k_rad)^(-1)
    ∇ = ∇ₐ + SA 
       
    dm = 0.5*(sm.props.dm[k + 1] + sm.props.dm[k])
    
    return ((lnT₊ - lnT₀) / dm + CGRAV * sm.props.m[k] / (4π * r₀^4 * Pface) * ∇) /
        (CGRAV * sm.props.m[k] / (4π * r₀^4 * Pface))

end


# arcsinh implementation of the omega variable 
# function gammaTurb_arcsin(sm::StellarModel, k::Int)
#     ###Constants###
#     C_d = 8/3 * sqrt(2/3)
#     # α_w = 0.25 

#     """
#     mixing term = 1/ dm(i,face) [ A(i+1) * (ω(i+1)- ω(i))/ dm(i+1,cell) - A(i) * (ω(i)- ω(i-1))/ dm(i,cell)]
#     A(i) = (4πr^2(c,i))^2 * ρ(i,c) * Λ(i,c)^2 * √ω(i,c) { X(i,c) = 0.5(X(i,f) + X(i-1,f))}
#     A(1) = A(nz+1)= 0
#     """

#     ## Outer boundary condition (k = 1)
#     if k == 1
#         γ_face_00 = get_00_dual(sm.props.gamma_turb[k]) #defined at the outer face 
#         ω_face_00 = (sinh(γ_face_00))^2
#         dm_cell_00 = sm.props.dm[k]  
#         m_cell_00 = sm.props.m[k]
#         r_face_00 = exp(get_00_dual(sm.props.lnr[k])) #defined at the outer face 
#         # PABLO: opacity is defined at inner face for first cell 
#         #κ_face_00  = get_00_dual(sm.props.κ[k])  #defined at outer face 
#         κ_inner_face_00  = get_00_dual(sm.props.κ[k])  #defined at outer face 
#         κ_cc_p1  = get_p1_dual(sm.props.κ[k+1])
#         L_face_00  = get_00_dual(sm.props.L[k]) * LSUN #defined at outer face 

#         P_inner_face_00 = get_00_dual(sm.props.eos_res[k].P) #defined at the inner face 
#         P_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].P)
#         ρ_inner_face_00 = get_00_dual(sm.props.eos_res[k].ρ) #defined at the inner face 
#         ρ_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].ρ)
#         T_inner_face_00 = get_00_dual(sm.props.eos_res[k].T) #defined at the inner face 
#         T_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].T) 
#         ∇ₐ_inner_face_00 = get_00_dual(sm.props.eos_res[k].∇ₐ) #defined at the inner face 
#         ∇ₐ_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].∇ₐ)
#         cₚ_inner_face_00 = get_00_dual(sm.props.eos_res[k].cₚ) #defined at the inner face 
#         cₚ_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].cₚ)

#         ## Interpolating the values to get the outer face values ##
#         ## Interpolating the values to get the outer face values ##
#         # PABLO: masses are wrong here, should be dm
#         #κ_face_00 = κ_inner_face_00 + (κ_cc_p1 - κ_inner_face_00)*(m_cell_00/(sm.props.m[k] + 0.5* sm.props.m[k+1]))
#         #P_face_00 = P_inner_face_00 + (P_cc_p1 - P_inner_face_00)*(m_cell_00/(sm.props.m[k] + 0.5* sm.props.m[k+1]))
#         #ρ_face_00 = ρ_inner_face_00 + (ρ_cc_p1 - ρ_inner_face_00)*(m_cell_00/(sm.props.m[k] + 0.5* sm.props.m[k+1]))
#         #T_face_00 = T_inner_face_00 + (T_cc_p1 - T_inner_face_00)*(m_cell_00/(sm.props.m[k] + 0.5* sm.props.m[k+1]))
#         #cₚ_face_00 = cₚ_inner_face_00 + (cₚ_cc_p1 - cₚ_inner_face_00)*(m_cell_00/(sm.props.m[k] + 0.5* sm.props.m[k+1]))
#         #∇ₐ_face_00 = ∇ₐ_inner_face_00 + (∇ₐ_cc_p1 - ∇ₐ_inner_face_00)*(m_cell_00/(sm.props.m[k] + 0.5* sm.props.m[k+1]))
#         κ_face_00 = κ_inner_face_00 + (κ_cc_p1 - κ_inner_face_00)*(m_cell_00/(sm.props.dm[k] + 0.5* sm.props.dm[k+1]))
#         P_face_00 = P_inner_face_00 + (P_cc_p1 - P_inner_face_00)*(m_cell_00/(sm.props.dm[k] + 0.5* sm.props.dm[k+1]))
#         ρ_face_00 = ρ_inner_face_00 + (ρ_cc_p1 - ρ_inner_face_00)*(m_cell_00/(sm.props.dm[k] + 0.5* sm.props.dm[k+1]))
#         T_face_00 = T_inner_face_00 + (T_cc_p1 - T_inner_face_00)*(m_cell_00/(sm.props.dm[k] + 0.5* sm.props.dm[k+1]))
#         cₚ_face_00 = cₚ_inner_face_00 + (cₚ_cc_p1 - cₚ_inner_face_00)*(m_cell_00/(sm.props.dm[k] + 0.5* sm.props.dm[k+1]))
#         ∇ₐ_face_00 = ∇ₐ_inner_face_00 + (∇ₐ_cc_p1 - ∇ₐ_inner_face_00)*(m_cell_00/(sm.props.dm[k] + 0.5* sm.props.dm[k+1]))
        
#         Hₚ_face_00 = P_face_00 / (ρ_face_00 * m_cell_00 * CGRAV / r_face_00^2)
#         Λ_face_00 = 1 / (1 / Hₚ_face_00 + 1 / r_face_00)
#         ∇ᵣ_face_00 = (3 * κ_face_00 * L_face_00 * P_face_00) / (16π * CRAD * CLIGHT * CGRAV * m_cell_00 * T_face_00^4)
#         τᵣ_face_00 = (cₚ_face_00 * κ_face_00 * ρ_face_00^2 * Λ_face_00^2) / (48 * SIGMA_SB * T_face_00^3)
#         c_s_face_00 = sqrt(P_face_00 / ρ_face_00)
#         k_rad_face_00 = (16 * SIGMA_SB * T_face_00^3) / (3 * κ_face_00 * ρ_face_00)  ##
#         α₂_face_00 = ρ_face_00 * cₚ_face_00 * 0.5 * sqrt(2 / 3) * Λ_face_00 * abs(sinh(γ_face_00))
#         α₁_face_00 = ∇ₐ_face_00 * T_face_00 * Λ_face_00 * 0.5 * sqrt(2 / 3) * cₚ_face_00 / Hₚ_face_00^2
#         SA_face_00 = (∇ᵣ_face_00 - ∇ₐ_face_00) * (1 + α₂_face_00 / k_rad_face_00)^(-1)
        
#         # CHANGED: exp -> sin
#         # dgammadt_face_00 = (ω_face_00- (get_value(sm.start_step_props.gamma_turb[k]))) / sm.props.dt * cosh(γ_face_00)
        
#         """
#         mixing term = 1/ dm(i,face) [ A(2) * (ω(2)- ω(1))/ dm(2,cell) ]

#         """

#         # ### p1 terms needed as F_1 = 0 and F_2 = 0 is the next cell so everything is wrt to i+1 so the p1 cell, also mixing term calcualted at cell centre ###
#         # P_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].P)
#         # ρ_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].ρ)
#         # T_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].T)
#         # dm_face_p1 = 0.5*(sm.props.dm[k] + sm.props.dm[k+1])##
#         # dm_cell_p1 = sm.props.dm[k+1] 
#         # m_face_p1 = 0.5*(sm.props.m[k] + sm.props.m[k+1])
#         # L_cc_p1 = 0.5*(get_p1_dual(sm.props.L[k+1]) + get_00_dual(sm.props.L[k])) * LSUN
#         # r_cc_p1 = 0.5*(exp(get_p1_dual(sm.props.lnr[k+1])) + exp(get_00_dual(sm.props.lnr[k])))
#         # Hₚ_cc_p1 = P_cc_p1 / (ρ_cc_p1 * m_face_p1 * CGRAV / r_cc_p1^2)
        
#         # # CHANGED: exp -> sin
#         # ω_cc_p1 = 0.5*((sinh(get_p1_dual(sm.props.gamma_turb[k+1])))^2 + (sinh(get_00_dual(sm.props.gamma_turb[k])))^2)
#         # μ_cc_p1 = 0.5 * (abs(sinh(get_p1_dual(sm.props.gamma_turb[k+1]))) + abs(sinh(get_00_dual(sm.props.gamma_turb[k]))))
#         # Λ_cc_p1 = 1 / (1 / Hₚ_cc_p1 + 1 / r_cc_p1)
#         # A_p1 = (4π  *ρ_cc_p1* r_cc_p1^2)^2 * Λ_cc_p1 * α_w_FACTOR[] * (μ_cc_p1)    #sqrt(ω_cc_p1)             

#         # # CHANGED: exp -> sin
#         # F_p1 = (A_p1 / dm_cell_p1) * ((sinh(get_p1_dual(sm.props.gamma_turb[k+1])))^2 - (sinh(get_00_dual(sm.props.gamma_turb[k])))^2) 

#         # # PABLO: this should be full mass of cell 1 + half mass of cell 2
#         # #mixing_term =  (F_p1/  dm_face_p1)  
#         # mixing_term =  (F_p1/  (sm.props.dm[k]+0.5*sm.props.dm[k+1]))  

#         F_00 = get_00_dual(sm.props.turb_flux_cc[k])
#         mixing_term = (F_00/  (sm.props.dm[k]+0.5*sm.props.dm[k+1]))  
        
#         # CHANGED: Chain rule application. Multiply by cos(γ) instead of ω.
#         omega_var_term = ((γ_face_00- get_value(sm.start_step_props.gamma_turb[k])) / sm.props.dt) * 2*sinh(γ_face_00)*cosh(γ_face_00) 
        
#         source_term = α₁_face_00 * SA_face_00 * abs(sinh(γ_face_00))
#         turb_dissipation_term = C_d * (abs(sinh(γ_face_00)))^3 / Λ_face_00
#         rad_dissipation_term = ω_face_00 / τᵣ_face_00
#         excess_term = C_d * (c_s_face_00 * EXCESS_FACTOR[])^3 / Λ_face_00
#         return  omega_var_term  + turb_dissipation_term + rad_dissipation_term  -  source_term - mixing_term - excess_term
#     end



#     ## Outer boundary condition (k = sm.props.nz)
#     if k == sm.props.nz

#         # ==============================================================================
#         # THERMODYNAMICS & GEOMETRY (From EOS Results = Face Values)
#         # ==============================================================================
#         # As per instruction: EOS results here are defined at the face
#         γ_face_00 = get_00_dual(sm.props.gamma_turb[k])
#         ω_face_00 = (sinh(γ_face_00))^2 # CHANGED: exp -> sin
        
#         # CHANGED: exp -> sin
#         ω_cc_00 = 0.5*((sinh(get_00_dual(sm.props.gamma_turb[k])))^2 + (sinh(get_m1_dual(sm.props.gamma_turb[k-1])))^2)
#         γ_cc_00 = 0.5* ((abs(sinh(get_00_dual(sm.props.gamma_turb[k])))) + abs(sinh(get_m1_dual(sm.props.gamma_turb[k-1]))))
#         dm_cell_00 = sm.props.dm[k]
#         m_face_00  = sm.props.m[k] # Mass at the outer boundary
#         m_cc_00 = 0.5*(sm.props.m[k] + sm.props.m[k-1])
#         r_face_00  = exp(get_00_dual(sm.props.lnr[k]))
#         r_cc_00 = 0.5*(exp(get_00_dual(sm.props.lnr[k])) + exp(get_m1_dual(sm.props.lnr[k-1])))
#         L_face_00  = get_00_dual(sm.props.L[k]) * LSUN
        
#         # Thermodynamics directly from EOS (treated as Face values)
#         P_face_00  = get_00_dual(sm.props.eos_res[k].P)
#         ρ_face_00  = get_00_dual(sm.props.eos_res[k].ρ)
#         T_face_00  = get_00_dual(sm.props.eos_res[k].T)
        
#         # Opacity and gradients (Assuming these are available in props or eos_res)
#         # If kappa/nabla are in eos_res, use those. If in props arrays, fetch index k.
#         κ_face_00  = get_00_dual(sm.props.κ[k]) 
#         ∇ₐ_face_00 = get_00_dual(sm.props.eos_res[k].∇ₐ)
#         cₚ_face_00 = get_00_dual(sm.props.eos_res[k].cₚ)

#         # PABLO: why are the ones below using cc values?
#         #Hₚ_face_00 = P_face_00 / (ρ_face_00 * m_cc_00 * CGRAV / r_cc_00^2)
#         #Λ_face_00  = 1 / (1 / Hₚ_face_00 + 1 / r_cc_00)
#         Hₚ_face_00 = P_face_00 / (ρ_face_00 * m_face_00 * CGRAV / r_face_00^2)
#         Λ_face_00  = 1 / (1 / Hₚ_face_00 + 1 / r_face_00)
        
#         ∇ᵣ_face_00 = (3 * κ_face_00 * L_face_00 * P_face_00) / (16π * CRAD * CLIGHT * CGRAV * m_face_00 * T_face_00^4)
#         τᵣ_face_00 = (cₚ_face_00 * κ_face_00 * ρ_face_00^2 * Λ_face_00^2) / (48 * SIGMA_SB * T_face_00^3)
#         c_s_face_00 = sqrt(P_face_00 / ρ_face_00)
#         k_rad_face_00 = (16 * SIGMA_SB * T_face_00^3) / (3 * κ_face_00 * ρ_face_00)
        
#         α₂_face_00 = ρ_face_00 * cₚ_face_00 * 0.5 * sqrt(2 / 3) * Λ_face_00 * abs(sinh(γ_face_00))
#         α₁_face_00 = ∇ₐ_face_00 * T_face_00 * Λ_face_00 * 0.5 * sqrt(2 / 3) * cₚ_face_00 / Hₚ_face_00^2
#         SA_face_00 = (∇ᵣ_face_00 - ∇ₐ_face_00) * (1 + α₂_face_00 / k_rad_face_00)^(-1)
        
#         dgammadt_face_00 = (γ_face_00- get_value(sm.start_step_props.gamma_turb[k])) / sm.props.dt
#         # # ==============================================================================
#         # #  FLUX CALCULATION (A_00)
#         # # ==============================================================================
#         # # We use the same Face values for A_00 as they are the definitive properties at k
#         # ρ_cc_00 = 0.5 * (exp(get_m1_dual(sm.props.lnρ_face[k-1])) + exp(get_00_dual(sm.props.lnρ_face[k])))
#         # P_cc_00 = 0.5 * (exp(get_m1_dual(sm.props.lnP_face[k-1])) + exp(get_00_dual(sm.props.lnP_face[k])))
#         # Hₚ_cc_00 = P_cc_00 / (ρ_cc_00 * m_cc_00 * CGRAV / r_cc_00^2)
#         # Λ_cc_00 = 1 / (1 / Hₚ_cc_00 + 1 / r_cc_00)
#         # A_00 = (4π * ρ_cc_00 * r_cc_00^2)^2 * Λ_face_00 * α_w_FACTOR[] * (γ_cc_00) #sqrt(ω_cc_00)
        
#         # # CHANGED: exp -> sin
#         # F_00 = (A_00 / sm.props.dm[k]) * ((sinh(get_00_dual(sm.props.gamma_turb[k])))^2 - (sinh(get_m1_dual(sm.props.gamma_turb[k-1])))^2)

#         # # ==============================================================================
#         # # RESIDUAL
#         # # ==============================================================================
#         # # PABLO: we should just use half of the surface cell mass
#         # #mixing_term = -(F_00 / sm.props.dm[k]) # Flux out (F_p1) is zero at surface
#         # mixing_term = -(F_00 / (0.5*sm.props.dm[k])) # Flux out (F_p1) is zero at surface

#         F_m1 = get_m1_dual(sm.props.turb_flux_cc[k-1])
#         mixing_term = -(F_m1/  (0.5*sm.props.dm[k]))  
        
#         # CHANGED: Chain rule. Multiply by cos(γ)
#         omega_var_term = dgammadt_face_00 * 2*sinh(γ_face_00)*cosh(γ_face_00) 
        
#         source_term = α₁_face_00 * SA_face_00 * abs(sinh(γ_face_00))
#         turb_dissipation_term = C_d * (abs(sinh(γ_face_00)))^3 / Λ_face_00
#         rad_dissipation_term = ω_face_00 / τᵣ_face_00
#         excess_term = C_d * (c_s_face_00 * EXCESS_FACTOR[])^3 / Λ_face_00
#         return  omega_var_term  + turb_dissipation_term + rad_dissipation_term - source_term - mixing_term - excess_term
#     end



#     ### Other Calculations : 1 < k < nz ###
#     begin 

#         ### face Values are required for all other terms except the mixing term ###
#         γ_face_00 = get_00_dual(sm.props.gamma_turb[k])
#         ω_face_00 = (sinh(γ_face_00))^2 ### CHANGED: exp -> sin
#         dm_cell_00 = sm.props.dm[k] 
#         m_cell_00 = sm.props.m[k]
#         r_face_00 = exp(get_00_dual(sm.props.lnr[k]))
#         P_face_00 = exp(get_00_dual(sm.props.lnP_face[k]))
#         ρ_face_00 = exp(get_00_dual(sm.props.lnρ_face[k]))
#         T_face_00 = exp(get_00_dual(sm.props.lnT_face[k]))
#         κ_face_00 = get_00_dual(sm.props.κ_face[k])
#         ∇ₐ_face_00 = get_00_dual(sm.props.∇ₐ_face[k])
#         cₚ_face_00 = get_00_dual(sm.props.cₚ_face[k])
#         L_face_00 = get_00_dual(sm.props.L[k]) * LSUN

#         Hₚ_face_00 = P_face_00 / (ρ_face_00 * m_cell_00 * CGRAV / r_face_00^2)
#         Λ_face_00 = 1 / (1 / Hₚ_face_00 + 1 / r_face_00)
#         ∇ᵣ_face_00 = (3 * κ_face_00 * L_face_00 * P_face_00) / (16π * CRAD * CLIGHT * CGRAV * m_cell_00 * T_face_00^4)
#         τᵣ_face_00 = (cₚ_face_00 * κ_face_00 * ρ_face_00^2 * Λ_face_00^2) / (48 * SIGMA_SB * T_face_00^3)
#         c_s_face_00 = sqrt(P_face_00 / ρ_face_00)
#         k_rad_face_00 = (16 * SIGMA_SB * T_face_00^3) / (3 * κ_face_00 * ρ_face_00)  ##
#         α₂_face_00 = ρ_face_00 * cₚ_face_00 * 0.5 * sqrt(2 / 3) * Λ_face_00 * abs(sinh(γ_face_00))
#         α₁_face_00 = ∇ₐ_face_00 * T_face_00 * Λ_face_00 * 0.5 * sqrt(2 / 3) * cₚ_face_00 / Hₚ_face_00^2
#         SA_face_00 = (∇ᵣ_face_00 - ∇ₐ_face_00) * (1 + α₂_face_00 / k_rad_face_00)^(-1)
#         dgammadt_face_00 = (γ_face_00- get_value(sm.start_step_props.gamma_turb[k])) / sm.props.dt

#         ## cell centre values required for mixing term ##
#         ## Both 00 and p1 terms are needed for F_i and F_i+1

#         # ## 00 values ##
#         # # CHANGED: exp -> sin
#         # ω_cc_00 = 0.5*((sinh(get_00_dual(sm.props.gamma_turb[k])))^2 + (sinh(get_m1_dual(sm.props.gamma_turb[k-1])))^2)
#         # γ_cc_00 = 0.5* ((abs(sinh(get_00_dual(sm.props.gamma_turb[k])))) + abs(sinh(get_m1_dual(sm.props.gamma_turb[k-1]))))
#         # dm_cell_00 = sm.props.dm[k]
#         # dm_face_00 = 0.5*(sm.props.dm[k] + sm.props.dm[k-1]) ##
#         # m_face_00 = 0.5*(sm.props.m[k] + sm.props.m[k-1])
#         # r_cc_00 = 0.5*(exp(get_00_dual(sm.props.lnr[k])) + exp(get_m1_dual(sm.props.lnr[k-1])))
#         # P_cc_00 = get_00_dual(sm.props.eos_res[k].P)
#         # ρ_cc_00 = get_00_dual(sm.props.eos_res[k].ρ)
#         # L_cc_00 = 0.5*(get_00_dual(sm.props.L[k]) + get_m1_dual(sm.props.L[k-1])) * LSUN
#         # Hₚ_cc_00 = P_cc_00 / (ρ_cc_00 * m_face_00 * CGRAV / r_cc_00^2)
#         # Λ_cc_00 = 1 / (1 / Hₚ_cc_00 + 1 / r_cc_00)
#         # A_00 = (4π  *ρ_cc_00 * r_cc_00^2)^2 * Λ_cc_00 * α_w_FACTOR[] *  (γ_cc_00) # sqrt(ω_cc_00)
        
#         # # CHANGED: exp -> sin
#         # F_00 = (A_00 / sm.props.dm[k]) * ((sinh(get_00_dual(sm.props.gamma_turb[k])))^2 - (sinh(get_m1_dual(sm.props.gamma_turb[k-1])))^2)

#         # ## P1 values ##
#         # P_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].P)
#         # ρ_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].ρ)
#         # T_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].T)
#         # dm_face_p1 = 0.5*(sm.props.dm[k] + sm.props.dm[k+1])
#         # m_face_p1 = 0.5*(sm.props.m[k] + sm.props.m[k+1]) ##
#         # L_cc_p1 = 0.5*(get_p1_dual(sm.props.L[k+1]) + get_00_dual(sm.props.L[k])) * LSUN
#         # r_cc_p1 = 0.5*(exp(get_p1_dual(sm.props.lnr[k+1])) + exp(get_00_dual(sm.props.lnr[k])))
#         # Hₚ_cc_p1 = P_cc_p1 / (ρ_cc_p1 * m_face_p1 * CGRAV / r_cc_p1^2)
#         # Λ_cc_p1 = 1 / (1 / Hₚ_cc_p1 + 1 / r_cc_p1)
        
#         # # CHANGED: exp -> sin
#         # ω_cc_p1 = 0.5*((sinh(get_p1_dual(sm.props.gamma_turb[k+1])))^2 + (sinh(get_00_dual(sm.props.gamma_turb[k])))^2)
#         # γ_cc_p1 = 0.5 * (abs(sinh(get_p1_dual(sm.props.gamma_turb[k+1]))) + abs(sinh(get_00_dual(sm.props.gamma_turb[k]))))
#         # A_p1 = (4π *ρ_cc_p1* r_cc_p1^2)^2 * Λ_cc_p1 * α_w_FACTOR[] * ( γ_cc_p1)  #sqrt(ω_cc_p1)    
        
#         # # CHANGED: exp -> sin
#         # F_p1 = (A_p1 / sm.props.dm[k+1]) * ((sinh(get_p1_dual(sm.props.gamma_turb[k+1])))^2 - (sinh(get_00_dual(sm.props.gamma_turb[k])))^2)

#         # # Calculation of all terms for residual 
#         # mixing_term = (F_p1 - F_00) / dm_face_p1


#         F_m1 = get_m1_dual(sm.props.turb_flux_cc[k-1])
#         F_00 = get_00_dual(sm.props.turb_flux_cc[k])
#         mixing_term = ((F_00 - F_m1)/  (0.5*sm.props.dm[k-1] + 0.5*sm.props.dm[k]))  
        
#         # CHANGED: Chain rule. Multiply by cos(γ)
#         omega_var_term = dgammadt_face_00 * 2*sinh(γ_face_00)*cosh(γ_face_00) 
        
#         source_term = α₁_face_00 * SA_face_00 * abs(sinh(γ_face_00))
#         turb_dissipation_term = C_d * (abs(sinh(γ_face_00)))^3 / Λ_face_00
#         rad_dissipation_term = ω_face_00 / τᵣ_face_00
#         excess_term = C_d * (c_s_face_00 * EXCESS_FACTOR[])^3 / Λ_face_00
#         return omega_var_term  + turb_dissipation_term + rad_dissipation_term  - source_term - mixing_term - excess_term
#     end 
# end

function gammaTurb_arcsin(sm::StellarModel, k::Int)
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
        mixing_term = (F_00/  (0.5*sm.props.dm[k]+0.5*sm.props.dm[k+1]))

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
        mixing_term = ((F_00 - F_m1)/  (0.5*sm.props.dm[k] + 0.5*sm.props.dm[k+1]))
    end 

    # variables derived from primary variables at face (00 cell)
    Hₚ_face_00 = P_face_00 / (ρ_face_00 * m_cell_00 * CGRAV / r_face_00^2)
    Λ_face_00 = 1 / (1 / (α_Λ_FACTOR[]*Hₚ_face_00) + 1 / r_face_00)
    ∇ᵣ_face_00 = (3 * κ_face_00 * L_face_00 * P_face_00) / (16π * CRAD * CLIGHT * CGRAV * m_cell_00 * T_face_00^4)
    τᵣ_face_00 = (cₚ_face_00 * κ_face_00 * ρ_face_00^2 * Λ_face_00^2) / (48 * SIGMA_SB * T_face_00^3)
    c_s_face_00 = sqrt(P_face_00 / ρ_face_00)
    k_rad_face_00 = (16 * SIGMA_SB * T_face_00^3) / (3 * κ_face_00 * ρ_face_00)
    α₂_face_00 = ρ_face_00 * cₚ_face_00 * 0.5 * sqrt(2 / 3) * Λ_face_00 * abs(sinh(γ_face_00))
    α₁_face_00 = ∇ₐ_face_00 * T_face_00 * Λ_face_00 * 0.5 * sqrt(2 / 3) * cₚ_face_00 / Hₚ_face_00^2
    SA_face_00 = (∇ᵣ_face_00 - ∇ₐ_face_00) * (1 + α₂_face_00 / k_rad_face_00)^(-1)

    # Terms for calculating the residual 
    #omega_var_term = ((γ_face_00- abs(get_value(sm.start_step_props.gamma_turb[k]))) / sm.props.dt) * 2*sinh(γ_face_00)*cosh(γ_face_00)
    gamma_old = abs(get_value(sm.start_step_props.gamma_turb[k]))

    omega_var_term = (sinh(γ_face_00)^2- sinh(gamma_old)^2) / sm.props.dt
    source_term = α₁_face_00 * SA_face_00 * abs(sinh(γ_face_00))
    turb_dissipation_term = C_d * (abs(sinh(γ_face_00)))^3 / Λ_face_00
    rad_dissipation_term = ω_face_00 / τᵣ_face_00
    excess_term = C_d * (c_s_face_00 * EXCESS_FACTOR[])^3 / Λ_face_00


    return (omega_var_term  + turb_dissipation_term + rad_dissipation_term  - source_term  - excess_term - mixing_term)
end 

function equationLuminosity_arcsin(sm::StellarModel, k::Int)
    L₀ = get_00_dual(sm.props.L[k]) * LSUN
    ρ₀ = get_00_dual(sm.props.eos_res[k].ρ)
    cₚ = get_00_dual(sm.props.eos_res[k].cₚ)
    δ = get_00_dual(sm.props.eos_res[k].δ)
    T₀ = get_00_dual(sm.props.eos_res[k].T)
    P₀ = get_00_dual(sm.props.eos_res[k].P)
    dTdt = (T₀ - get_value(sm.start_step_props.eos_res[k].T)) / sm.props.dt
    dPdt = (P₀ - get_value(sm.start_step_props.eos_res[k].P)) / sm.props.dt
    L_surf = get_value(sm.start_step_props.L[sm.props.nz]) * LSUN
    #L_surf = get_00_dual(sm.start_step_props.L[sm.props.nz]) * LSUN
    ϵnuc::typeof(L₀) = 0
    for i in eachindex(sm.network.reactions)
        ϵnuc += get_00_dual(sm.props.rates[k,i])*sm.network.reactions[i].Qvalue
    end
    if k > 1
        L₋ = get_m1_dual(sm.props.L[k-1]) * LSUN
        return ((L₀ - L₋) / sm.props.dm[k] - ϵnuc + cₚ * dTdt - (δ / ρ₀) * dPdt)/(L_surf/sm.props.m[sm.props.nz])  # no neutrinos
    else
        return (L₀ / sm.props.dm[k] - ϵnuc + cₚ * dTdt - (δ / ρ₀) * dPdt)/(L_surf/sm.props.m[sm.props.nz])  # no neutrinos
    end
end

@inline function smooth_step_func(x::T, floor::Float64, ceil::Float64) where T
    if x <= floor 
        return zero(T)
    elseif x >= ceil
        return one(T)
    else
        t = (x - floor)/(ceil - floor)

        return t * t * (3.0 -2.0 * t)
    end 

end
function gammaTurb_blend(sm::StellarModel, k::Int)
    m_solar = sm.props.m[k] / MSUN
    m_total_solar = sm.props.m[sm.props.nz] / MSUN
    m_low = 0.9 * m_total_solar
    m_high = 0.92 * m_total_solar
    w = smooth_step_func(m_solar, m_low , m_high)

    # #  Pure TDC region
    # if m_solar <= m_low
    #     return gammaTurb_arcsin(sm, k)
    # end

    #Get the target values (duals stripped)
    # if k == sm.props.nz
    #     v_turb = get_value(sm.props.turb_res[k-1].v_turb)
    #     P_face = exp(get_value(sm.props.lnP_face[k-1]))
    #     ρ_face = exp(get_value(sm.props.lnρ_face[k-1]))
    # else
    #     v_turb = get_value(sm.props.turb_res[k].v_turb)
    #     P_face = exp(get_value(sm.props.lnP_face[k]))
    #     ρ_face = exp(get_value(sm.props.lnρ_face[k]))
    # end
    
    # P_face = exp(get_value(sm.props.lnP_face[k]))
    # ρ_face = exp(get_value(sm.props.lnρ_face[k]))
    # c_s = sqrt(P_face / ρ_face)

    if k == sm.props.nz
        v_turb = get_m1_dual(sm.props.turb_res[k-1].v_turb)
        P_face = exp(get_m1_dual(sm.props.lnP_face[k-1]))
        ρ_face = exp(get_m1_dual(sm.props.lnρ_face[k-1]))
    else
        v_turb = get_00_dual(sm.props.turb_res[k].v_turb)
        P_face = exp(get_00_dual(sm.props.lnP_face[k]))
        ρ_face = exp(get_00_dual(sm.props.lnρ_face[k]))
    end
    
   
    c_s = sqrt(P_face / ρ_face)
    v_floor = 1e-10 * c_s
    v_turb_floored = max(v_turb, v_floor)
    
    γ_target = ForwardDiff.value(asinh(abs((v_turb_floored) / sqrt(2.0))))
    γ_face_active = abs(get_00_dual(sm.props.gamma_turb[k]))

    V_target = (v_turb_floored^2) / 2.0
    scale = max(ForwardDiff.value(V_target) / sm.props.dt, 1.0)

    # if w == 1.0
    #     return (γ_face_active - γ_target) * scale
    # end
    
    # V_target = (v_turb_floored^2) / 2.0
    # scale = max(ForwardDiff.value(V_target) / sm.props.dt, 1.0)
    
    res_fixed = (γ_face_active - γ_target) * scale
    
    return (1.0 - w) * gammaTurb_arcsin(sm, k) + w * res_fixed
end
function equationTDC_blend_temp(sm::StellarModel, k::Int)
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

    # 1. Fetch the alternative nabla from turb_res
    ∇_turb = get_00_dual(sm.props.turb_res[k].∇)

    # Calculating TDC ∇
    L = get_00_dual(sm.props.L[k]) * LSUN
    γ₀ = abs(get_00_dual(sm.props.gamma_turb[k]))
    ω = (sinh(γ₀))^2
    ρ_face = exp(get_00_dual(sm.props.lnρ_face[k]))
    P_face = exp(get_00_dual(sm.props.lnP_face[k]))
    r₀ = exp(get_00_dual(sm.props.lnr[k]))
    T_face = exp(get_00_dual(sm.props.lnT_face[k]))
    ∇ₐ = get_00_dual(sm.props.∇ₐ_face[k])
    cₚ =  get_00_dual(sm.props.cₚ_face[k])
    κ = get_00_dual(sm.props.κ_face[k])
    m₀ = sm.props.m[k]
    Hₚ = P_face / (ρ_face * CGRAV * m₀ / r₀^2) #defined at face 
    Λ = 1/(1/(Hₚ*α_Λ_FACTOR[]) + 1/r₀) 
    k_rad = 16 * SIGMA_SB * T_face^3 / (3 * κ * ρ_face)
    α₂ = ρ_face*cₚ*0.5*sqrt(2/3)*Λ*abs(sinh(γ₀))
    ∇ᵣ = 3 * κ * L * P_face / (16π * CRAD * CLIGHT * CGRAV * m₀ * T_face^4)
    SA = (∇ᵣ - ∇ₐ)*(1 + α₂/k_rad)^(-1)
    ∇_tdc = ∇ₐ + SA 
       
   
    m_solar = m₀ / MSUN 
    m_total_solar = sm.props.m[sm.props.nz] / MSUN
    m_low = 0.99 * m_total_solar
    m_high = 0.995 * m_total_solar
    w = smooth_step_func(m_solar, m_low , m_high)
    

    ∇ = (1.0 - w) * ∇_tdc + w * ∇_turb
       
    dm = 0.5*(sm.props.dm[k + 1] + sm.props.dm[k])
    
    return ((lnT₊ - lnT₀) / dm + CGRAV * sm.props.m[k] / (4π * r₀^4 * Pface) * ∇) /
        (CGRAV * sm.props.m[k] / (4π * r₀^4 * Pface))
end

function Jems.Evolution.eval_cell_eqs!(sm::StellarModel, ::TDCEquationSet, k::Int)
    sm.solver_data.eqs_duals[k, 1] = Evolution.equationHSE(sm, k)
    sm.solver_data.eqs_duals[k, 2] = equationTDC_temp_arcsin(sm,k)
    sm.solver_data.eqs_duals[k, 3] = Evolution.equationContinuity(sm, k)
    sm.solver_data.eqs_duals[k, 4] = equationLuminosity_arcsin(sm,k)
    sm.solver_data.eqs_duals[k, 5] = gammaTurb_arcsin(sm,k)
    # evaluate all composition equations
    for i = 1:(sm.network.nspecies)
        sm.solver_data.eqs_duals[k, sm.nvars - sm.network.nspecies + i] = Evolution.equation_composition(sm, k, sm.network.species_names[i])
    end
end
