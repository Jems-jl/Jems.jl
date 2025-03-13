using ForwardDiff
using Jems.Turbulence

@kwdef mutable struct StellarModelProperties{TN, TDual, TDualFace, TCellDualData, TFaceDualData} <: AbstractModelProperties
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

    # turbulence (i.e. convection, face valued)
    turb_res_dual::Vector{TurbResults{TDualFace}}
    turb_res::Vector{TurbResults{TFaceDualData}}

    # flux term for mixing equations (4πr^2ρ)^2 D / dm
    flux_term::Vector{TFaceDualData}

    ϵ_nuc::Vector{TN}

    mixing_type::Vector{Symbol}
end

function StellarModelProperties(nvars::Int, nz::Int, nextra::Int, nrates::Int, nspecies::Int, vari::Dict{Symbol,Int},
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
    xa = Matrix{CDDTYPE}(undef,nz+nextra, nspecies)
    for k in 1:(nz+nextra)
        for i in 1:nspecies
            xa[k,i] = CellDualData(nvars, TN;
                        is_ind_var=true, ind_var_i=4+i) # 4 in here is the number of non-composition variables being solved
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
    flux_term = Vector{FDDTYPE}(undef, nz+nextra)#zeros(FDDTYPE, nz+nextra)
    mixing_type::Vector{Symbol} = repeat([:no_mixing], nz+nextra)
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
        flux_term[k] = FaceDualData(nvars, TN)
    end

    rates_dual = zeros(TD, nz + nextra, nrates)
    rates = Matrix{CDDTYPE}(undef, nz + nextra, nrates)
    for k = 1:(nz + nextra)
        for i = 1:nrates
            rates[k, i] = CellDualData(nvars, TN)
        end
    end

    return StellarModelProperties(; ind_vars=ind_vars, model_number=zero(Int),
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
                                  κ=κ,
                                  rates=rates,
                                  ϵ_nuc=zeros(TN, nz + nextra),
                                  rates_dual=rates_dual,
                                  mixing_type=mixing_type)
end

@inline function eval_face_property_log!(prop_00, prop_p1, dm_00, dm_p1, face_prop)
    val00 = get_face_00_dual(prop_00)
    valp1 = get_face_p1_dual(prop_p1)
    valface_dual = exp((dm_p1 * log(val00) + dm_00 * log(valp1)) / (dm_00 + dm_p1))
    update_face_dual_data!(face_prop, valface_dual)
end

@inline function eval_face_property!(prop_00, prop_p1, dm_00, dm_p1, face_prop)
    val00 = get_face_00_dual(prop_00)
    valp1 = get_face_p1_dual(prop_p1)
    valface_dual = (dm_p1 * val00 + dm_00 * valp1) / (dm_00 + dm_p1)
    update_face_dual_data!(face_prop, valface_dual)
end

@generated function update_struct_cell_dual_data(obj::T1,obj_dual::T2) where{T1,T2}
    names = fieldnames(T1)
    lines::Vector{Expr} = [:(update_cell_dual_data!(obj.$name, obj_dual.$name)) for name in names]
    return quote
        $(lines...)
    end
end

@generated function update_struct_face_dual_data(obj::T1,obj_dual::T2) where{T1,T2}
    names = fieldnames(T1)
    lines::Vector{Expr} = [:(update_face_dual_data!(obj.$name, obj_dual.$name)) for name in names]
    return quote
        $(lines...)
    end
end

"""
    function evaluate_stellar_model_properties!(sm, props::StellarModelProperties{TDual, TCellDualData}) where
        {TDual <: ForwardDiff.Dual, TCellDualData}

Evaluates the stellar model properties `props` from the `ind_vars` array. The goal is to save the 'state' of the
StellarModel so we can easily get properties like rates, eos, opacity values, and retrace if a retry is called.
This does _not_ update the mesh/ind_vars arrays.
"""
function evaluate_stellar_model_properties!(sm, props::StellarModelProperties{TN,TDual,TCellDualData}) where
         {TN<:Real,TDual<:ForwardDiff.Dual,TCellDualData}
    lnT_i = sm.vari[:lnT]
    lnρ_i = sm.vari[:lnρ]
    lnr_i = sm.vari[:lnr]
    L_i = sm.vari[:lum]

    Threads.@threads for i = 1:(props.nz)
        # update independent variables
        update_cell_dual_data_value!(props.lnT[i], props.ind_vars[(i-1)*(sm.nvars)+lnT_i])
        update_cell_dual_data_value!(props.lnρ[i], props.ind_vars[(i-1)*(sm.nvars)+lnρ_i])
        update_cell_dual_data_value!(props.lnr[i], props.ind_vars[(i-1)*(sm.nvars)+lnr_i])
        update_cell_dual_data_value!(props.L[i], props.ind_vars[(i-1)*(sm.nvars)+L_i])
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
        eval_face_property_log!(props.κ[i], props.κ[i + 1], props.dm[i], props.dm[i+1], props.κ_face[i])
        eval_face_property!(props.eos_res[i].lnP, props.eos_res[i+1].lnP, props.dm[i], props.dm[i+1], props.lnP_face[i])
        eval_face_property!(props.eos_res[i].lnρ, props.eos_res[i+1].lnρ, props.dm[i], props.dm[i+1], props.lnρ_face[i])
        eval_face_property!(props.eos_res[i].lnT, props.eos_res[i+1].lnT, props.dm[i], props.dm[i+1], props.lnT_face[i])
        eval_face_property!(props.eos_res[i].∇ₐ, props.eos_res[i+1].∇ₐ, props.dm[i], props.dm[i+1], props.∇ₐ_face[i])
        eval_face_property!(props.eos_res[i].δ, props.eos_res[i+1].δ, props.dm[i], props.dm[i+1], props.δ_face[i])
        eval_face_property!(props.eos_res[i].cₚ, props.eos_res[i+1].cₚ, props.dm[i], props.dm[i+1], props.cₚ_face[i])

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

        flux_term_dual = (4π*r_dual^2*ρ_face_dual)^2*props.turb_res_dual[i].D_turb/
                            (0.5*(sm.props.dm[i]+sm.props.dm[i+1]))
        update_face_dual_data!(props.flux_term[i], flux_term_dual)

        if get_value(props.turb_res[i].∇) < get_value(props.turb_res[i].∇ᵣ)
            props.mixing_type[i] = :convection
        else
            props.mixing_type[i] = :no_mixing
        end
    end
end

"""
    function copy_mesh!(sm, props_in, props_out)

Copies over the mesh quantities, i.e., `nz`, `m`, `dm`, and `ind_vars` from `props_in` into `props_out`.
"""
function copy_mesh_properties!(sm, props_out, props_in)
    props_out.nz = props_in.nz
    Threads.@threads for i = 1:(props_in.nz)
        props_out.m[i] = props_in.m[i]
        props_out.dm[i] = props_in.dm[i]
        for j = 1:(sm.nvars)
            props_out.ind_vars[(i - 1) * sm.nvars + j] = props_in.ind_vars[(i - 1) * sm.nvars + j]
        end
    end
end

function copy_scalar_properties!(props_out, props_in)
    props_out.time = props_in.time
    props_out.dt = props_in.dt
    props_out.dt_next = props_in.dt
    props_out.model_number = props_in.model_number
    props_out.mstar = props_in.mstar
end