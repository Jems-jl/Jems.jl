using ForwardDiff

abstract type AbstractStellarModelProperties end

@kwdef mutable struct StellarModelProperties{TN, TDual, TCellDualData} <: AbstractStellarModelProperties
    nz::Int
    dt::TN

    m::Vector{TN}
    dm::Vector{TN}

    eos_res_dual::Vector{EOSResults{TDual}}
    eos_res::Vector{EOSResults{TCellDualData}}

    # independent variables
    lnT::Vector{TCellDualData}
    lnρ::Vector{TCellDualData}
    lnr::Vector{TCellDualData}
    L::Vector{TCellDualData}
    xa::Matrix{TCellDualData}
    xa_dual::Matrix{TDual}

    # opacity
    κ::Vector{TCellDualData}

    #rates
    rates::Matrix{TCellDualData}
    rates_dual::Matrix{TDual}
end

function StellarModelProperties(nvars::Int, nz::Int, nextra::Int,
                                nrates::Int, nspecies::Int, vari::Dict{Symbol, Int},
                                ::Type{TN}) where {TN<:Real}

    CDDTYPE = CellDualData{nvars+1,3*nvars+1,TN}  # full dual arrays
    TD = typeof(ForwardDiff.Dual(zero(TN), (zeros(TN, nvars))...))  # only the cell duals
    
    eos_res_dual = [EOSResults{TD}() for i in 1:(nz+nextra)]
    eos_res = [EOSResults{CDDTYPE}() for i in 1:(nz+nextra)]

    lnT = [CellDualData(nvars, TN; is_ind_var=true, ind_var_i=vari[:lnT]) for i in 1:(nz+nextra)]
    lnρ = [CellDualData(nvars, TN; is_ind_var=true, ind_var_i=vari[:lnρ]) for i in 1:(nz+nextra)]
    lnr = [CellDualData(nvars, TN; is_ind_var=true, ind_var_i=vari[:lnr]) for i in 1:(nz+nextra)]
    L = [CellDualData(nvars, TN; is_ind_var=true, ind_var_i=vari[:lum]) for i in 1:(nz+nextra)]
    xa = Matrix{CDDTYPE}(undef,nz+nextra, nspecies)
    for k in 1:(nz+nextra)
        for i in 1:nspecies
            xa[k,i] = CellDualData(nvars, TN;
                        is_ind_var=true, ind_var_i=4+i)
        end
    end
    xa_dual = zeros(TD, nz+nextra, nspecies)
    rates_dual = zeros(TD, nz+nextra, nrates)
    m = zeros(TN, nz+nextra)
    dm = zeros(TN, nz+nextra)

    # for some reason using zeros just creates a bunch of instances of the same object
    # so we just initialize a vector of undef
    κ = Vector{CDDTYPE}(undef, nz+nextra)  # zeros(CDDTYPE, nz+nextra)
    for k in 1:(nz+nextra)
        κ[k] = CellDualData(nvars, TN)
    end

    rates = Matrix{CDDTYPE}(undef, nz+nextra, nrates)
    for k in 1:(nz+nextra)
        for i in 1:nrates
            rates[k,i] = CellDualData(nvars, TN)
        end
    end

    return StellarModelProperties(nz=nz,
                                  m=m,
                                  dm=dm,
                                  dt=zero(TN),
                                  eos_res_dual=eos_res_dual,
                                  eos_res=eos_res,
                                  lnT=lnT,
                                  lnρ=lnρ,
                                  lnr=lnr,
                                  L=L,
                                  xa=xa,
                                  xa_dual=xa_dual,
                                  κ=κ,
                                  rates=rates,
                                  rates_dual=rates_dual)
end

"""
    function update_stellar_model_properties!(sm, props::StellarModelProperties{TDual, TCellDualData}) where
        {TDual <: ForwardDiff.Dual, TCellDualData}

Update the stellar model properties `props` to the contents of `sm`. The goal is to save the 'state' of the StellarModel
so we can easily get properties like rates, eos, opacity values, and retrace if a retry is called.
"""
function update_stellar_model_properties!(sm, props::StellarModelProperties{TN, TDual, TCellDualData}) where
                                                {TN<:Real, TDual<:ForwardDiff.Dual, TCellDualData}
    lnT_i = sm.vari[:lnT]
    lnρ_i = sm.vari[:lnρ]
    lnr_i = sm.vari[:lnr]
    L_i = sm.vari[:lum]
    props.dt = sm.dt
    props.nz = sm.nz
    Threads.@threads for i in 1:props.nz
        # update mesh
        props.m[i] = sm.m[i]
        props.dm[i] = sm.dm[i]
        # update independent variables
        update_cell_dual_data_value!(props.lnT[i], sm.ind_vars[(i-1)*(sm.nvars)+lnT_i])
        update_cell_dual_data_value!(props.lnρ[i], sm.ind_vars[(i-1)*(sm.nvars)+lnρ_i])
        update_cell_dual_data_value!(props.lnr[i], sm.ind_vars[(i-1)*(sm.nvars)+lnr_i])
        update_cell_dual_data_value!(props.L[i], sm.ind_vars[(i-1)*(sm.nvars)+L_i])
        for j in 1:sm.network.nspecies
            update_cell_dual_data_value!(props.xa[i,j], 
                            sm.ind_vars[(i-1)*(sm.nvars)+(sm.nvars - sm.network.nspecies + j)])
            props.xa_dual[i,j] = get_cell_dual(props.xa[i,j])
        end

        lnT = get_cell_dual(props.lnT[i])
        lnρ = get_cell_dual(props.lnρ[i])
        xa = @view props.xa_dual[i,:]

        # Get EOS
        set_EOS_resultsTρ!(sm.eos, props.eos_res_dual[i], lnT, lnρ,
                            xa, sm.network.species_names)
        #names = fieldnames(EOSResults)
        #for name in names
        #    dual = getfield(props.eos_res_dual[i], name)
        #    dual_cell_data = getfield(props.eos_res[i], name)
        #    update_cell_dual_data!(dual_cell_data, dual)
        #end
        update_cell_dual_data!(props.eos_res[i].T, props.eos_res_dual[i].T)
        update_cell_dual_data!(props.eos_res[i].P, props.eos_res_dual[i].P)
        update_cell_dual_data!(props.eos_res[i].ρ, props.eos_res_dual[i].ρ)
        update_cell_dual_data!(props.eos_res[i].lnT, props.eos_res_dual[i].lnT)
        update_cell_dual_data!(props.eos_res[i].lnP, props.eos_res_dual[i].lnP)
        update_cell_dual_data!(props.eos_res[i].lnρ, props.eos_res_dual[i].lnρ)
        update_cell_dual_data!(props.eos_res[i].Prad, props.eos_res_dual[i].Prad)
        update_cell_dual_data!(props.eos_res[i].μ, props.eos_res_dual[i].μ)
        update_cell_dual_data!(props.eos_res[i].α, props.eos_res_dual[i].α)
        update_cell_dual_data!(props.eos_res[i].β, props.eos_res_dual[i].β)
        update_cell_dual_data!(props.eos_res[i].δ, props.eos_res_dual[i].δ)
        update_cell_dual_data!(props.eos_res[i].χ_ρ, props.eos_res_dual[i].χ_ρ)
        update_cell_dual_data!(props.eos_res[i].χ_T, props.eos_res_dual[i].χ_T)
        update_cell_dual_data!(props.eos_res[i].u, props.eos_res_dual[i].u)
        update_cell_dual_data!(props.eos_res[i].cₚ, props.eos_res_dual[i].cₚ)
        update_cell_dual_data!(props.eos_res[i].∇ₐ, props.eos_res_dual[i].∇ₐ)
        update_cell_dual_data!(props.eos_res[i].Γ₁, props.eos_res_dual[i].Γ₁)

        # get opacity
        κ_dual = get_opacity_resultsTρ(sm.opacity, lnT, lnρ,
                    xa, sm.network.species_names)
        update_cell_dual_data!(props.κ[i], κ_dual)

        #get rates
        rates = @view props.rates_dual[i,:]
        set_rates_for_network!(rates, sm.network, props.eos_res_dual[i], xa)
        for j in eachindex(rates)
            update_cell_dual_data!(props.rates[i,j], rates[j])
        end
    end
end