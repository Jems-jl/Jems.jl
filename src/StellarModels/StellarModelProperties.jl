@kwdef struct StellarModelProperties{TDual, TCellDualData}
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
                                    nrates::Int, nspecies::Int, vari::Dict{Symbol, Int}, ::Type{TN}) where{TN<:Real}

    CDDTYPE = CellDualData{nvars,3*nvars,TN}
    TDSC = typeof(ForwardDiff.Dual(zero(TN), (zeros(TN, nvars))...))
    
    eos_res_dual = [EOSResults{TDSC}() for i in 1:(nz+nextra)]
    eos_res = [EOSResults{CDDTYPE}() for i in 1:(nz+nextra)]

    lnT = [CellDualData(nvars, TN;
                            is_ind_var=true, ind_var_i=vari[:lnT]) for i in 1:(nz+nextra)]
    lnρ = [CellDualData(nvars, TN;
                            is_ind_var=true, ind_var_i=vari[:lnρ]) for i in 1:(nz+nextra)]
    lnr = [CellDualData(nvars, TN;
                            is_ind_var=true, ind_var_i=vari[:lnr]) for i in 1:(nz+nextra)]
    L = [CellDualData(nvars, TN;
                            is_ind_var=true, ind_var_i=vari[:lum]) for i in 1:(nz+nextra)]
    xa = Matrix{CDDTYPE}(undef,nz+nextra, nspecies)
    for k in 1:(nz+nextra)
        for i in 1:nspecies
            xa[k,i] = CellDualData(nvars, TN;
                        is_ind_var=true, ind_var_i=4+i)
        end
    end
    xa_dual = zeros(TDSC, nz+nextra, nspecies)

    κ = zeros(CDDTYPE, nz+nextra)
    rates = zeros(CDDTYPE, nz+nextra, nrates)
    rates_dual = zeros(TDSC, nz+nextra, nrates)

    return StellarModelProperties(eos_res_dual=eos_res_dual, eos_res=eos_res,
                                  lnT=lnT, lnρ=lnρ, lnr=lnr, L=L, xa=xa, xa_dual=xa_dual,
                                  κ=κ, rates=rates, rates_dual=rates_dual)
end

function update_stellar_model_properties!(sm)
    Threads.@threads for i in 1:sm.nz
        lnT_i = sm.vari[:lnT]
        lnρ_i = sm.vari[:lnρ]
        lnr_i = sm.vari[:lnr]
        L_i = sm.vari[:lum]
        # update independent variables
        update_cell_dual_data_value!(sm.props.lnT[i], 
                                        sm.ind_vars[(i-1)*(sm.nvars)+lnT_i])
        update_cell_dual_data_value!(sm.props.lnρ[i], 
                                        sm.ind_vars[(i-1)*(sm.nvars)+lnρ_i])
        update_cell_dual_data_value!(sm.props.lnr[i], 
                                        sm.ind_vars[(i-1)*(sm.nvars)+lnr_i])
        update_cell_dual_data_value!(sm.props.L[i], 
                                        sm.ind_vars[(i-1)*(sm.nvars)+L_i])
        for j in 1:sm.network.nspecies
            update_cell_dual_data_value!(sm.props.xa[i,j], 
                            sm.ind_vars[(i-1)*(sm.nvars)+(sm.nvars - sm.network.nspecies + j)])
            sm.props.xa_dual[i,j] = get_cell_dual(sm.props.xa[i,j])
        end

        lnT = get_cell_dual(sm.props.lnT[i])
        lnρ = get_cell_dual(sm.props.lnρ[i])
        xa = @view sm.props.xa_dual[i,:]

        # Get EOS
        set_EOS_resultsTρ!(sm.eos, sm.props.eos_res_dual[i], lnT, lnρ,
                            xa, sm.network.species_names)
        #names = fieldnames(EOSResults)
        #for name in names
        #    dual = getfield(props.eos_res_dual[i], name)
        #    dual_cell_data = getfield(props.eos_res[i], name)
        #    update_cell_dual_data!(dual_cell_data, dual)
        #end
        update_cell_dual_data!(sm.props.eos_res[i].T, sm.props.eos_res_dual[i].T)
        update_cell_dual_data!(sm.props.eos_res[i].P, sm.props.eos_res_dual[i].P)
        update_cell_dual_data!(sm.props.eos_res[i].ρ, sm.props.eos_res_dual[i].ρ)
        update_cell_dual_data!(sm.props.eos_res[i].lnT, sm.props.eos_res_dual[i].lnT)
        update_cell_dual_data!(sm.props.eos_res[i].lnP, sm.props.eos_res_dual[i].lnP)
        update_cell_dual_data!(sm.props.eos_res[i].lnρ, sm.props.eos_res_dual[i].lnρ)
        update_cell_dual_data!(sm.props.eos_res[i].Prad, sm.props.eos_res_dual[i].Prad)
        update_cell_dual_data!(sm.props.eos_res[i].μ, sm.props.eos_res_dual[i].μ)
        update_cell_dual_data!(sm.props.eos_res[i].α, sm.props.eos_res_dual[i].α)
        update_cell_dual_data!(sm.props.eos_res[i].β, sm.props.eos_res_dual[i].β)
        update_cell_dual_data!(sm.props.eos_res[i].δ, sm.props.eos_res_dual[i].δ)
        update_cell_dual_data!(sm.props.eos_res[i].χ_ρ, sm.props.eos_res_dual[i].χ_ρ)
        update_cell_dual_data!(sm.props.eos_res[i].χ_T, sm.props.eos_res_dual[i].χ_T)
        update_cell_dual_data!(sm.props.eos_res[i].u, sm.props.eos_res_dual[i].u)
        update_cell_dual_data!(sm.props.eos_res[i].cₚ, sm.props.eos_res_dual[i].cₚ)
        update_cell_dual_data!(sm.props.eos_res[i].∇ₐ, sm.props.eos_res_dual[i].∇ₐ)
        update_cell_dual_data!(sm.props.eos_res[i].Γ₁, sm.props.eos_res_dual[i].Γ₁)

        # get opacity
        κ_dual = get_opacity_resultsTρ(sm.opacity, lnT, lnρ,
                    xa, sm.network.species_names)
        update_cell_dual_data!(sm.props.κ[i], κ_dual)

        #get rates
        rates = @view sm.props.rates_dual[i,:]
        set_rates_for_network!(rates, sm.network, sm.props.eos_res_dual[i], xa)
        for j in eachindex(rates)
            update_cell_dual_data!(sm.props.rates[i,j], rates[j])
        end
    end
end