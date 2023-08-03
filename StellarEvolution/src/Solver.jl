"""
    eval_cell_eqs(sm::StellarModel, k::Int, ind_vars_view::Vector{<:TT}) where{TT<:Real}

Evaluates the stellar structure equations of the stellar model, `sm`, at cell `k`, given the view of
the independent variables, `ind_vars_view`.
""" 
function eval_cell_eqs(sm::StellarModel, k::Int, ind_vars_view::Vector{<:TT}) where{TT<:Real}
    result = Vector{TT}(undef,sm.nvars)
    # initialize as undefined, since m1 and p1 values are not defined at edges
    eosm1 = Vector{TT}(undef,7)
    eos00 = Vector{TT}(undef,7)
    eosp1 = Vector{TT}(undef,7)
    κm1::TT = NaN
    κ00::TT = NaN
    κp1::TT = NaN
    varm1::Vector{TT} = []
    var00::Vector{TT} = []
    varp1::Vector{TT} = []
    species_names = sm.varnames[sm.nvars-sm.nspecies+1:end]
    if k>1 && k<sm.nz
        varm1 = ind_vars_view[1:sm.nvars]
        var00 = ind_vars_view[sm.nvars+1:2*sm.nvars]
        varp1 = ind_vars_view[2*sm.nvars+1:3*sm.nvars]
    elseif k==1
        var00 = ind_vars_view[1:sm.nvars]
        varp1 = ind_vars_view[sm.nvars+1:2*sm.nvars]
    elseif k==sm.nz
        varm1 = ind_vars_view[1:sm.nvars]
        var00 = ind_vars_view[sm.nvars+1:2*sm.nvars]
    end

    eos00 = get_EOS_resultsTP(sm.eos, sm.isotope_data, var00[sm.vari[:lnT]],var00[sm.vari[:lnP]],
        var00[sm.nvars-sm.nspecies+1:sm.nvars],species_names)
    κ00 = get_opacity_resultsTP(sm.opacity, sm.isotope_data, var00[sm.vari[:lnT]],var00[sm.vari[:lnP]],
        var00[sm.nvars-sm.nspecies+1:sm.nvars],species_names)
    if k!=1
        eosm1 = get_EOS_resultsTP(sm.eos, sm.isotope_data, varm1[sm.vari[:lnT]],varm1[sm.vari[:lnP]],
            varm1[sm.nvars-sm.nspecies+1:sm.nvars],species_names)
        κm1 = get_opacity_resultsTP(sm.opacity, sm.isotope_data, varm1[sm.vari[:lnT]],varm1[sm.vari[:lnP]],
            varm1[sm.nvars-sm.nspecies+1:sm.nvars],species_names)
    end
    if k!=sm.nz
        eosp1 = get_EOS_resultsTP(sm.eos, sm.isotope_data, varp1[sm.vari[:lnT]],varp1[sm.vari[:lnP]],
            varp1[sm.nvars-sm.nspecies+1:sm.nvars],species_names)
        κp1 = get_opacity_resultsTP(sm.opacity, sm.isotope_data, varp1[sm.vari[:lnT]],varp1[sm.vari[:lnP]],
            varp1[sm.nvars-sm.nspecies+1:sm.nvars],species_names)
    end

    for i in 1:sm.nvars
        result[i] = sm.structure_equations[i](sm, k, varm1, var00, varp1, 
                                                     eosm1, eos00, eosp1,
                                                     κm1, κ00, κp1)
    end
    return result
end

function eval_eqs!(sm::StellarModel)
    for k in 1:sm.nz
        ki = 0
        kf = 0
        if k==1
            ki = sm.nvars*(k-1)+1
            kf = sm.nvars*(k+1)
        elseif k==sm.nz
            ki = sm.nvars*(k-2)+1
            kf = sm.nvars*(k)
        else
            ki = sm.nvars*(k-2)+1
            kf = sm.nvars*(k+1)
        end
        ind_vars_view = sm.ind_vars[ki:kf]
        sm.eqs[sm.nvars*(k-1)+1:sm.nvars*k] = eval_cell_eqs(sm, k, ind_vars_view)
    end
end

function eval_jacobian_row!(sm::StellarModel, k::Int)
    # ranges of ind_vars vector that needs to be considered, needs special cases for k=1 or nz
    # Jacobian has tridiagonal block structure:
    # x x - - - -
    # x x x - - -
    # - x x x - -
    # - - x x x -
    # - - - x x x
    # - - - - x x
    # In the first row, the first block corresponds to the derivatives of the structure equations
    # with respect to the variables at the first cell. The block to the right corresponds to the
    # derivatives with respect to the variables at the second cell.
    # TODO: complete comment
    ki = 0
    kf = 0
    if k==1
        ki = sm.nvars*(k-1)+1
        kf = sm.nvars*(k+1)
    elseif k==sm.nz
        ki = sm.nvars*(k-2)+1
        kf = sm.nvars*(k)
    else
        ki = sm.nvars*(k-2)+1
        kf = sm.nvars*(k+1)
    end
    ind_vars_view = view(sm.ind_vars,ki:kf)

    eval_eqs_wrapper = x->eval_cell_eqs(sm, k, x)
    
    jac_view = view(sm.jac,sm.nvars*(k-1)+1:sm.nvars*k,ki:kf)
    #sm.jac[sm.nvars*(k-1)+1:sm.nvars*k,ki:kf] = ForwardDiff.jacobian(eval_eqs_wrapper, ind_vars_view)
    # ForwardDiff computes things in chunks of variables. Perhaps one way to precompute everything is
    # to enforce calculations are done using a chunk equal to the total variable count (essentially, do
    # all derivatives in one go). Larger chunks are more memory intensive (go as chunk_size^2), but faster.
    #cfg = ForwardDiff.JacobianConfig(eval_eqs_wrapper, ind_vars_view, ForwardDiff.Chunk{10}());
    ForwardDiff.jacobian!(jac_view, eval_eqs_wrapper, ind_vars_view)#, cfg)
end

function eval_jacobian!(sm::StellarModel)
    Threads.@threads for k in 1:sm.nz
        eval_jacobian_row!(sm, k)
    end
end