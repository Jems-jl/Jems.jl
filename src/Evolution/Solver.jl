"""
    eval_cell_eqs(sm::StellarModel, k::Int, ind_vars_view::Vector{<:TT}) where{TT<:Real}

Evaluates the stellar structure equations of the stellar model, `sm`, at cell `k`, given the view of the independent
variables, `ind_vars_view`.
"""
function eval_cell_eqs(sm::StellarModel, k::Int, ind_vars_view::Vector{<:TT}) where {TT<:Real}
    result = Vector{TT}(undef, sm.nvars)
    # initialize as undefined, since m1 and p1 values are not defined at edges
    eosm1 = Vector{TT}(undef, 7)
    eos00 = Vector{TT}(undef, 7)
    eosp1 = Vector{TT}(undef, 7)
    κm1::TT = NaN
    κ00::TT = NaN
    κp1::TT = NaN
    varm1::Vector{TT} = []
    var00::Vector{TT} = []
    varp1::Vector{TT} = []
    species_names = sm.varnames[(sm.nvars - sm.nspecies + 1):end]

    # collect required variables
    if k > 1 && k < sm.nz
        varm1 = ind_vars_view[1:(sm.nvars)]
        var00 = ind_vars_view[(sm.nvars + 1):(2 * sm.nvars)]
        varp1 = ind_vars_view[(2 * sm.nvars + 1):(3 * sm.nvars)]
    elseif k == 1
        var00 = ind_vars_view[1:(sm.nvars)]
        varp1 = ind_vars_view[(sm.nvars + 1):(2 * sm.nvars)]
    elseif k == sm.nz
        varm1 = ind_vars_view[1:(sm.nvars)]
        var00 = ind_vars_view[(sm.nvars + 1):(2 * sm.nvars)]
    end

    # collect eos and κ info (could be sped up by doing this before eval. the Jacobian!)
    eos00 = get_EOS_resultsTP(sm.eos, sm.isotope_data, var00[sm.vari[:lnT]], var00[sm.vari[:lnP]],
                              var00[(sm.nvars - sm.nspecies + 1):(sm.nvars)], species_names)
    κ00 = get_opacity_resultsTP(sm.opacity, sm.isotope_data, var00[sm.vari[:lnT]], var00[sm.vari[:lnP]],
                                var00[(sm.nvars - sm.nspecies + 1):(sm.nvars)], species_names)
    if k != 1
        eosm1 = get_EOS_resultsTP(sm.eos, sm.isotope_data, varm1[sm.vari[:lnT]], varm1[sm.vari[:lnP]],
                                  varm1[(sm.nvars - sm.nspecies + 1):(sm.nvars)], species_names)
        κm1 = get_opacity_resultsTP(sm.opacity, sm.isotope_data, varm1[sm.vari[:lnT]], varm1[sm.vari[:lnP]],
                                    varm1[(sm.nvars - sm.nspecies + 1):(sm.nvars)], species_names)
    end
    if k != sm.nz
        eosp1 = get_EOS_resultsTP(sm.eos, sm.isotope_data, varp1[sm.vari[:lnT]], varp1[sm.vari[:lnP]],
                                  varp1[(sm.nvars - sm.nspecies + 1):(sm.nvars)], species_names)
        κp1 = get_opacity_resultsTP(sm.opacity, sm.isotope_data, varp1[sm.vari[:lnT]], varp1[sm.vari[:lnP]],
                                    varp1[(sm.nvars - sm.nspecies + 1):(sm.nvars)], species_names)
    end

    # evaluate all equations!
    for i = 1:(sm.nvars)
        result[i] = sm.structure_equations[i](sm, k, varm1, var00, varp1, eosm1, eos00, eosp1, κm1, κ00, κp1)
    end
    return result
end

"""
    eval_eqs(sm::StellarModel)

Evaluates the stellar structure equations of the stellar model, `sm`, for all cells.
"""
function eval_eqs!(sm::StellarModel)
    for k = 1:(sm.nz)
        ki = 0
        kf = 0
        if k == 1
            ki = sm.nvars * (k - 1) + 1
            kf = sm.nvars * (k + 1)
        elseif k == sm.nz
            ki = sm.nvars * (k - 2) + 1
            kf = sm.nvars * (k)
        else
            ki = sm.nvars * (k - 2) + 1
            kf = sm.nvars * (k + 1)
        end
        ind_vars_view = sm.ind_vars[ki:kf]
        sm.eqs[(sm.nvars * (k - 1) + 1):(sm.nvars * k)] = eval_cell_eqs(sm, k, ind_vars_view)
    end
end

"""
    eval_jacobian_row!(sm::StellarModel, k::int)

Evaluates row `k` of the Jacobian matrix of the given StellarModel `sm`.
"""
function eval_jacobian_row!(sm::StellarModel, k::Int)
    # ranges of ind_vars vector that needs to be considered, needs special cases for k=1 or nz
    # Jacobian has tridiagonal block structure:
    # x x - - - -
    # x x x - - -
    # - x x x - -
    # - - x x x -
    # - - - x x x
    # - - - - x x
    # In the first row, the first block corresponds to the derivatives of the structure equations with respect to the
    # variables at the first cell. The block to the right corresponds to the derivatives with respect to the variables
    # at the second cell. In the second row until row nz-1, the first block contains the derivatives of the structure
    # equations wrt. the variables of the previous cell, the second block wrt. to the varibles in the current cell, and
    # the third wrt. the variables of the next cell. The last row only contains the derivatives wrt. the penultimate
    # cell and the last cell.
    ki = 0
    kf = 0
    if k == 1
        ki = sm.nvars * (k - 1) + 1
        kf = sm.nvars * (k + 1)
    elseif k == sm.nz
        ki = sm.nvars * (k - 2) + 1
        kf = sm.nvars * (k)
    else
        ki = sm.nvars * (k - 2) + 1
        kf = sm.nvars * (k + 1)
    end
    ind_vars_view = view(sm.ind_vars, ki:kf)

    eval_eqs_wrapper = x -> eval_cell_eqs(sm, k, x)

    jac_view = view(sm.jacobian, (sm.nvars * (k - 1) + 1):(sm.nvars * k), ki:kf)
    #sm.jacobian[sm.nvars*(k-1)+1:sm.nvars*k,ki:kf] = ForwardDiff.jacobian(eval_eqs_wrapper, ind_vars_view)
    # ForwardDiff computes things in chunks of variables. Perhaps one way to precompute everything is
    # to enforce calculations are done using a chunk equal to the total variable count (essentially, do
    # all derivatives in one go). Larger chunks are more memory intensive (go as chunk_size^2), but faster.
    #cfg = ForwardDiff.JacobianConfig(eval_eqs_wrapper, ind_vars_view, ForwardDiff.Chunk{10}());
    return ForwardDiff.jacobian!(jac_view, eval_eqs_wrapper, ind_vars_view)#, cfg)
end

"""
    eval_jacobian!(sm::StellarModel)

Evaluates the whole Jacobian matrix of the given StellarModel `sm`.
"""
function eval_jacobian!(sm::StellarModel)
    Threads.@threads for k = 1:(sm.nz)
        eval_jacobian_row!(sm, k)
    end
end
