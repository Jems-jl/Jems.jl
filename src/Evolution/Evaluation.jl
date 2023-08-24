"""
    eval_cell_eqs(sm::StellarModel, k::Int, ind_vars_view::Vector{<:TT}) where{TT<:Real}

Evaluates the stellar structure equations of the stellar model, `sm`, at cell `k`, given the view of the independent
variables, `ind_vars_view`.
"""
function eval_cell_eqs!(sm::StellarModel, k::Int)

    var00 = get_tmp(sm.diff_caches[k, :][2], sm.eqs_duals[k, 1])
    eos00 = get_EOS_resultsTP(sm.eos, var00[sm.vari[:lnT]], var00[sm.vari[:lnP]],
                              var00[(sm.nvars - sm.nspecies + 1):(sm.nvars)], sm.species_names)
    κ00 = get_opacity_resultsTP(sm.opacity, var00[sm.vari[:lnT]], var00[sm.vari[:lnP]],
                                var00[(sm.nvars - sm.nspecies + 1):(sm.nvars)], sm.species_names)
    if k != 1
        varm1 = get_tmp(sm.diff_caches[k, :][1], sm.eqs_duals[k, 1])
        eosm1 = get_EOS_resultsTP(sm.eos, varm1[sm.vari[:lnT]], varm1[sm.vari[:lnP]],
                                  varm1[(sm.nvars - sm.nspecies + 1):(sm.nvars)], sm.species_names)
        κm1 = get_opacity_resultsTP(sm.opacity, varm1[sm.vari[:lnT]], varm1[sm.vari[:lnP]],
                                    varm1[(sm.nvars - sm.nspecies + 1):(sm.nvars)], sm.species_names)
    else
        varm1 = Vector{eltype(var00)}(undef, length(var00[1]))
        fill!(varm1, eltype(var00)(NaN))
        eosm1 = Vector{eltype(eos00)}(undef, length(eos00[1]))
        fill!(eosm1, eltype(eos00)(NaN))
        κm1 = typeof(κ00)(NaN)
    end
    if k != sm.nz
        varp1 = get_tmp(sm.diff_caches[k, :][3], sm.eqs_duals[k, 1])
        eosp1 = get_EOS_resultsTP(sm.eos, varp1[sm.vari[:lnT]], varp1[sm.vari[:lnP]],
                                  varp1[(sm.nvars - sm.nspecies + 1):(sm.nvars)], sm.species_names)
        κp1 = get_opacity_resultsTP(sm.opacity, varp1[sm.vari[:lnT]], varp1[sm.vari[:lnP]],
                                    varp1[(sm.nvars - sm.nspecies + 1):(sm.nvars)], sm.species_names)
    else
        varp1 = Vector{eltype(var00)}(undef, length(var00[1]))
        fill!(varp1, eltype(var00)(NaN))
        eosp1 = Vector{eltype(eos00)}(undef, length(eos00[1]))
        fill!(eosp1, eltype(eos00)(NaN))
        κp1 = typeof(κ00)(NaN)
    end

    # evaluate all equations!
    for i = 1:(sm.nvars)
        sm.eqs_duals[k, i] = sm.structure_equations[i](sm, k, varm1, var00, varp1, eosm1, eos00, eosp1, κm1, κ00, κp1)
    end
end

"""
    eval_jacobian_eqs_row!(sm::StellarModel, k::int)

Evaluates row `k` of the Jacobian matrix of the given StellarModel `sm`.
"""
function eval_jacobian_eqs_row!(sm::StellarModel, k::Int)
    #=
    Jacobian has tridiagonal block structure:
    x x - - - -
    x x x - - -
    - x x x - -
    - - x x x -
    - - - x x x
    - - - - x x
    In the first row, the first block corresponds to the derivatives of the structure equations with respect to the
    variables at the first cell. The block to the right corresponds to the derivatives with respect to the variables
    at the second cell. In the second row until row nz-1, the first block contains the derivatives of the structure
    equations wrt. the variables of the previous cell, the second block wrt. to the varibles in the current cell, and
    the third wrt. the variables of the next cell. The last row only contains the derivatives wrt. the penultimate
    cell and the last cell.

    Because this function acts on duals `eqs_duals`, we immediately evaluate the equations also.
    =#
    set_diff_cache!(sm, k)  # set diff_caches to hold 1s where needed
    eval_cell_eqs!(sm, k)  # evaluate equations on the duals, so we have jacobian also
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
    jac_view = view(sm.jacobian, (sm.nvars * (k - 1) + 1):(sm.nvars * k), ki:kf)
    for i = 1:(sm.nvars)
        # populate the jacobian with the relevant entries
        for j = 1:(k == 1 || k == sm.nz ? 2 * sm.nvars : 3 * sm.nvars)
            if (k == 1)  # for k==1 the correct derivatives are displaced!
                jac_view[i, j] = sm.eqs_duals[k, i].partials[j + sm.nvars]
            else
                jac_view[i, j] = sm.eqs_duals[k, i].partials[j]
            end
        end
        # populate the eqs_numbers with relevant entries (will be RHS for linear solver)
        sm.eqs_numbers[(k - 1) * sm.nvars + i] = sm.eqs_duals[k, i].value
    end
end

"""
    eval_jacobian_eqs!(sm::StellarModel)

Evaluates the whole Jacobian matrix and equations of the given StellarModel `sm`.
"""
function eval_jacobian_eqs!(sm::StellarModel)
    Threads.@threads for k = 1:(sm.nz)
        eval_jacobian_eqs_row!(sm, k)
    end
end

"""
    set_diff_cache!(sm::StellarModel, k::Int)

Sets the diff_caches to the values of the independent variables
"""
function set_diff_cache!(sm::StellarModel, k::Int)
    for i = 1:(sm.nvars)
        # set variable values in du and dual_du
        if k != 1
            sm.diff_caches[k, 1].du[i] = sm.ind_vars[sm.nvars * (k - 2) + i]
            sm.diff_caches[k, 1].dual_du[(i - 1) * (3 * sm.nvars + 1) + 1] = sm.ind_vars[sm.nvars * (k - 2) + i]
        end
        sm.diff_caches[k, 2].du[i] = sm.ind_vars[sm.nvars * (k - 1) + i]
        sm.diff_caches[k, 2].dual_du[(i - 1) * (3 * sm.nvars + 1) + 1] = sm.ind_vars[sm.nvars * (k - 1) + i]
        if k != sm.nz
            sm.diff_caches[k, 3].du[i] = sm.ind_vars[sm.nvars * k + i]
            sm.diff_caches[k, 3].dual_du[(i - 1) * (3 * sm.nvars + 1) + 1] = sm.ind_vars[sm.nvars * k + i]
        end
    end
end
