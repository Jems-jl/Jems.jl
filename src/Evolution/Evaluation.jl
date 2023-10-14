"""
    eval_cell_eqs(sm::StellarModel, k::Int, ind_vars_view::Vector{<:TT}) where{TT<:Real}

Evaluates the stellar structure equations of the stellar model, `sm`, at cell `k`, given the view of the independent
variables, `ind_vars_view`.
"""
function eval_cell_eqs!(sm::StellarModel, k::Int)
    sm.var00[k, :] .= get_tmp(view(sm.diff_caches, k, :)[2], sm.eqs_duals[k, 1])
    κ00 = get_opacity_resultsTP(sm.opacity, sm.var00[k, sm.vari[:lnT]], sm.var00[k, sm.vari[:lnP]],
                                view(sm.var00, k, (sm.nvars - sm.network.nspecies + 1):(sm.nvars)), sm.network.species_names)

    set_EOS_resultsTP!(sm.eos, sm.eos_res[k, 2], sm.var00[k, sm.vari[:lnT]],
                       sm.var00[k, sm.vari[:lnP]],
                       view(sm.var00, k, (sm.nvars - sm.network.nspecies + 1):(sm.nvars)), sm.network.species_names)

    set_rates_for_network!(view(sm.rates_res, k, :), sm.network, sm.eos_res[k,2], 
                           view(sm.var00, k, (sm.nvars - sm.network.nspecies + 1):(sm.nvars)))

    if k != 1
        sm.varm1[k, :] .= get_tmp(view(sm.diff_caches, k, :)[1], sm.eqs_duals[k, 1])
        κm1 = get_opacity_resultsTP(sm.opacity, sm.varm1[k, sm.vari[:lnT]], sm.varm1[k, sm.vari[:lnP]],
                                    view(sm.varm1, k, (sm.nvars - sm.network.nspecies + 1):(sm.nvars)), sm.network.species_names)
        set_EOS_resultsTP!(sm.eos, sm.eos_res[k, 1], sm.varm1[k, sm.vari[:lnT]], sm.varm1[k, sm.vari[:lnP]],
                           view(sm.varm1, k, (sm.nvars - sm.network.nspecies + 1):(sm.nvars)), sm.network.species_names)
    else
        sm.varm1[k, :] .= (eltype(sm.varm1))(NaN)
        κm1 = κ00
    end
    if k != sm.nz
        sm.varp1[k, :] .= get_tmp(view(sm.diff_caches, k, :)[3], sm.eqs_duals[k, 1])
        κp1 = get_opacity_resultsTP(sm.opacity, sm.varp1[k, sm.vari[:lnT]], sm.varp1[k, sm.vari[:lnP]],
                                    view(sm.varp1, k, (sm.nvars - sm.network.nspecies + 1):(sm.nvars)), sm.network.species_names)
        set_EOS_resultsTP!(sm.eos, sm.eos_res[k, 3], sm.varp1[k, sm.vari[:lnT]], sm.varp1[k, sm.vari[:lnP]],
                           view(sm.varp1, k, (sm.nvars - sm.network.nspecies + 1):(sm.nvars)), sm.network.species_names)
    else
        sm.varp1[k, :] .= (eltype(sm.varp1))(NaN)
        κp1 = κ00
    end

    # evaluate all equations! (except composition)
    for i = 1:(sm.nvars - sm.network.nspecies)
        sm.eqs_duals[k, i] = sm.structure_equations[i].func(sm, k,
                                                            sm.varm1, sm.var00, sm.varp1,
                                                            sm.eos_res[k, 1], sm.eos_res[k, 2], sm.eos_res[k, 3],
                                                            sm.rates_res,
                                                            κm1, κ00, κp1)
    end
    # evaluate all composition equations
    for i = 1:sm.network.nspecies
        sm.eqs_duals[k, sm.nvars - sm.network.nspecies + i] = equation_composition(sm, k, sm.network.species_names[i],
                                                            sm.varm1, sm.var00, sm.varp1,
                                                            sm.eos_res[k, 1], sm.eos_res[k, 2], sm.eos_res[k, 3],
                                                            sm.rates_res,
                                                            κm1, κ00, κp1)
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
    init_diff_cache!(sm, k)  # set diff_caches to hold 1s where needed
    eval_cell_eqs!(sm, k)  # evaluate equations on the duals, so we have jacobian also
    # populate the jacobian with the relevant entries
    jacobian_Lk = sm.jacobian_L[k]
    jacobian_Dk = sm.jacobian_D[k]
    jacobian_Uk = sm.jacobian_U[k]
    for i = 1:sm.nvars
        # for the solver we normalize all rows of the Jacobian so they don't have crazy values
        if k==1
            for j=1:sm.nvars
                jacobian_Lk[i,j] = 0
                jacobian_Dk[i,j] = sm.eqs_duals[k, i].partials[j + sm.nvars]
                jacobian_Uk[i,j] = sm.eqs_duals[k, i].partials[j + 2*sm.nvars]
            end              
        elseif k==sm.nz      
            for j=1:sm.nvars  
                jacobian_Lk[i,j] = sm.eqs_duals[k, i].partials[j]
                jacobian_Dk[i,j] = sm.eqs_duals[k, i].partials[j + sm.nvars]
                jacobian_Uk[i,j] = 0
            end              
        else                 
            for j=1:sm.nvars  
                jacobian_Lk[i,j] = sm.eqs_duals[k, i].partials[j]
                jacobian_Dk[i,j] = sm.eqs_duals[k, i].partials[j + sm.nvars]
                jacobian_Uk[i,j] = sm.eqs_duals[k, i].partials[j + 2*sm.nvars]
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

Initializes the diff_caches to the values of the independent variables, and sets 1s in the correct spots where the
dx_i^k/dx_i^k entries lie.
"""
function init_diff_cache!(sm::StellarModel, k::Int)
    # set all partials to 0 for the moment
    sm.diff_caches[k, 1].dual_du[:] .= 0.0
    sm.diff_caches[k, 2].dual_du[:] .= 0.0
    sm.diff_caches[k, 3].dual_du[:] .= 0.0
    #= these indices are headache inducing...
    diff_caches[k, 2].dual_du has structure:
    (x_1, dx_1^k/dx_1^k-1, ..., dx_1^k/dx_n^k-1, dx_1^k/dx_1^k, ..., dx_1^k/dx_n^k, dx_1^k/dx_1^k+1, ..., dx_1^k/dx_n^k+1,  # subsize 3*nvars+1
     ...
     x_n, dx_n^k/dx_1^k-1, ..., dx_n^k/dx_n^k-1, dx_n^k/dx_1^k, ..., dx_n^k/dx_n^k, dx_n^k/dx_1^k+1, ..., dx_n^k/dx_n^k+1)
    diff_caches[k, 3].dual_du has numerators k -> k+1
    diff_caches[k, 1].dual_du has numerators k -> k-1
    =#
    for i = 1:(sm.nvars)
        # set variable values in du, dual_du and its corresponding non-zero derivative
        if k != 1
            sm.diff_caches[k, 1].du[i] = sm.ind_vars[sm.nvars * (k - 2) + i]
            sm.diff_caches[k, 1].dual_du[(i - 1) * (3 * sm.nvars + 1) + 1] = sm.ind_vars[sm.nvars * (k - 2) + i]
            sm.diff_caches[k, 1].dual_du[(i - 1) * (3 * sm.nvars + 1) + 1 + i] = 1.0  # dx^k-1_i/dx^k-1_i = 1!!
        end
        sm.diff_caches[k, 2].du[i] = sm.ind_vars[sm.nvars * (k - 1) + i]
        sm.diff_caches[k, 2].dual_du[(i - 1) * (3 * sm.nvars + 1) + 1] = sm.ind_vars[sm.nvars * (k - 1) + i]
        sm.diff_caches[k, 2].dual_du[(i - 1) * (3 * sm.nvars + 1) + 1 + sm.nvars + i] = 1.0  # dx^k_i/dx^k_i = 1!!
        if k != sm.nz
            sm.diff_caches[k, 3].du[i] = sm.ind_vars[sm.nvars * k + i]
            sm.diff_caches[k, 3].dual_du[(i - 1) * (3 * sm.nvars + 1) + 1] = sm.ind_vars[sm.nvars * k + i]
            sm.diff_caches[k, 3].dual_du[(i - 1) * (3 * sm.nvars + 1) + 1 + 2 * sm.nvars + i] = 1.0  # dx^k+1_i/dx^k+1_i = 1!!
        end
    end
end
