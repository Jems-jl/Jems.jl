"""
    eval_cell_eqs(sm::StellarModel, k::Int, ind_vars_view::Vector{<:TT}) where{TT<:Real}

Evaluates the stellar structure equations of the stellar model, `sm`, at cell `k`, given the view of the independent
variables, `ind_vars_view`.
"""
function eval_cell_eqs!(sm::StellarModel, k::Int)
    # evaluate all equations! (except composition)
    for i = 1:(sm.props.nvars - sm.network.nspecies)
        sm.solver_data.eqs_duals[k, i] = sm.structure_equations[i].func(sm.props, sm.opt, k)
    end
    # evaluate all composition equations
    for i = 1:sm.network.nspecies
        sm.solver_data.eqs_duals[k, sm.props.nvars - sm.network.nspecies + i] = equation_composition(sm, k, sm.network.species_names[i])
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
    eval_cell_eqs!(sm, k)  # evaluate equations on the duals, so we have jacobian also
    # populate the jacobian with the relevant entries
    jacobian_Lk = sm.solver_data.jacobian_L[k]
    jacobian_Dk = sm.solver_data.jacobian_D[k]
    jacobian_Uk = sm.solver_data.jacobian_U[k]
    eqs_duals = sm.solver_data.eqs_duals
    for i = 1:sm.nvars
        # for the solver we normalize all rows of the Jacobian so they don't have crazy values
        if k==1
            for j=1:sm.nvars
                jacobian_Lk[i,j] = 0
                jacobian_Dk[i,j] = eqs_duals[k, i].partials[j + sm.nvars]
                jacobian_Uk[i,j] = eqs_duals[k, i].partials[j + 2*sm.nvars]
            end              
        elseif k==sm.nz      
            for j=1:sm.nvars  
                jacobian_Lk[i,j] = eqs_duals[k, i].partials[j]
                jacobian_Dk[i,j] = eqs_duals[k, i].partials[j + sm.nvars]
                jacobian_Uk[i,j] = 0
            end              
        else                 
            for j=1:sm.nvars  
                jacobian_Lk[i,j] = eqs_duals[k, i].partials[j]
                jacobian_Dk[i,j] = eqs_duals[k, i].partials[j + sm.nvars]
                jacobian_Uk[i,j] = eqs_duals[k, i].partials[j + 2*sm.nvars]
            end
        end
        # populate the eqs_numbers with relevant entries (will be RHS for linear solver)
        sm.solver_data.eqs_numbers[(k - 1) * sm.nvars + i] = eqs_duals[k, i].value
    end
end

"""
    eval_jacobian_eqs!(sm::StellarModel)

Evaluates the whole Jacobian matrix and equations of the given StellarModel `sm`.
"""
function eval_jacobian_eqs!(sm::StellarModel)
    Threads.@threads for k = 1:(sm.props.nz)
        eval_jacobian_eqs_row!(sm, k)
    end
end
